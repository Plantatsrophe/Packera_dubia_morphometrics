"""
Multi-portal voucher connector for iDigBio Search API v2 and Symbiota REST endpoints
with Darwin Core normalization, GBIF-first deduplication, and missing-media backfilling.
"""

from __future__ import annotations

import json
import logging
import re
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple
import urllib.parse

from scripts.core.io_utils import default_http_get, to_text

logger = logging.getLogger("PortalHarvester")
IDIGBIO_SEARCH_URL = "https://search.idigbio.org/v2/search/records"
IDIGBIO_MEDIA_VIEW_URL = "https://search.idigbio.org/v2/view/mediarecords"
IDIGBIO_MEDIA_IMG_URL = "https://api.idigbio.org/v2/media"
_COUNTRY_TO_IDIGBIO = {"US": "united states", "USA": "united states", "CA": "canada", "CAN": "canada", "MX": "mexico", "MEX": "mexico"}
_COUNTRY_NAME_TO_ISO2 = {
    "united states": "US", "united states of america": "US", "usa": "US", "us": "US",
    "canada": "CA", "can": "CA", "ca": "CA", "mexico": "MX", "méxico": "MX", "mex": "MX", "mx": "MX",
}
_URL_KEYS = (
    "accessuri", "accessURI", "ac:accessURI", "originalUrl", "largeUrl", "identifier",
    "http://rs.tdwg.org/ac/terms/accessURI", "http://purl.org/dc/terms/identifier", "url", "image_url", "thumbnailUrl",
)
_DWC_FIELD_ALIASES: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("catalogNumber", ("catalogNumber", "catalognumber")), ("institutionCode", ("institutionCode", "institutioncode", "collectionCode", "collectioncode")),
    ("scientificName", ("scientificName", "sciname", "scientificname")), ("recordedBy", ("recordedBy", "recordedby", "collector")),
    ("recordNumber", ("recordNumber", "recordnumber")), ("identifiedBy", ("identifiedBy", "identifiedby")),
    ("typeStatus", ("typeStatus", "typestatus")), ("stateProvince", ("stateProvince", "stateprovince")),
    ("county", ("county",)), ("locality", ("locality", "verbatimLocality")), ("habitat", ("habitat",)),
    ("sex", ("sex",)), ("reproductiveCondition", ("reproductiveCondition", "reproductivecondition")),
    ("dynamicProperties", ("dynamicProperties", "dynamicproperties")), ("occurrenceRemarks", ("occurrenceRemarks", "occurrenceremarks")),
    ("fieldNotes", ("fieldNotes", "fieldnotes")),
)


def _normalize_portal_sources(portal_sources: Optional[Sequence[str] | str]) -> List[str]:
    if portal_sources is None:
        return ["gbif"]
    raw = [portal_sources] if isinstance(portal_sources, str) else list(portal_sources)
    norm, seen = [], set()
    for item in raw:
        for tok in re.split(r"[\s,;]+", str(item or "").strip().lower()):
            if tok and tok not in seen:
                seen.add(tok)
                norm.append(tok)
    return norm or ["gbif"]


def _clean_str(val: Any) -> str:
    if val is None or (isinstance(val, float) and val != val):
        return ""
    s = str(val).strip()
    return "" if s.lower() in {"none", "nan", "null"} else s


def _to_num(val: Any, as_int: bool = False) -> Optional[Any]:
    s = _clean_str(val)
    if not s:
        return None
    try:
        f = float(s)
        return None if f != f else (int(f) if as_int else f)
    except (TypeError, ValueError):
        return None


def _http_get_json(url: str, http_fn: Optional[Callable[[str], Any]] = None) -> Any:
    payload = http_fn(url) if http_fn is not None else default_http_get(url, timeout=20.0)
    if payload is None or isinstance(payload, (dict, list)):
        return payload
    txt = to_text(payload).strip()
    return json.loads(txt) if txt else None


def _strip_dwc_prefixes(mapping: Mapping[str, Any]) -> Dict[str, Any]:
    out = dict(mapping)
    for k, v in mapping.items():
        if ":" in str(k):
            out.setdefault(str(k).split(":", 1)[1], v)
    return out


def _first_nonempty(sources: Sequence[Mapping[str, Any]], *keys: str) -> Any:
    for src in sources:
        for key in keys:
            val = src.get(key)
            if val is not None and str(val).strip() != "":
                return val
    return None


def _normalize_country_code(raw_cc: Any, raw_country: Any) -> str:
    cc = _clean_str(raw_cc).upper()
    if len(cc) == 2 and cc.isalpha():
        return cc
    cname = _clean_str(raw_country).lower()
    return _COUNTRY_NAME_TO_ISO2.get(cname, cc[:2] if cc in {"USA", "CAN", "MEX"} else cc)


def _extract_urls_from_any(val: Any) -> List[str]:
    urls: List[str] = []
    if isinstance(val, str):
        return [p for p in re.split(r"[\s|;,]+", val.strip()) if p.startswith(("http://", "https://"))]
    if isinstance(val, Mapping):
        for k in _URL_KEYS:
            urls.extend(_extract_urls_from_any(val.get(k)))
    elif isinstance(val, (list, tuple)):
        for elem in val:
            urls.extend(_extract_urls_from_any(elem))
    return urls


def _extract_existing_record_image_url(rec: Mapping[str, Any]) -> Optional[str]:
    for fld in ("image_url", "media", "images", "mediarecords", "associatedMedia", "ac:accessURI", "accessuri"):
        found = _extract_urls_from_any(rec.get(fld))
        if found:
            return found[0]
    exts = rec.get("extensions")
    if isinstance(exts, Mapping):
        for ext_rows in exts.values():
            found = _extract_urls_from_any(ext_rows)
            if found:
                return found[0]
    return None


def _resolve_idigbio_media_url(
    item: Mapping[str, Any], data_norm: Mapping[str, Any], idx_terms: Mapping[str, Any], http_fn: Optional[Callable[[str], Any]]
) -> Optional[str]:
    for cand in (
        item.get("image_url"), item.get("accessuri"), item.get("accessURI"), data_norm.get("accessURI"),
        data_norm.get("accessuri"), data_norm.get("associatedMedia"), idx_terms.get("accessuri"),
        item.get("media"), item.get("mediarecords"), idx_terms.get("mediarecords"),
    ):
        found = _extract_urls_from_any(cand)
        if found:
            return found[0]
    for m_src in (item.get("mediarecords"), idx_terms.get("mediarecords"), data_norm.get("mediarecords")):
        if not isinstance(m_src, (list, tuple)):
            continue
        for elem in m_src:
            if not isinstance(elem, str) or not elem.strip() or elem.strip().startswith(("http://", "https://")):
                continue
            muuid = elem.strip()
            try:
                vj = _http_get_json(f"{IDIGBIO_MEDIA_VIEW_URL}/{urllib.parse.quote(muuid)}", http_fn)
                if isinstance(vj, Mapping):
                    v_data = _strip_dwc_prefixes(vj.get("data") or {})
                    v_idx = vj.get("indexTerms") if isinstance(vj.get("indexTerms"), Mapping) else {}
                    for cu in (vj.get("accessuri"), vj.get("accessURI"), v_data.get("accessURI"), v_data.get("accessuri"), v_idx.get("accessuri")):
                        found = _extract_urls_from_any(cu)
                        if found:
                            return found[0]
            except Exception:
                pass
            return f"{IDIGBIO_MEDIA_IMG_URL}/{urllib.parse.quote(muuid)}?size=fullsize"
    return None


def _build_dwc_record(
    sources: Sequence[Mapping[str, Any]], portal_source: str, image_url: Optional[str],
    lat_val: Optional[float], lon_val: Optional[float], fallback_key: Any = None,
) -> Dict[str, Any]:
    event_date = _clean_str(_first_nonempty(sources, "eventDate", "eventdate", "datecollected"))
    y_val, m_val, d_val = (_to_num(_first_nonempty(sources, k), as_int=True) for k in ("year", "month", "day"))
    if event_date and (y_val is None or m_val is None or d_val is None):
        m_dt = re.match(r"^(\d{4})(?:-(\d{1,2})(?:-(\d{1,2}))?)?", event_date)
        if m_dt:
            y_val = y_val if y_val is not None else int(m_dt.group(1))
            m_val = m_val if m_val is not None else (int(m_dt.group(2)) if m_dt.group(2) else None)
            d_val = d_val if d_val is not None else (int(m_dt.group(3)) if m_dt.group(3) else None)
    occ_id = _clean_str(_first_nonempty(sources, "occurrenceID", "occurrenceid", "uuid", "occid"))
    rec: Dict[str, Any] = {
        "portal_source": portal_source,
        "key": _clean_str(_first_nonempty(sources, "key", "uuid", "occid")) or occ_id or _clean_str(fallback_key),
        "occurrenceID": occ_id, "decimalLatitude": lat_val, "decimalLongitude": lon_val,
        "coordinateUncertaintyInMeters": _to_num(_first_nonempty(sources, "coordinateUncertaintyInMeters", "coordinateuncertainty")),
        "eventDate": event_date, "year": y_val, "month": m_val, "day": d_val,
        "countryCode": _normalize_country_code(_first_nonempty(sources, "countryCode", "countrycode"), _first_nonempty(sources, "country")),
        "image_url": image_url or "", "media": [{"identifier": image_url, "type": "StillImage", "format": "image/jpeg"}] if image_url else [],
    }
    for out_k, aliases in _DWC_FIELD_ALIASES:
        rec[out_k] = _clean_str(_first_nonempty(sources, *aliases))
    return rec


def _normalize_idigbio_item(item: Mapping[str, Any], http_fn: Optional[Callable[[str], Any]]) -> Dict[str, Any]:
    data_norm = _strip_dwc_prefixes(item.get("data") if isinstance(item.get("data"), Mapping) else {})
    idx_terms = item.get("indexTerms") if isinstance(item.get("indexTerms"), Mapping) else {}
    geopoint = idx_terms.get("geopoint") if isinstance(idx_terms.get("geopoint"), Mapping) else {}
    sources = (item, data_norm, idx_terms)
    lat_v, lon_v = _to_num(_first_nonempty(sources, "decimalLatitude", "decimallatitude")), _to_num(_first_nonempty(sources, "decimalLongitude", "decimallongitude"))
    return _build_dwc_record(
        sources, "iDigBio", _resolve_idigbio_media_url(item, data_norm, idx_terms, http_fn),
        lat_v if lat_v is not None else _to_num(geopoint.get("lat")), lon_v if lon_v is not None else _to_num(geopoint.get("lon")), fallback_key=item.get("uuid"),
    )


def _normalize_symbiota_item(item: Mapping[str, Any]) -> Dict[str, Any]:
    norm = _strip_dwc_prefixes(item)
    return _build_dwc_record(
        (norm,), "Symbiota", _extract_existing_record_image_url(norm),
        _to_num(_first_nonempty((norm,), "decimalLatitude", "decimallatitude", "lat")),
        _to_num(_first_nonempty((norm,), "decimalLongitude", "decimallongitude", "lon")), fallback_key=norm.get("occid"),
    )


def _record_dedup_keys(rec: Mapping[str, Any]) -> Tuple[Optional[Tuple[str, str]], Optional[str]]:
    inst = _clean_str(rec.get("institutionCode") or rec.get("collectionCode")).upper()
    cat = _clean_str(rec.get("catalogNumber")).upper()
    return ((inst, cat) if (inst and cat) else None), (_clean_str(rec.get("occurrenceID")).lower() or None)


def _build_idigbio_search_url(query: str, countries: Sequence[str], limit: int, offset: int) -> str:
    q_clean = re.sub(r"\s+", " ", query.strip())
    tokens = q_clean.split(" ") if q_clean else []
    rq: Dict[str, Any] = {"hasImage": True}
    if len(tokens) == 1:
        rq["genus"] = q_clean
    elif len(tokens) >= 2:
        rq["scientificname"] = q_clean
    mapped = [_COUNTRY_TO_IDIGBIO.get(str(c).strip().upper(), str(c).strip().lower()) for c in countries if str(c).strip()]
    if mapped:
        rq["country"] = mapped[0] if len(mapped) == 1 else mapped
    return f"{IDIGBIO_SEARCH_URL}?{urllib.parse.urlencode({'rq': json.dumps(rq, separators=(',', ':')), 'limit': int(limit), 'offset': int(offset)})}"


def _build_symbiota_search_url(base_url: str, query: str, limit: int, offset: int) -> str:
    parsed = urllib.parse.urlsplit(base_url.strip())
    q_params = urllib.parse.parse_qsl(parsed.query, keep_blank_values=True)
    q_params.extend([("scientificName", query.strip()), ("limit", str(int(limit))), ("offset", str(int(offset)))])
    return urllib.parse.urlunsplit((parsed.scheme, parsed.netloc, parsed.path, urllib.parse.urlencode(q_params), parsed.fragment))


def _fetch_portal_pages(
    taxa_queries: Sequence[str], max_records_per_taxon: int, url_builder: Callable[[str, int, int], str],
    http_fn: Optional[Callable[[str], Any]], item_normalizer: Callable[[Mapping[str, Any]], Dict[str, Any]], portal_label: str,
) -> List[Dict[str, Any]]:
    collected: List[Dict[str, Any]] = []
    cap = max(1, int(max_records_per_taxon))
    page_limit = min(500, cap)
    for raw_q in taxa_queries:
        query = _clean_str(raw_q)
        if not query:
            continue
        offset, taxon_count, seen_sigs = 0, 0, set()
        while taxon_count < cap:
            req_limit = min(page_limit, cap - taxon_count)
            try:
                payload = _http_get_json(url_builder(query, req_limit, offset), http_fn)
            except Exception as exc:
                logger.warning(f"{portal_label} query failed for '{query}' (offset={offset}): {exc}")
                break
            if isinstance(payload, list):
                raw_items, total_hint = payload, None
            elif isinstance(payload, Mapping):
                raw_items = payload.get("items") or payload.get("results") or payload.get("records") or payload.get("occurrences") or payload.get("data") or []
                total_hint = _to_num(payload.get("itemCount") or payload.get("count") or payload.get("total"), as_int=True)
            else:
                break
            if not isinstance(raw_items, list) or not raw_items:
                break
            new_in_page = 0
            for raw_item in raw_items:
                if not isinstance(raw_item, Mapping):
                    continue
                norm_rec = item_normalizer(raw_item)
                ic, oi = _record_dedup_keys(norm_rec)
                sig = (ic, oi, _clean_str(norm_rec.get("key")))
                if sig in seen_sigs:
                    continue
                seen_sigs.add(sig)
                collected.append(norm_rec)
                taxon_count += 1
                new_in_page += 1
                if taxon_count >= cap:
                    break
            offset += len(raw_items)
            if new_in_page == 0 or (total_hint is not None and offset >= total_hint) or (total_hint is None and len(raw_items) < req_limit):
                break
    return collected


def fetch_and_deduplicate_portal_records(
    taxa_queries: Sequence[str],
    gbif_records: Sequence[Mapping[str, Any]],
    portal_sources: Sequence[str] = ("gbif",),
    symbiota_api_urls: Sequence[str] = (),
    max_records_per_taxon: int = 200,
    countries: Sequence[str] = ("US", "CA"),
    idigbio_http_fn: Optional[Callable[[str], Any]] = None,
    symbiota_http_fn: Optional[Callable[[str], Any]] = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, int]]:
    """Query iDigBio/Symbiota portals, backfill missing GBIF media, and append net-new sheets."""
    norm_sources = _normalize_portal_sources(portal_sources)
    merged: List[Dict[str, Any]] = []
    by_inst_cat: Dict[Tuple[str, str], int] = {}
    by_occ_id: Dict[str, int] = {}

    def _register(idx: int, rec: Mapping[str, Any]) -> None:
        ic, oi = _record_dedup_keys(rec)
        if ic is not None and ic not in by_inst_cat: by_inst_cat[ic] = idx
        if oi is not None and oi not in by_occ_id: by_occ_id[oi] = idx

    for g_rec in gbif_records:
        rec_copy = dict(g_rec)
        rec_copy["portal_source"] = _clean_str(rec_copy.get("portal_source")) or "GBIF"
        _register(len(merged), rec_copy)
        merged.append(rec_copy)

    candidates: List[Dict[str, Any]] = []
    idigbio_fetched = symbiota_fetched = 0
    if "idigbio" in norm_sources and taxa_queries:
        idig_recs = _fetch_portal_pages(taxa_queries, max_records_per_taxon, lambda q, lim, off: _build_idigbio_search_url(q, countries, lim, off), idigbio_http_fn, lambda item: _normalize_idigbio_item(item, idigbio_http_fn), "iDigBio")
        idigbio_fetched = len(idig_recs)
        candidates.extend(idig_recs)
    eff_symbiota = [symbiota_api_urls] if isinstance(symbiota_api_urls, str) else list(symbiota_api_urls or ())
    if "symbiota" in norm_sources and taxa_queries and eff_symbiota:
        for base_url in eff_symbiota:
            if _clean_str(base_url):
                sym_recs = _fetch_portal_pages(taxa_queries, max_records_per_taxon, lambda q, lim, off, b=str(base_url): _build_symbiota_search_url(b, q, lim, off), symbiota_http_fn, _normalize_symbiota_item, "Symbiota")
                symbiota_fetched += len(sym_recs)
                candidates.extend(sym_recs)

    dedup_skipped = media_backfilled = net_new_added = 0
    for cand in candidates:
        ic, oi = _record_dedup_keys(cand)
        match_idx = by_inst_cat.get(ic) if ic is not None else None
        if match_idx is None and oi is not None:
            match_idx = by_occ_id.get(oi)
        if match_idx is not None:
            dedup_skipped += 1
            existing = merged[match_idx]
            if not _extract_existing_record_image_url(existing):
                cand_img = _extract_existing_record_image_url(cand)
                if cand_img:
                    existing["image_url"] = cand_img
                    if not existing.get("media"):
                        existing["media"] = [{"identifier": cand_img, "type": "StillImage", "format": "image/jpeg"}]
                    media_backfilled += 1
            continue
        _register(len(merged), cand)
        merged.append(cand)
        net_new_added += 1

    return merged, {
        "gbif_input": len(gbif_records), "idigbio_fetched": idigbio_fetched, "symbiota_fetched": symbiota_fetched,
        "deduplicated_skipped": dedup_skipped, "media_backfilled": media_backfilled,
        "net_new_added": net_new_added, "total_merged": len(merged),
    }

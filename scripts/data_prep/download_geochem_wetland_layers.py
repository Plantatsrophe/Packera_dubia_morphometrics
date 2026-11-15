#!/usr/bin/env python3
"""USGS DS-801 soil geochemistry, USFWS NWI wetland/seep, and USACE NWPL provisioner."""
from __future__ import annotations

import argparse, csv, http.client, io, json, math, os, re, sqlite3, sys, threading, zipfile
import urllib.parse
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path: sys.path.insert(0, str(PROJECT_ROOT))
from scripts.core.io_utils import atomic_write_csv_rows, default_http_get, to_text

_UA = "Plantmatics/1.0 (https://github.com/plantmatics; mailto:plantmatics@users.noreply.github.com)"
_DS801_INDEX_URL = "https://pubs.usgs.gov/ds/801/downloads/"
_DS801_SB_URL = "https://www.sciencebase.gov/catalog/item/514203dce4b0dc0b39f414dd?format=json"
_NWI_QUERY_URL = "https://fwspublicservices.wim.usgs.gov/wetlandsmapservice/rest/services/Wetlands/MapServer/0/query"
_NWPL_STATIC_BASE = "https://wetland-plants.sec.usace.army.mil/nwpl_static/v34"
_EASTERN_NWPL_REGIONS = {"EMP", "NCNE", "AGCP", "MW"}
_NWPL_SCORES = {"OBL": 1.00, "FACW": 0.75, "FAC": 0.50, "FACU": 0.25, "UPL": 0.00}
_NWI_TLS = threading.local()
_USDA_REGIONS = {"eastern mountains and piedmont": "EMP", "northcentral and northeast": "NCNE", "atlantic and gulf coastal plain": "AGCP", "midwest": "MW"}
_DS801_COLS = ("site_id", "latitude", "longitude", "horizon", "ca_pct", "mg_pct", "ni_ppm", "cr_ppm", "c_inorg_pct", "c_tot_pct")
_NWPL_COLS = ("species_standardized", "nwpl_indicator_code", "nwpl_wetland_score", "nwpl_region")
_ENV_GEOCHEM_DDL = (("soil_ca_mg_molar_ratio", "REAL"), ("soil_ultramafic_ni_cr_index", "REAL"), ("soil_carbonate_c_pct", "REAL"), ("is_nwi_wetland", "INTEGER"), ("nwi_wetland_system", "TEXT"), ("nwi_wetland_dist_m", "REAL"))


def _resolve_target_genus() -> str:
    try:
        from scripts.core import config as _cfg
        fn = getattr(_cfg, "resolve_target_genus", None)
        if callable(fn) and (val := str(fn()).strip()):
            return val
    except Exception:
        pass
    return (os.environ.get("PLANTMATICS_TARGET_GENUS") or os.environ.get("PACKERA_TARGET_GENUS") or "Packera").strip() or "Packera"


def _default_db_path() -> Path:
    try:
        from scripts.core.db import resolve_db_path
        return Path(resolve_db_path())
    except Exception:
        env = os.environ.get("PLANTMATICS_DB_PATH") or os.environ.get("PACKERA_DB_PATH")
        return Path(env) if env else (PROJECT_ROOT / "data" / "plantmatics.sqlite")


def _lat_lon_key(lat: float, lon: float) -> str:
    return f"{float(lat):.5f}_{float(lon):.5f}"


def _default_http_get(url: str, timeout: int = 30) -> bytes:
    return default_http_get(url, timeout=float(timeout), user_agent=_UA)


def _parse_geochem_num(val: Any) -> Optional[float]:
    s = str(val or "").strip()
    if not s or s.upper() in {"NA", "ND", "N.D.", "NONE", "NULL", "INS", "N.A."}:
        return None
    half = s.startswith("<")
    try:
        v = float(s[1:].strip() if half else s)
    except ValueError:
        return None
    if not math.isfinite(v) or v <= -999.0:
        return None
    return abs(v) * 0.5 if (half or v < 0.0) else v


def _haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    p1, p2, dp, dl = map(math.radians, (lat1, lat2, lat2 - lat1, lon2 - lon1))
    a = math.sin(dp / 2.0) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2.0) ** 2
    return 12742.0 * math.asin(min(1.0, math.sqrt(max(0.0, a))))


def compute_idw_soil_geochem(
    lat: float, lon: float, ds801_sites: Sequence[Dict[str, Any]], max_dist_km: float = 75.0, k: int = 5,
) -> Dict[str, Optional[float]]:
    """Interpolate USGS DS-801 Ca/Mg molar ratio, ultramafic Ni+Cr index, and carbonate C%."""
    empty: Dict[str, Optional[float]] = {"soil_ca_mg_molar_ratio": None, "soil_ultramafic_ni_cr_index": None, "soil_carbonate_c_pct": None}
    if not (math.isfinite(lat) and math.isfinite(lon) and ds801_sites):
        return empty
    nearby: List[Tuple[float, Dict[str, Any]]] = []
    for s in ds801_sites:
        try:
            slat = float(s.get("latitude") if s.get("latitude") is not None else s.get("lat"))
            slon = float(s.get("longitude") if s.get("longitude") is not None else s.get("lon"))
        except (TypeError, ValueError):
            continue
        if math.isfinite(slat) and math.isfinite(slon) and (d := _haversine_km(lat, lon, slat, slon)) <= max_dist_km:
            nearby.append((d, s))
    if not nearby:
        return empty
    nearest = sorted(nearby, key=lambda item: item[0])[:max(1, int(k))]

    def _idw(key: str) -> Optional[float]:
        w_sum, v_sum = 0.0, 0.0
        for d, rec in nearest:
            if (raw := rec.get(key)) in (None, ""):
                continue
            try:
                val = float(raw)
            except (TypeError, ValueError):
                continue
            if math.isfinite(val):
                w = 1.0 / (max(d, 0.25) ** 2)
                w_sum += w
                v_sum += w * val
        return (v_sum / w_sum) if w_sum > 0.0 else None

    ca, mg, ni, cr, c_inorg = _idw("ca_pct"), _idw("mg_pct"), _idw("ni_ppm"), _idw("cr_ppm"), _idw("c_inorg_pct")
    return {
        "soil_ca_mg_molar_ratio": round((ca / 40.078) / max(mg / 24.305, 1e-4), 6) if (ca is not None and mg is not None) else None,
        "soil_ultramafic_ni_cr_index": round(math.log1p(max(ni or 0.0, 0.0) + max(cr or 0.0, 0.0)), 6) if (ni is not None or cr is not None) else None,
        "soil_carbonate_c_pct": round(max(c_inorg, 0.0), 6) if c_inorg is not None else None,
    }


def _parse_ds801_table(text: str, default_horizon: str = "A") -> List[Dict[str, Any]]:
    lines = text.splitlines()
    hdr_idx = next((i for i, ln in enumerate(lines[:25]) if "latitude" in ln.lower() and "longitude" in ln.lower() and ("\t" in ln or "," in ln)), 0)
    lines = lines[hdr_idx:]
    delim = "\t" if (lines and "\t" in lines[0]) else ","
    out: List[Dict[str, Any]] = []
    for row in csv.DictReader(lines, delimiter=delim):
        raw_low = {str(k).strip().lower(): v for k, v in row.items() if k is not None}
        low = {**{re.sub(r"^(?:top5|a|c)_", "", k): v for k, v in raw_low.items()}, **raw_low}
        sid = next((str(low[k]).strip() for k in ("site_id", "siteid", "lab_id", "labid", "field_id", "id") if low.get(k)), "")
        try:
            lat = float(next((low[k] for k in ("latitude", "lat_nad83", "lat") if low.get(k) not in (None, "")), ""))
            lon = float(next((low[k] for k in ("longitude", "lon_nad83", "long_nad83", "long", "lon") if low.get(k) not in (None, "")), ""))
        except (TypeError, ValueError):
            continue
        if not (sid and math.isfinite(lat) and math.isfinite(lon)):
            continue
        hor = str(low.get("horizon") or low.get("layer") or default_horizon).strip() or default_horizon
        _col = lambda prefs: next((_parse_geochem_num(cv) for p in prefs for ck, cv in low.items() if ck == p or ck.startswith(f"{p}_")), None)
        out.append({
            "site_id": sid, "latitude": lat, "longitude": lon, "horizon": hor,
            "ca_pct": _col(("ca_pct", "ca_icp40", "ca")), "mg_pct": _col(("mg_pct", "mg_icp40", "mg")),
            "ni_ppm": _col(("ni_ppm", "ni_icp40", "ni")), "cr_ppm": _col(("cr_ppm", "cr_icp40", "cr")),
            "c_inorg_pct": _col(("c_inorg_pct", "c_inorg", "inorg_c")), "c_tot_pct": _col(("c_tot_pct", "c_tot", "tot_c", "c_total")),
        })
    return out


def _horizon_rank(horizon: str) -> int:
    h = str(horizon or "").strip().lower()
    return 0 if (h.startswith("a") or "top5" in h) else (2 if "c" in h else 1)


def ensure_usgs_ds801_soils(
    env_dir: Path, force: bool = False, http_get_fn: Optional[Callable[..., Union[bytes, str]]] = None,
) -> List[Dict[str, Any]]:
    """Download or load cached USGS DS-801 soil geochemistry into `usgs_ds801_soils.csv`."""
    csv_path = Path(env_dir) / "usgs_ds801_soils.csv"
    if csv_path.is_file() and csv_path.stat().st_size > 0 and not force:
        return _parse_ds801_table(csv_path.read_text(encoding="utf-8", errors="replace"))
    getter = http_get_fn or _default_http_get
    links: List[Tuple[str, str]] = []
    try:
        sb = json.loads(to_text(getter(_DS801_SB_URL, timeout=25)))
        links.extend((str(f.get("name", "")), str(f["url"])) for f in (sb.get("files") or []) if isinstance(f, dict) and f.get("url") and str(f.get("name", "")).lower().endswith((".csv", ".txt", ".zip")))
    except Exception:
        pass
    try:
        html = to_text(getter(_DS801_INDEX_URL, timeout=25))
        links.extend((Path(h).name, urllib.parse.urljoin(_DS801_INDEX_URL, h)) for h in re.findall(r'href=["\']([^"\']+\.(?:csv|txt|zip))["\']', html, re.I))
    except Exception:
        pass
    links.sort(key=lambda item: (0 if any(k in item[0].lower() for k in ("ahorizon", "a_horizon")) else (1 if "top5" in item[0].lower() else (3 if any(k in item[0].lower() for k in ("chorizon", "c_horizon")) else 2)), item[0]))
    by_site: Dict[str, Dict[str, Any]] = {}
    for name, url in links[:6]:
        try:
            raw = getter(url, timeout=60)
            payload = raw.encode("utf-8") if isinstance(raw, str) else raw
        except Exception: continue
        texts: List[Tuple[str, str]] = []
        if payload[:2] == b"PK":
            with zipfile.ZipFile(io.BytesIO(payload)) as zf:
                for m in sorted(zf.namelist()):
                    if m.lower().endswith((".csv", ".txt")) and "__MACOSX" not in m:
                        texts.append(("C" if "c_horizon" in m.lower() else ("Top5cm" if "top5" in m.lower() else "A"), zf.read(m).decode("utf-8", errors="replace")))
        else:
            texts.append(("C" if "c_horizon" in name.lower() else ("Top5cm" if "top5" in name.lower() else "A"), payload.decode("utf-8", errors="replace")))
        for hor, txt in texts:
            for rec in _parse_ds801_table(txt, default_horizon=hor):
                prev = by_site.get(rec["site_id"])
                if prev is None or _horizon_rank(rec["horizon"]) < _horizon_rank(prev["horizon"]): by_site[rec["site_id"]] = rec
        if by_site: break
    if not by_site:
        print("[WARN] USGS DS-801 geochemistry unavailable offline; skipping DS-801 staging.", file=sys.stderr)
        return []
    atomic_write_csv_rows(csv_path, _DS801_COLS, ({k: ("" if rec.get(k) is None else rec[k]) for k in _DS801_COLS} for rec in by_site.values()))
    return list(by_site.values())


def _point_to_rings_dist_m(lat: float, lon: float, rings: Sequence[Sequence[Sequence[float]]]) -> float:
    cos_lat, inside = max(math.cos(math.radians(lat)), 0.1), False
    for ring in rings:
        if len(ring) < 3:
            continue
        j = len(ring) - 1
        for i in range(len(ring)):
            xi, yi, xj, yj = float(ring[i][0]), float(ring[i][1]), float(ring[j][0]), float(ring[j][1])
            if ((yi > lat) != (yj > lat)) and (lon < (xj - xi) * (lat - yi) / ((yj - yi) or 1e-12) + xi):
                inside = not inside
            j = i
    if inside:
        return 0.0
    min_d = float("inf")
    for ring in rings:
        for i in range(len(ring) - 1):
            ax, ay = (float(ring[i][0]) - lon) * 111320.0 * cos_lat, (float(ring[i][1]) - lat) * 111320.0
            bx, by = (float(ring[i + 1][0]) - lon) * 111320.0 * cos_lat, (float(ring[i + 1][1]) - lat) * 111320.0
            dx, dy = bx - ax, by - ay
            denom = dx * dx + dy * dy
            t = max(0.0, min(1.0, -(ax * dx + ay * dy) / denom)) if denom > 0.0 else 0.0
            min_d = min(min_d, math.hypot(ax + t * dx, ay + t * dy))
    return min_d


def _nwi_get_json(q_str: str) -> Dict[str, Any]:
    parts = urllib.parse.urlsplit(_NWI_QUERY_URL)
    for attempt in range(2):
        qs = q_str if attempt == 0 else q_str.replace("returnGeometry=true", "returnGeometry=false")
        conn = getattr(_NWI_TLS, "conn", None)
        if conn is None:
            conn = _NWI_TLS.conn = http.client.HTTPSConnection(parts.netloc, timeout=8.0)
        try:
            conn.request("GET", f"{parts.path}?{qs}", headers={"User-Agent": _UA, "Connection": "keep-alive"})
            resp = conn.getresponse(); body = resp.read()
            if resp.status >= 400: raise OSError(f"HTTP {resp.status}")
            return json.loads(to_text(body))
        except Exception:
            conn.close(); _NWI_TLS.conn = None
            if attempt == 1: raise
    return {}


def query_nwi_point_or_cache(
    lat: float, lon: float, cache_map: Dict[str, Dict[str, Any]], cache_path: Path,
    cache_lock: Optional[threading.Lock] = None, skip_network: bool = False,
    http_get_fn: Optional[Callable[..., Union[bytes, str]]] = None,
) -> Dict[str, Any]:
    """Query USFWS NWI REST MapServer within 250 m or reuse `nwi_point_cache.jsonl`."""
    key = _lat_lon_key(lat, lon)
    if key in cache_map:
        return cache_map[key]
    if not (24.0 <= lat <= 50.0 and -125.0 <= lon <= -66.0) or skip_network:
        return {"lat_lon_key": key, "is_nwi_wetland": None, "nwi_wetland_system": None, "nwi_wetland_dist_m": None}
    dlat, dlon = 250.0 / 111320.0, 250.0 / (111320.0 * max(math.cos(math.radians(lat)), 0.1))
    params = urllib.parse.urlencode({
        "geometry": f"{lon - dlon:.6f},{lat - dlat:.6f},{lon + dlon:.6f},{lat + dlat:.6f}", "where": "Shape_Area < 50000000",
        "geometryType": "esriGeometryEnvelope", "inSR": "4326", "outSR": "4326",
        "spatialRel": "esriSpatialRelIntersects", "outFields": "ATTRIBUTE,WETLAND_TYPE", "returnGeometry": "true", "f": "json",
    })
    try:
        payload = json.loads(to_text(http_get_fn(f"{_NWI_QUERY_URL}?{params}", timeout=20))) if http_get_fn is not None else _nwi_get_json(params)
    except Exception as exc:
        print(f"[WARN] NWI REST query failed at {key}: {exc}", file=sys.stderr)
        return {"lat_lon_key": key, "is_nwi_wetland": None, "nwi_wetland_system": None, "nwi_wetland_dist_m": None}
    best_dist, best_sys = 250.0, "UPL"
    for feat in payload.get("features") or []:
        rings = (feat.get("geometry") or {}).get("rings") or []
        dist_m = _point_to_rings_dist_m(lat, lon, rings) if rings else 250.0
        if dist_m <= best_dist:
            attrs = feat.get("attributes") or {}
            attr = str(attrs.get("ATTRIBUTE") or attrs.get("WETLAND_TYPE") or "WET").strip().upper()
            m_code = re.match(r"^[A-Z]{1,3}", attr)
            best_dist, best_sys = dist_m, (m_code.group(0) if m_code else attr[:8])
    res = {
        "lat_lon_key": key, "is_nwi_wetland": 1 if (best_sys != "UPL" and best_dist <= 30.0) else 0,
        "nwi_wetland_system": best_sys, "nwi_wetland_dist_m": round(min(best_dist, 250.0), 2),
    }
    with (cache_lock if cache_lock is not None else nullcontext()):
        cache_map[key] = res
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(res) + "\n")
    return res


def _fetch_usda_nwpl_species(genus: str, by_sp: Dict[str, List[Tuple[str, float]]]) -> None:
    try:
        from scripts.core.usda_plants_ranges import load_or_fetch_plantlst_index
        sp_map = (load_or_fetch_plantlst_index(allow_network=True) or {}).get("species_map") or {}
    except Exception:
        return
    syms = [(sp, syms[0]) for sp, syms in sorted(sp_map.items()) if sp.startswith(f"{genus.lower()} ") and syms]
    def _one(item: Tuple[str, str]) -> List[Tuple[str, str, str]]:
        sp_low, sym = item
        sp_std = f"{sp_low.split()[0].capitalize()} {sp_low.split()[1].lower()}"
        try:
            prof = json.loads(to_text(_default_http_get(f"https://plantsservices.sc.egov.usda.gov/api/PlantProfile?symbol={urllib.parse.quote(sym)}", timeout=15)))
            if not (isinstance(prof, dict) and prof.get("HasWetlandData") and prof.get("Id")): return []
            wlist = json.loads(to_text(_default_http_get(f"https://plantsservices.sc.egov.usda.gov/api/PlantWetland/{prof['Id']}", timeout=15)))
            desigs = (wlist[0] if isinstance(wlist, list) and wlist else (wlist if isinstance(wlist, dict) else {})).get("WetlandDesignations") or []
            return [(sp_std, _USDA_REGIONS.get(str(d.get("Region") or "").replace("&", "and").strip().lower(), str(d.get("Region") or "EAST").strip().upper()), str(d.get("WetlandCode") or "").strip().upper().rstrip("+-*")) for d in desigs]
        except Exception:
            return []
    with ThreadPoolExecutor(max_workers=12) as pool:
        for res in pool.map(_one, syms):
            for sp_std, reg, code in res:
                if code in _NWPL_SCORES: by_sp.setdefault(sp_std, []).append((reg, _NWPL_SCORES[code]))


def fetch_nwpl_wetland_ratings(
    env_dir: Path, target_genus: str, conn: Optional[sqlite3.Connection] = None,
    force: bool = False, http_get_fn: Optional[Callable[..., Union[bytes, str]]] = None,
) -> Dict[str, Dict[str, Any]]:
    """Fetch USACE NWPL wetland indicator codes/scores and sync optional SQLite backbone profiles."""
    out_csv = Path(env_dir) / "usace_nwpl_ratings.csv"
    genus = str(target_genus or _resolve_target_genus()).strip() or "Packera"
    ratings: Dict[str, Dict[str, Any]] = {}
    if out_csv.is_file() and out_csv.stat().st_size > 0 and not force:
        with open(out_csv, encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r.get("species_standardized"):
                    sp_norm = str(r["species_standardized"]).replace("_", " ").strip()
                    ratings[sp_norm] = {**r, "species_standardized": sp_norm, "nwpl_wetland_score": float(r["nwpl_wetland_score"])}
    if not ratings:
        getter, q = http_get_fn or _default_http_get, urllib.parse.quote(genus)
        by_sp: Dict[str, List[Tuple[str, float]]] = {}
        for url in (f"{_NWPL_STATIC_BASE}/api/species?scientificName={q}", f"https://plantsservices.sc.egov.usda.gov/api/CharacteristicsSearch?scientificName={q}"):
            try:
                obj = json.loads(to_text(getter(url, timeout=25)))
            except Exception:
                continue
            for item in (obj if isinstance(obj, list) else (obj.get("data") or obj.get("species") or obj.get("records") or [])):
                if not isinstance(item, dict): continue
                toks = str(item.get("scientificName") or item.get("ScientificName") or item.get("species") or "").replace("_", " ").split()
                if len(toks) < 2 or toks[0].lower() != genus.lower(): continue
                reg = str(item.get("region") or item.get("regionCode") or item.get("nwplRegion") or "EAST").strip().upper()
                code = str(item.get("indicatorStatus") or item.get("wetlandIndicator") or item.get("rating") or "").strip().upper().rstrip("+-")
                if code in _NWPL_SCORES:
                    by_sp.setdefault(f"{toks[0].capitalize()} {toks[1].lower()}", []).append((reg or "EAST", _NWPL_SCORES[code]))
            if by_sp: break
        if not by_sp and http_get_fn is None:
            _fetch_usda_nwpl_species(genus, by_sp)
        for sp, entries in sorted(by_sp.items()):
            east = [sc for r, sc in entries if r in _EASTERN_NWPL_REGIONS]
            vals = east or [sc for _, sc in entries]
            mean_sc = round(sum(vals) / len(vals), 4)
            ratings[sp] = {"species_standardized": sp, "nwpl_indicator_code": min(_NWPL_SCORES.items(), key=lambda kv: abs(kv[1] - mean_sc))[0], "nwpl_wetland_score": mean_sc, "nwpl_region": ";".join(sorted({r for r, _ in entries if not east or r in _EASTERN_NWPL_REGIONS})) or "EAST"}
        if ratings:
            atomic_write_csv_rows(out_csv, _NWPL_COLS, ratings.values())
        else:
            print(f"[WARN] USACE NWPL wetland ratings unavailable for {genus} (continuing).", file=sys.stderr)
    if conn is not None and ratings and conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='taxon_backbone_profiles'").fetchone():
        cols = {r[1] for r in conn.execute("PRAGMA table_info(taxon_backbone_profiles)").fetchall()}
        if {"species_standardized", "nwpl_indicator_code", "nwpl_wetland_score"}.issubset(cols):
            conn.executemany("UPDATE taxon_backbone_profiles SET nwpl_indicator_code = ?, nwpl_wetland_score = ? WHERE species_standardized = ?", [(r["nwpl_indicator_code"], r["nwpl_wetland_score"], r["species_standardized"]) for r in ratings.values()])
            conn.commit()
    return ratings


def _collect_points(vouchers_csv: Optional[Path], conn: Optional[sqlite3.Connection]) -> Dict[str, Tuple[float, float]]:
    pts: Dict[str, Tuple[float, float]] = {}
    def _add(lat_v: Any, lon_v: Any) -> None:
        try:
            lat, lon = float(lat_v), float(lon_v)
            if math.isfinite(lat) and math.isfinite(lon): pts[_lat_lon_key(lat, lon)] = (lat, lon)
        except (TypeError, ValueError): pass
    if vouchers_csv is not None and Path(vouchers_csv).is_file():
        with open(vouchers_csv, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh): _add(r.get("decimalLatitude") or r.get("latitude"), r.get("decimalLongitude") or r.get("longitude"))
    if conn is not None:
        for tbl in ("vouchers", "environmental_cache"):
            if conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (tbl,)).fetchone():
                cols = {r[1] for r in conn.execute(f"PRAGMA table_info({tbl})").fetchall()}
                lc, rc = ("decimalLatitude" if "decimalLatitude" in cols else "latitude"), ("decimalLongitude" if "decimalLongitude" in cols else "longitude")
                if lc in cols and rc in cols:
                    for row in conn.execute(f"SELECT {lc}, {rc} FROM {tbl} WHERE {lc} IS NOT NULL AND {rc} IS NOT NULL"): _add(row[0], row[1])
    return pts


def download_geochem_wetland_layers(
    env_dir: Path, vouchers_csv: Optional[Path] = None, db_path: Optional[Path] = None,
    genus: Optional[str] = None, skip_nwi_network: bool = False, force: bool = False,
    http_get_fn: Optional[Callable[..., Union[bytes, str]]] = None,
) -> Dict[str, Any]:
    """Ensure DS-801 geochemistry, NWI wetland/seep cache, NWPL scores, and SQLite upsert."""
    env_dir = Path(env_dir)
    env_dir.mkdir(parents=True, exist_ok=True)
    target_genus = str(genus or _resolve_target_genus()).strip() or "Packera"
    conn = sqlite3.connect(str(db_path)) if (db_path is not None and Path(db_path).exists()) else None
    try:
        ds801_sites = ensure_usgs_ds801_soils(env_dir, force=force, http_get_fn=http_get_fn)
        nwpl = fetch_nwpl_wetland_ratings(env_dir, target_genus, conn=conn, force=force, http_get_fn=http_get_fn)
        cache_path, cache_map = env_dir / "nwi_point_cache.jsonl", {}
        if cache_path.is_file() and not force:
            for ln in cache_path.read_text(encoding="utf-8", errors="replace").splitlines():
                if ln.strip().startswith("{"):
                    try:
                        obj = json.loads(ln)
                        if obj.get("lat_lon_key"): cache_map[str(obj["lat_lon_key"])] = obj
                    except Exception: pass
        pts, nwi_results = _collect_points(vouchers_csv, conn), {}
        if pts:
            lock = threading.Lock()
            with ThreadPoolExecutor(max_workers=20) as pool:
                nwi_results = {k: fut.result() for k, fut in {k: pool.submit(query_nwi_point_or_cache, lat, lon, cache_map, cache_path, lock, skip_nwi_network, http_get_fn) for k, (lat, lon) in pts.items()}.items()}
        n_upserted = 0
        if conn is not None and pts and conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='environmental_cache'").fetchone():
            cols = {r[1] for r in conn.execute("PRAGMA table_info(environmental_cache)").fetchall()}
            for col, typ in _ENV_GEOCHEM_DDL:
                if col not in cols:
                    conn.execute(f"ALTER TABLE environmental_cache ADD COLUMN {col} {typ}")
            rows_sql = [
                (k, lat, lon, (g := compute_idw_soil_geochem(lat, lon, ds801_sites))["soil_ca_mg_molar_ratio"], g["soil_ultramafic_ni_cr_index"], g["soil_carbonate_c_pct"],
                 (n := nwi_results.get(k, {})).get("is_nwi_wetland"), n.get("nwi_wetland_system"), n.get("nwi_wetland_dist_m"))
                for k, (lat, lon) in pts.items()
            ]
            set_sql = ", ".join(f"{c}=COALESCE(excluded.{c}, environmental_cache.{c})" for c, _ in _ENV_GEOCHEM_DDL)
            conn.executemany(
                f"INSERT INTO environmental_cache (lat_lon_key, latitude, longitude, {', '.join(c for c, _ in _ENV_GEOCHEM_DDL)}) "
                f"VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT(lat_lon_key) DO UPDATE SET {set_sql}",
                rows_sql,
            )
            conn.commit()
            n_upserted = len(rows_sql)
        return {"ds801_sites": len(ds801_sites), "nwpl_species": len(nwpl), "points_processed": len(pts), "sqlite_upserted": n_upserted}
    finally:
        if conn is not None:
            conn.close()


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Provision USGS DS-801 geochemistry, USFWS NWI, and USACE NWPL layers.")
    for flg, kw in (("--env-dir", {"type": Path, "default": PROJECT_ROOT / "data" / "environmental"}), ("--vouchers", {"type": Path, "default": None}), ("--db", {"type": Path, "default": _default_db_path()}), ("--genus", {"default": _resolve_target_genus()}), ("--skip-nwi-network", {"action": "store_true", "default": False}), ("--force", {"action": "store_true", "default": False})): p.add_argument(flg, **kw)
    a = p.parse_args(list(argv) if argv is not None else None)
    print(f"[INFO] Geochem/wetland summary: {json.dumps(download_geochem_wetland_layers(a.env_dir, a.vouchers, a.db, a.genus, bool(a.skip_nwi_network), bool(a.force)))}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

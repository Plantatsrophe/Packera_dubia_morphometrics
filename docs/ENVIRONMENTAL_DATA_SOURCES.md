# Environmental, Pedological, and Geospatial Data Sources (`UNC-BOT-SOP-2026-04-REV5`)

**Provisioning Script:** `scripts/data_prep/download_environmental_layers.py`
**Extraction & Harmonization Modules:** `scripts/analysis/edaphic_fine_grain_utils.R`, `scripts/core/epa_ecoregions.py`, and `scripts/core/usda_plants_ranges.py`

## 1. Why Rhizosphere Scale and Uncertainty Gates Matter in *Packera*

Ecological divergence across North American *Packera* frequently operates at the scale of meters rather than kilometers. In the southeastern Piedmont and Ridge and Valley, narrow endemics and edaphic ecotypes colonize Piedmont granitic flatrocks (*Packera tomentosa* and *Packera anonyma* contact margins), Iredell and Mecklenburg diabase glades (*Packera schweinitziana* and *Packera paupercula* s.l. segregates), serpentine barrens (*Packera serpenticola*), and shallow Ordovician limestone glades, whereas adjacent forest matrices less than $100\text{ m}$ away support widespread mesic congeners such as *Packera aurea* or *Packera obovata*. Coarse $1\text{ km}$ global soil grids average those rock outcrops into the surrounding Ultisol or Alfisol forest matrix, erasing the exact rhizosphere chemistry and shallow lithic contact ($< 25\text{ cm}$) that maintain reproductive and habitat isolation in the field.

At first glance, pairing $10\text{ to }100\text{ m}$ soil and bedrock layers with historical herbarium records seems risky because older specimen labels often carry county-centroid or town-distance georeferences (`"5 mi NW of Chapel Hill"`). That tension resolves once georeferencing precision is treated as an explicit physical filter rather than an ignored metadata field. I gate all fine-grain rhizosphere and $30\text{ m}$ topographic extractions strictly at `coordinateUncertaintyInMeters <= 2000` m, whereas $1\text{ km}$ macro-climatic rasters permit `<= 10000` m because temperature and precipitation turn over gradually across flat regional gradients. Specimens with coarser or missing uncertainty still enter morphological and phenological evaluation, and they project passively or fall back cleanly without polluting the active Spatial Random Forest training folds.

## 2. Coordinate-Uncertainty Gates, `50 m` Outcrop Buffers, and *Locus Classicus* Controls

Step 01 (`scripts/core/harvester.py`) and Step 05 (`scripts/analysis/edaphic_fine_grain_utils.R` and `05_multimodal_spatial_rf.R`) enforce four spatial provenance rules before any environmental value enters `environmental_cache`:

1. **Dual Eligibility Gates (`is_edaphic_eligible_coord` vs. `is_clim_eligible_coord`):** Point and buffer extractions from SSURGO, POLARIS, SOLUS100, USGS SGMC, and $30\text{ m}$ 3DEP topography require finite coordinates with `coordinateUncertainty <= 2000` m (`is_edaphic_eligible_coord = 1`), and continuous weights (`edaphic_coord_weight`) down-weight observations approaching the $2000\text{ m}$ ceiling. WorldClim v2.1 bioclimatic rasters allow `coordinateUncertainty <= 10000` m (`is_clim_eligible_coord = 1`).
2. **Sub-Grid `50 m` Buffer Quantiles (`buffer_radius_m = 50`):** Even GPS-determined herbarium vouchers collected on a narrow flatrock margin or roadcut ledge can sit one pixel off a mapped lithic contact. Around each eligible coordinate, I extract both the point value and `50 m` circular buffer extremes (`ph_q10_buffer`, `ph_q90_buffer`, `min_bedrock_depth_buffer_cm`, `max_cfvo_buffer_pct`, and `azonal_lithology_presence_buffer`) so narrow outcrop exposures register cleanly.
3. **Shallow Fibrous Rhizosphere Integration (`0 to 15 cm`):** *Packera* individuals root via shallow fibrous stolons and horizontal rhizomes concentrated in the upper $15\text{ cm}$ of mineral soil rather than deep taproots. All multi-depth soil rasters integrate strictly across the $0\text{ to }15\text{ cm}$ root zone using trapezoidal horizon weights:
   $$X_{0\text{ to }15\text{ cm}} = \frac{5 \cdot X_{0\text{ to }5\text{ cm}} + 10 \cdot X_{5\text{ to }15\text{ cm}}}{15}$$
   where discrete SOLUS100 point-depth predictions ($0\text{ cm}$, $5\text{ cm}$, and $15\text{ cm}$) first average adjacent boundaries ($X_{0\text{ to }5} = 0.5(X_0 + X_5)$ and $X_{5\text{ to }15} = 0.5(X_5 + X_{15})$).
4. **Historical *Locus Classicus* Fallback:** Eighteenth- and nineteenth-century nomenclatural types rarely carry GPS precision, yet their protologue localities anchor species names under the Shenzhen Code (Turland et al. 2018). When a verified type specimen lacks `<= 2000 m` coordinates, `compute_type_locality_ecological_control()` estimates its *locus classicus* baseline from `Tier_1_Gold` conspecific anchors located within $150\text{ km}$ or inside the identical US EPA Level III Ecoregion (`regional_group`).

## 3. Automated Layer Provisioning (`download_environmental_layers.py`)

Before Step 05 begins, `run_pipeline.sh` and `main.py synthesis` automatically invoke `scripts/data_prep/download_environmental_layers.py` (unless `--skip-env-download` or `--demo` is set) to verify and cache all open-access continental soil, bedrock, and cohort-targeted elevation layers on disk:

```bash
# Provision all default layers (SOLUS100 + USGS SGMC + cohort-targeted 3DEP 30 m DEM tiles)
.venv/bin/python scripts/data_prep/download_environmental_layers.py \
  --vouchers runs/latest/data/tables/curated_vouchers.csv \
  --env-dir data/environmental \
  --workers 4
```

Every download streams into a temporary `.part` file, validates GDAL/GeoPackage readability, and renames atomically into `data/environmental/`. If network access is unavailable or a tile directory is absent when `--demo` is off, Step 05 logs a clear notice and records `NA` for the missing layer rather than fabricating synthetic raster values.

## 4. Layer-by-Layer Provenance and Harmonization Specifications

### 4.1 USDA-NRCS SSURGO via Soil Data Access (`10 to 30 m` Vector Pedology)

* **Provider & License:** USDA Natural Resources Conservation Service (NRCS) Soil Survey Geographic Database (`SSURGO`) queried via the REST Soil Data Access (`SDA`) PostGIS endpoint (`https://SDMDataAccess.sc.egov.usda.gov/Tabular/post.rest`); U.S. federal public domain (`CC0 1.0`).
* **Local Cache:** `data/environmental/ssurgo_cache/` and SQLite `environmental_cache`.
* **Extracted Variables & Outcrop Classification:** Component-weighted map unit symbol (`musym`), map unit name (`muname`), Soil Taxonomy suborder and great group (`taxsuborder`, `taxgreatgroup`), depth to lithic or paralithic bedrock (`ssurgo_bedrock_depth_cm` from `corestrictions`), rock outcrop component percentage (`ssurgo_rock_outcrop_pct`), available water storage (`ssurgo_aws_0_50cm`), and $0\text{ to }15\text{ cm}$ horizon-weighted sand, clay, $\text{pH}_{\text{H}_2\text{O}}$, organic matter, and $\text{CEC}_7$. Map unit names and component parent materials classify into six discrete `outcrop_class` states: `Granitic flatrock`, `Limestone glade`, `Sandstone/Siliceous`, `Serpentine`, `Diabase/Mafic glade`, and `Non-outcrop matrix`.

### 4.2 USDA-NRCS SOLUS100 (`100 m` CONUS Soil & Lithic Contact Grids)

* **Provider & Citation:** Nauman, T. W., S. Kienast-Brown, S. M. Roecker, J. A. Brungard, D. White, J. Philippe, and J. A. Thompson. 2024. Soil Landscapes of the United States 100-meter (`SOLUS100`) soil property maps. Ag Data Commons, USDA-NRCS (`https://storage.googleapis.com/solus100pub`), `doi: 10.15482/USDA.ADC/25033856` (`CC0 1.0`).
* **Local Cache:** `data/environmental/solus_100m/` (20 Cloud-Optimized GeoTIFFs, ~4 to 8 GB total).
* **Cached Raster Inventory & Integer Storage Divisors:** SOLUS100 stores predicted rasters as scaled integers to reduce file size, so `edaphic_fine_grain_utils.R` divides raw pixel values by each layer's exact scalar divisor before trapezoidal $0\text{ to }15\text{ cm}$ integration:

| Property Group | Cached SOLUS100 GeoTIFF Filenames | Raw Divisor | Physical Units | Harmonized Target Column |
| :--- | :--- | :---: | :---: | :--- |
| **Lithic & Restriction Depth** | `anylithicdpt_cm_p.tif`, `resdept_all_cm_p.tif` | `1` | $\text{cm}$ | `solus_resdept_cm`, `fine_bedrock_depth_cm` |
| **Coarse Rock Fragments** | `fragvol_0_cm_p.tif`, `fragvol_5_cm_p.tif`, `fragvol_15_cm_p.tif` | `1` | $\% \text{ vol}$ | `solus_fragvol_pct`, `fine_cfvo_pct` |
| **Soil Reaction ($\text{pH}$)** | `ph1to1h2o_0_cm_p.tif`, `ph1to1h2o_5_cm_p.tif`, `ph1to1h2o_15_cm_p.tif` | `100` | $1:1\ \text{H}_2\text{O}\ \text{pH}$ | `fine_soil_ph` |
| **Sand Fraction** | `sandtotal_0_cm_p.tif`, `sandtotal_5_cm_p.tif`, `sandtotal_15_cm_p.tif` | `1` | $\% \text{ wt}$ | `solus_sand_pct`, `fine_soil_sand` |
| **Clay Fraction** | `claytotal_0_cm_p.tif`, `claytotal_5_cm_p.tif`, `claytotal_15_cm_p.tif` | `1` | $\% \text{ wt}$ | `solus_clay_pct`, `fine_soil_clay` |
| **Silt Fraction** | `silttotal_0_cm_p.tif`, `silttotal_5_cm_p.tif`, `silttotal_15_cm_p.tif` | `1` | $\% \text{ wt}$ | `solus_silt_pct`, `soil_silt_pct` |
| **Cation Exchange Capacity** | `cec7_0_cm_p.tif`, `cec7_5_cm_p.tif`, `cec7_15_cm_p.tif` | `10` | $\text{cmol}_+\cdot\text{kg}^{-1}$ | `solus_cec`, `fine_soil_cec` |

* **Isometric Log-Ratio (`ILR`) Compositional Soil Texture Transform:** Because `% sand + % silt + % clay = 100%` lies on a closed $D = 3$ simplex, raw percentages induce spurious negative collinearity in linear and tree-based models. After multiplicative zero replacement ($\delta = 0.05\%$) and closure, Step 05 maps the three particle fractions onto two unconstrained Euclidean coordinates (Egozcue et al. 2003, `doi: 10.1023/A:1023818214614`):
  $$\text{ilr\_texture\_1} = \sqrt{\frac{2}{3}} \ln\left(\frac{\text{sand}}{\sqrt{\text{silt} \cdot \text{clay}}}\right), \qquad \text{ilr\_texture\_2} = \sqrt{\frac{1}{2}} \ln\left(\frac{\text{silt}}{\text{clay}}\right)$$

### 4.3 POLARIS Probabilistic Soil Series (`30 m` / `1 arcsec` Optional Local Tiles)

* **Provider & Citation:** Chaney, N. W., B. Minasny, J. D. Herman, T. W. Nauman, C. W. Brungard, C. L. S. Morgan, A. B. McBratney, E. F. Wood, and Y. Yimam. 2019. POLARIS soil properties: 30-m probabilistic maps of soil properties over the contiguous United States. *Water Resources Research* 55: 2916 to 2938, `doi: 10.1029/2018WR022797`.
* **Local Cache & Provisioning Policy:** `data/environmental/polaris_30m/`. Because full CONUS POLARIS mirrors span ~1,800 tiles per property-depth combination and SOLUS100 + SSURGO SDA + SGMC + $30\text{ m}$ 3DEP already cover continental chemistry and outcrop lithology, `download_environmental_layers.py` intentionally skips bulk POLARIS downloads by default. Whenever `1 x 1 deg` POLARIS GeoTIFFs are staged locally in `data/environmental/polaris_30m/`, `edaphic_fine_grain_utils.R` detects them automatically and gives them top priority in the harmonized cascade (`POLARIS 30 m -> SOLUS100 100 m -> SSURGO SDA -> SoilGrids`).

### 4.4 USGS State Geologic Map Compilation (`SGMC`, Mapped Bedrock Lithology)

* **Provider & Citation:** Horton, J. D., C. A. San Juan, and D. B. Stoeser. 2017. The State Geologic Map Compilation (`SGMC`) geodatabase of the conterminous United States (ver. 1.1). U.S. Geological Survey Data Series 1052, ScienceBase item `5888bf4fe4b05ccb964bab9d`, `doi: 10.5066/F7WH2N65` (`CC0 1.0`).
* **Local Cache:** `data/environmental/sgmc/sgmc_conus_lithology.gpkg` (pre-joined GeoPackage created by joining `SGMC_Geology` polygons with `SGMC_Units` on `UNIT_LINK` in `WGS84`).
* **Six-Class Azonal & Zonal Bedrock Taxonomy (`sgmc_lithology_class`):**
  1. `Ultramafic_Serpentine` (`is_azonal_lithology = 1`): serpentinite, dunite, peridotite, soapstone, and talc-tremolite bodies (elevated $\text{Mg}:\text{Ca}$ and $\text{Ni}/\text{Cr}$ stress).
  2. `Mafic_Igneous_Metamorphic` (`is_azonal_lithology = 1`): diabase dikes/sills, gabbro, basalt, greenstone, and amphibolite (high $\text{Ca}/\text{Mg}$ shrink-swell Iredell clays).
  3. `Carbonate_Calcareous` (`is_azonal_lithology = 1`): limestone, dolostone, marble, travertine, marl, and chalk (circumneutral to alkaline shallow glades).
  4. `Felsic_Granitic_Gneiss` (`is_azonal_lithology = 0`): granite, granodiorite, monzonite, rhyolite, and felsic gneiss (acidic Piedmont flatrocks and domes).
  5. `Siliciclastic_Consolidated` (`is_azonal_lithology = 0`): sandstone, quartzite, conglomerate, arkose, shale, siltstone, slate, phyllite, and mica schist.
  6. `Unconsolidated_Sedimentary` (`is_azonal_lithology = 0`): Coastal Plain sands, alluvial terraces, loess, glacial till, and outwash.

### 4.5 USGS 3DEP / AWS Terrain Skadi `1-arcsec` (`~30 m`) Elevation Tiles

* **Provider & Access:** U.S. Geological Survey 3D Elevation Program (`3DEP`) $1\text{ arcsec}$ ($\sim 30\text{ m}$) seamless DEM served via AWS Open Data Terrain Tiles (`https://s3.amazonaws.com/elevation-tiles-prod/skadi`, SRTM HGT format converted to DEFLATE-compressed GeoTIFFs `N35W079.tif`) with direct fallback to the USGS National Map `3DEPElevation/ImageServer/exportImage` endpoint (`CC0 1.0`).
* **Local Cache:** `data/environmental/dem_3dep_30m/` (cohort-targeted: downloads strictly the $1^\circ \times 1^\circ$ tiles intersected by `curated_vouchers.csv`, typically 40 to 180 MB per cohort).
* **Derived Micro-Topographic Indices:**
  * `elevation_m` and `topo_slope_deg`: Point elevation ($\text{m}$) and Horn $3 \times 3$ terrain slope ($\text{deg}$).
  * `topo_tpi_30m`: Topographic Position Index ($z_0 - \bar{z}_{3 \times 3}$), separating exposed convex rock domes and bluff rims ($\text{TPI} > 0$) from sheltered seepage coves ($\text{TPI} < 0$).
  * `topo_hli_30m`: Potential direct incident solar radiation and southwest afternoon thermal stress (`[0, 1]`) folded around $225^\circ$ aspect (McCune and Keon 2002, `doi: 10.1111/j.1654-1103.2002.tb02087.x`).
  * `topo_twi_30m`: D8 flow-accumulation Topographic Wetness Index $\ln(a / \tan\beta)$ separating xeric outcrops from saturated fen and bottomland habitats (Beven and Kirkby 1979, `doi: 10.1080/02626667909491834`).

### 4.6 WorldClim v2.1 Bioclimatic Grids (`30 arcsec` / `~1 km`)

* **Provider & Citation:** Fick, S. E., and R. J. Hijmans. 2017. WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology* 37: 4302 to 4315, `doi: 10.1002/joc.5086`.
* **Local Cache:** `data/environmental/` (`wc2.1_30s_bio_1.tif` through `wc2.1_30s_bio_19.tif`), automatically cropped on first download to the North American bounding box `[-130, -50, 14, 66]` with `DEFLATE` compression to minimize disk and memory overhead during MaxEnt (`step05_sdm_maxent.R`) raster projection.

### 4.7 US EPA Level III Ecoregions and USDA PLANTS County Distributions

* **US EPA Level III Ecoregions (`scripts/core/epa_ecoregions.py`):** Omernik, J. M., and G. E. Griffith. 2014. Ecoregions of the conterminous United States: evolution of a hierarchical spatial framework. *Environmental Management* 54: 1249 to 1266, `doi: 10.1007/s00267-014-0364-1`. Cached in `data/environmental/epa_ecoregions_l3/`, assigning `regional_group` via exact point-in-polygon intersection for georeferenced sheets or $\ge 95\%$ single-ecoregion county-area dominance for un-georeferenced historical vouchers.
* **USDA PLANTS County Ranges (`scripts/core/usda_plants_ranges.py`):** USDA-NRCS PLANTS Database (`https://plants.usda.gov`, `CC0 1.0`). Resolves accepted symbols and evaluates each georeferenced collection against county distributions with a `50 km` county-border buffer (`usda_county_status` in `{'in_range', 'buffered_50km', 'out_of_range', 'state_only_match', 'unverified'}`). By default (`strict_usda_range_filter: false`), out-of-range vouchers remain in the cohort as high-value curatorial audit candidates rather than being silently dropped, as novel range extensions and misfiled folder sheets both appear outside historical county maps.

## 5. Literature Cited

1. Beven, K. J., and M. J. Kirkby. 1979. A physically based, variable contributing area model of basin hydrology. *Hydrological Sciences Bulletin* 24: 43 to 69, `doi: 10.1080/02626667909491834`.
2. Chaney, N. W., et al. 2019. POLARIS soil properties: 30-m probabilistic maps of soil properties over the contiguous United States. *Water Resources Research* 55: 2916 to 2938, `doi: 10.1029/2018WR022797`.
3. Egozcue, J. J., V. Pawlowsky-Glahn, G. Mateu-Figueras, and C. Barceló-Vidal. 2003. Isometric logratio transformations for compositional data analysis. *Mathematical Geology* 35: 279 to 300, `doi: 10.1023/A:1023818214614`.
4. Fick, S. E., and R. J. Hijmans. 2017. WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology* 37: 4302 to 4315, `doi: 10.1002/joc.5086`.
5. Horton, J. D., C. A. San Juan, and D. B. Stoeser. 2017. The State Geologic Map Compilation (`SGMC`) geodatabase of the conterminous United States (ver. 1.1). U.S. Geological Survey Data Series 1052, `doi: 10.5066/F7WH2N65`.
6. McCune, B., and D. Keon. 2002. Equations for potential annual direct incident radiation and heat load. *Journal of Vegetation Science* 13: 603 to 606, `doi: 10.1111/j.1654-1103.2002.tb02087.x`.
7. Nauman, T. W., et al. 2024. Soil Landscapes of the United States 100-meter (`SOLUS100`) soil property maps. Ag Data Commons, USDA-NRCS, `doi: 10.15482/USDA.ADC/25033856`.
8. Omernik, J. M., and G. E. Griffith. 2014. Ecoregions of the conterminous United States: evolution of a hierarchical spatial framework. *Environmental Management* 54: 1249 to 1266, `doi: 10.1007/s00267-014-0364-1`.
9. Turland, N. J., et al. 2018. *International Code of Nomenclature for algae, fungi, and plants (Shenzhen Code)*. Regnum Vegetabile 159. Koeltz Botanical Books, Glashütten, `doi: 10.12705/Code.2018`.

# PEF pipeline (Brightway 2.5 + wurst)

Applies rules of the PEF method (EF Recommendation (EU) 2021/2279, Annex I and IV, PEFCR dairy) to
Brightway databases and compares EF 3.1 results with and without them.

| Script | What it does | Data |
|---|---|---|
| `pef_allocation.py` | PEF allocation rules, CFF, single score (library) | wurst-format dicts |
| `run_demo.py` | dairy chain on a toy background | own project `pef_pipeline_demo` |
| `run_bafu.py` | dairy foreground on BAFU | `bafu test`: writes `pef-dairy-foreground` |
| `bafu_pef_compare.py` | BAFU vs. PEF version for 70 food and feed products | `bafu test` (read only) |
| `ecoinvent_pef_realloc.py` | PEF re-allocation of dairy co-products, 877 food markets | `salton sea pef` (read only) |
| `bafu_workbook_check.py` | our BAFU baseline vs. BAFU's published EF 3.1 results | `reports/bafu_pef_compare.csv` + release workbook |
| `hestia_db.py` | standalone HESTIA crop database on BAFU, PEF allocation | `bafu test`: writes `hestia-2026-04-30` |
| `export_hestia_db.py` | own project `hestia-bafu-pef` + portable files (Brightway write format, Excel) | writes `reports/hestia_bw/` |
| `import_hestia_db.py` | writes the portable files into any project, Brightway 2.5 or legacy 2 | any project |
| `ef31_nw.csv` | EF 3.1 normalisation and weighting factors | normalisation factors not yet verified |

Results go to `reports/`, which stays local because it holds numbers derived from licensed data
(BAFU, ecoinvent).

## PEF rules: implemented and not implemented

### Implemented

| PEF rule | Code | Applied to | Notes |
|---|---|---|---|
| Pesticides to 90 % soil / 9 % air / 1 % fresh water | `bafu_pef_compare.py` | bafu-2026: 232 active ingredients, 189 datasets | Pyraclostrobin stays in soil: no air or water flow in the EF biosphere. |
| No negative heavy-metal balances to soil | `bafu_pef_compare.py` | bafu-2026: 107 datasets | Judgement call; PEF does not state it explicitly. |
| Water use with country-level factors | `bafu_pef_compare.py` (layer `water`) | bafu-2026: 650 datasets | Country factors from the JRC EF 3.1 CF file; regions without one (RER, GLO …) keep the global factor. |
| Dairy farm, IDF biophysical allocation: AF_milk = 1 − 6.04 · M_meat / M_milk | `rule_dairy_farm`, `ecoinvent_pef_realloc.py` | ecoinvent milk production (RoW, CA-QC, ZA); foreground demos | ecoinvent milk is raw milk, not FPCM. |
| Dairy plant, allocation by dry matter (PEFCR dairy) | `rule_physical`, `ecoinvent_pef_realloc.py` | ecoinvent cheese, yogurt, butter | — |
| EF 3.1 characterisation, 16 categories | all scripts | sentier EF v3.1 (BAFU), ecoinvent EF v3.1 | Sub-indicators are excluded so nothing is counted twice. |
| Normalisation, weighting, single score | `single_score`, `ef31_nw.csv` | BAFU comparison | Weights match published values; normalisation factors still to be checked against the JRC EF 3.1 package. |
| Slaughterhouse, economic allocation | `rule_economic` | foreground demo only | PEF default prices and mass fractions not entered; no multi-output slaughterhouse in BAFU or ecoinvent. |
| Crop / residue, economic allocation (LEAP) | `rule_economic` | — | Not needed: BAFU and ecoinvent are already economic. |
| Circular Footprint Formula | `cff_material`, `cff_energy_disposal`, `add_cff_exchanges` | PET tray in the dairy foreground | Only A = 0.5 is a PEF default; R1, R2, R3 and Qs/Qp are placeholders. |

### Not implemented

| PEF rule | Reason |
|---|---|
| N field emissions, PEF default (IPCC Tier 1: N₂O 1 % of N, NH₃ 10 % / 20 %, NO₃⁻ 30 %) | Feasible. On BAFU it needs assumed N contents for slurry and manure. |
| Small-ruminant default allocation | ecoinvent "sheep production, for meat" (meat vs. wool) not covered yet. |
| Slaughterhouse defaults on real data | No multi-output slaughterhouse in BAFU or ecoinvent. |
| Re-allocation of BAFU | bafu-2026 has no prices, production volumes or properties; its allocation rules are in the LCI reports of the documentation component. |
| Country-specific factors for acidification and eutrophication | The sentier EF 3.1 method ships only global factors; water use is regionalised by the `water` layer. |
| Long-term emissions | Whether PEF excludes them still needs checking. |
| Rice CH₄, drained peat, LUC within 20 years, biodiversity reporting | Need new activity data, not a rule. |
| Climate change split into fossil, biogenic and LUC | Reporting only; the sub-indicators exist but are not reported. |
| Data-quality rating, compliance class per dataset (transitional guidance, 14 July 2026) | Neither BAFU nor ecoinvent is EF-compliant, so a study on them needs the disclaimer. |
| Electricity country mix, cut-off ≤ 3 %, one reference flow, data needs matrix | Rules for building a product study, not for re-modelling a database. |

## How the pipeline uses wurst

[wurst](https://github.com/polca/wurst), the package premise builds on, holds Brightway databases as
plain lists of dicts. Here it is used for the **foreground only**:

- `split_multioutput` and the CFF helpers read and write wurst-format dicts.
- `run_demo.py` and `run_bafu.py` write the allocated foreground with `wurst.write_brightway2_database`,
  which also relinks internal exchanges and checks that every link resolves.
- `wurst.extract_brightway2_databases` (also copied into bw2data as `bw2data.backends.wurst_extraction`)
  assumes single-output datasets and drops custom exchange fields. The multi-output foreground is
  therefore built as dicts, not extracted from Brightway.

The BAFU and ecoinvent comparisons do **not** use wurst. Both changes are linear, so they are applied
to the matrices of one LCA:

- BAFU (pesticides, heavy metals): change in score = Σ_j s_j · Σ_f ΔB_fj · CF_f, on the unchanged
  supply vector s.
- ecoinvent (re-allocation): each co-product column of the technosphere and biosphere matrices is
  multiplied by AF_PEF / AF_economic, and the production entry is kept.

This is exact, fast, and writes no copy of a licensed database into the project. A persistent
"PEF version" database (e.g. for Activity Browser) would use wurst: extract, apply the same rules as
dict transforms, write.

`bw-functional` (`multifunctional`) is an alternative for the foreground: it stores co-products and
does property-based allocation inside Brightway, but cannot do the IDF formula or the CFF, and it is
still beta.

## Running it

```bash
python run_demo.py
```

```bash
python run_bafu.py
```

```bash
python bafu_pef_compare.py
```

```bash
python ecoinvent_pef_realloc.py
```

```bash
python bafu_workbook_check.py "<path>/BAFU-2026 v1 LCIA Results_corrected.xlsx"
```

- `run_demo.py` builds its own project with a GWP-only method; all numbers are illustrative.
- `run_bafu.py` (re)writes `pef-dairy-foreground` in `bafu test` with illustrative farm data.
- `ecoinvent_pef_realloc.py` needs a Brightway 2.5 project. `salton sea pef` is a migrated copy of
  `salton sea`, so the original project stays unchanged.
- The BAFU workbook comes from <https://nexus.openlca.org/database/BAFU> after accepting the terms of
  use; it cannot be downloaded by script.

## Main results

**BAFU vs. PEF version (70 food and feed products, against the corrected baseline).** Only toxicity
and water use change. Freshwater ecotoxicity rises for 63 products (median ×1.32, 20 products more
than doubled, asparagus ×40–119 from deltamethrin in the seedlings). Human toxicity changes for 29–35
products, mainly through the heavy-metal layer, and all 23 negative toxicity results disappear. The
`water` layer changes water use for 65 products (median ×1.04; green beans CH ×0.035, palm MY ×0.08,
sorghum CN ×0.09). Single score: median +0.9 %, 12 products up more than 10 %, 4 more than doubled,
and down to ×0.56 where water use falls.

**ecoinvent re-allocation.** About 30 of 893 products change, in all 16 categories together.
Climate change per kg: whey ×10.0, fresh cheese ×0.30, buttermilk ×0.40, butter ×1.20, cow milk ×0.97,
cattle for slaughter ×1.02. Knock-on: red meat ×1.01–1.10, feed and chicken ×0.97.

## HESTIA database on BAFU (`hestia_db.py`)

`hestia-2026-04-30` holds the HESTIA aggregated release 2026-04-30 (CC-BY 4.0) as Brightway
datasets: 367 crop and processed-crop cycles (latest period per product and country, 48 countries),
683 datasets after allocation. It links only to `bafu-2026` and `ef-3.1-biosphere`, so the system
contains no ecoinvent.

| Part | Modelling |
|---|---|
| Products | All products with an economic value share; 190 multi-output cycles (grain + straw, oil + meal, rice milling) split by `split_multioutput`, rule `crop` (economic, LEAP) or `economic` for processing. Co-products under 0.5 % of the main product's mass (intercrop artefacts of the aggregation) are dropped. |
| Inputs | Fertilisers per nutrient, diesel/petrol/LPG/fuel oil/natural gas as supply only (HESTIA carries the combustion emissions), grid electricity by country ("Electricity mix", GLO fallback), machinery, lime, tap water, processing chemicals → bafu-2026. Seed and crop inputs of processing cycles → HESTIA datasets. Manure and compost: no burden (residues). |
| Emissions | HESTIA on-farm emissions → EF 3.1 (non-urban air, fresh water). Fossil vs. biogenic vs. land-use-change CO₂/CH₄/CO split by emission source; dLUC, drained organic soils and vegetation burning count as land use change. HESTIA background emissions are dropped. |
| Pesticides | PEF default: 90 % of the applied active ingredient to agricultural soil, 9 % to air, 1 % to fresh water, matched by CAS; plus "Pesticide unspecified" production per kg AI. |
| Land, water | Occupation (site + unused duration) as Arable, Arable irrigated or Permanent crops; irrigation withdrawals as ground, river or lake water, taken as consumed. |

**Gaps** (`reports/hestia_unmapped.csv`): pesticides reported only as unspecified, herbicide,
insecticide or fungicide (about 1 t AI over all cycles), 75 named active ingredients without an EF
flow, PM10/TSP/black carbon (no EF factor), coal, wood fuel, sulphur, gypsum, branded NPK products,
saplings and seedlings (counted in pieces).

**Data issue:** HESTIA's almond cycle for Spain reports 4 kg/ha deltamethrin, roughly 200 times a
normal dose, which dominates its freshwater ecotoxicity.

Scores per product (`--score`) go to `reports/hestia_ef31_scores.csv`, allocation factors to
`reports/hestia_allocation.csv`.

**Using it in Activity Browser.** `export_hestia_db.py` creates the project `hestia-bafu-pef`
(Brightway 2.5, Activity Browser 3) with only the HESTIA database, BAFU, the EF 3.1 biosphere and the 25
EF 3.1 methods, and writes portable files to `reports/hestia_bw/`. `import_hestia_db.py <project>` loads
them into any project; run it with the legacy environment for Activity Browser 2 (`hestia-bafu-pef-AB2`
was made that way). BW2Packages are not used: Brightway 2.5 packages store characterisation factors by
internal id and do not open in legacy Brightway. A round trip into both versions gives identical scores.

## Validation

### Characterisation factors vs. JRC

The sentier EF v3.1 method in `bafu test` was compared with the JRC file
`EF-LCIAMethod_CF(EF-v3.1).xlsx`: all 25 indicators match every global factor (about 280,000) with no
deviation above 1 %. JRC also publishes about 300 country-specific factors per country for
acidification, eutrophication and water use, which the sentier method does not include.

### Our BAFU baseline vs. BAFU's published EF 3.1 results

`bafu_workbook_check.py` compares the baseline of `bafu_pef_compare.py` with
`BAFU-2026 v1 LCIA Results_corrected.xlsx` for the 70 products (all present in the workbook). The
baseline includes one correction of the Brightway import (`fix`): unspecified chromium emitted to soil
and water is mapped to Cr(III), as BAFU's own EF mapping does (731 datasets).

| Category | Within 1 % / 10 %, import as is | With the chromium fix | Median ours / published | Remaining cause |
|---|---|---|---|---|
| Climate change, acidification, eutrophication, PM, ozone, radiation, photochemical ozone, resources | 67–69 / 68–69 of 69 | unchanged | 1.000 | — |
| Human toxicity, non-cancer | 64 / 69 | unchanged | 1.000 | small |
| Human toxicity, cancer | 1 / 7 | **4 / 31** | 1.70 → **1.08** | Other flows (sunflower still ×8.6; sorghum CN sign differs). |
| Freshwater ecotoxicity | 31 / 47 | **40 / 53** | 1.00 → 0.99 | Flows characterised in BAFU's mapping but not in the import (rape oil ×0.10). |
| Land use | 43 / 51 | unchanged | 1.00 | Land transformation "To / From Arable" nearly cancel; ×5–6 for organic and extensive cereals. |
| Water use | 6 / 63 | unchanged | 1.05 | Turbine water and released water both characterised (±43 m³/m³) and nearly cancel; green beans CH ×29. |

Whey is left out: BAFU gives it no burden, so both values are about zero.

**Water use is not part of the fix.** Country factors for every dataset (the `water` layer) make
asparagus PE and MX and green beans KE match BAFU exactly, but lower agreement overall from 63 to 27
products within 10 % (palm MY, sorghum CN drop to ×0.08). BAFU regionalises water only in some
datasets, which the Brightway import does not show. Country-level water use is a PEF requirement, so
it is a PEF layer instead.

**What it means.** The baselines behind the main PEF findings match: asparagus (CH, MX, PE) and
soybeans BR within 0.2 % for freshwater ecotoxicity. The BAFU and PEF versions use the same import,
so the PEF ratios hold. Land use for organic and extensive cereals and water use for some vegetables
still differ from BAFU's published numbers.

**Still open in the import** (would also apply to the sentier loader): turbine water, the remaining
ecotoxicity flows (rape oil), and land transformation.

## Open points

- Verify the EF 3.1 normalisation factors in `ef31_nw.csv` against the JRC EF 3.1 reference package
  (not in the CF file, and not shipped by `sentier-methods`).
- Implement the PEF N field model.
- Enter PEF default values: slaughterhouse prices and mass fractions; CFF R1, R2, R3, Qs/Qp from
  Annex C.
- Characterise ecoinvent with the same EF 3.1 factors as BAFU (e.g. a flowmapper bridge) so both
  databases can be compared directly.

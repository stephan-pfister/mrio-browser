# PEF pipeline (Brightway 2.5 + wurst)

Applies rules of the PEF method (Recommendation (EU) 2021/2279, Annex I: 4.4.1 agricultural production,
4.4.8 circular footprint formula, 4.5 multifunctionality; default values in Part C of Annex II) and of the
dairy PEFCR (February 2025) to Brightway databases and compares EF 3.1 results with and without them.

| Script | What it does | Data |
|---|---|---|
| `pef_allocation.py` | PEF allocation rules, CFF, single score (library) | wurst-format dicts |
| `run_demo.py` | dairy chain on a toy background | own project `pef_pipeline_demo` |
| `run_bafu.py` | dairy foreground on BAFU | `bafu test`: writes `pef-dairy-foreground` |
| `bafu_pef_compare.py` | BAFU vs. PEF version for 100 food and feed products (list in `reports/bafu_food_products.csv`) | `bafu test` (read only) |
| `ecoinvent_pef_realloc.py` | PEF re-allocation of dairy co-products, 877 agricultural, food and beverage markets | `salton sea pef` (read only) |
| `bafu_workbook_check.py` | our BAFU baseline vs. BAFU's published EF 3.1 results | `reports/bafu_pef_compare.csv` + release workbook |
| `hestia_db.py` | standalone HESTIA crop database on BAFU, PEF allocation | `bafu test`: writes `hestia-2026-04-30` |
| `export_hestia_db.py` | own project `hestia-bafu-pef` + portable files (Brightway write format, Excel) | writes `reports/hestia_bw/` |
| `import_hestia_db.py` | writes the portable files into any project, Brightway 2.5 or legacy 2 | any project |
| `ef31_nw.csv` | EF 3.1 normalisation and weighting factors | JRC `Normalisation_Weighting_Factors_EF_3.1.xlsx` |

Results go to `reports/`, which stays local because it holds numbers derived from licensed data
(BAFU, ecoinvent).

## PEF rules: implemented and not implemented

### Implemented

| PEF rule | Code | Applied to | Notes |
|---|---|---|---|
| Pesticides to 90 % agricultural soil / 9 % air / 1 % fresh water | `bafu_pef_compare.py` | bafu-2026: 232 active ingredients, 169 datasets | PEF default (Annex I 4.4.1.4). A PEFCR may set otherwise: the dairy PEFCR keeps 100 % to soil, as BAFU does. Active ingredients that BAFU emits to non-agricultural soil in agriculture datasets go 90 % to agricultural soil. Sulfur and sulfuric acid count as active ingredients only in agriculture datasets (fungicide, haulm desiccant). Pyraclostrobin stays in soil (no air or water flow in the EF biosphere). EF duplicates (Mecoprop-P, Paraquat) are resolved by CAS number and UUID series. |
| No negative heavy-metal balances to soil | `bafu_pef_compare.py` | bafu-2026: 106 datasets (net per dataset and flow) | Partial version of Annex I 4.4.1.6 approach (a): for food crops, crop uptake shall not be modelled at all. Positive balances keep their uptake, and feed crops (approach b, uptake counted when the manure release is modelled) are treated the same. |
| Water use with country-level factors | `bafu_pef_compare.py` (layer `water`) | bafu-2026: 650 datasets | Country factors from the JRC EF 3.1 CF file; regions without one (RER, GLO …) keep the global factor. |
| Dairy farm, IDF biophysical allocation: AF_milk = 1 − 6.04 · M_meat / M_milk | `rule_dairy_farm`, `ecoinvent_pef_realloc.py` | ecoinvent milk production (RoW, CA-QC, ZA); foreground demos | ecoinvent milk is raw milk, not FPCM. |
| Dairy plant, allocation by dry matter (PEFCR dairy) | `rule_physical`, `ecoinvent_pef_realloc.py` | ecoinvent cheese, yogurt, butter | — |
| EF 3.1 characterisation, 16 categories | all scripts | sentier EF v3.1 (BAFU), ecoinvent EF v3.1 | Sub-indicators are excluded so nothing is counted twice. |
| Normalisation, weighting, single score | `single_score`, `ef31_nw.csv` | BAFU comparison | Factors from the JRC EF 3.1 file, identical to Tables 61 and 62 of the dairy PEFCR. Until 9 Nov 2026 the two human toxicity normalisation factors were EF 3.0 values. |
| Slaughterhouse, economic allocation | `rule_economic` | foreground demo only | PEF default prices and mass fractions not entered; no multi-output slaughterhouse in BAFU or ecoinvent. |
| Crop / residue, economic allocation (LEAP) | `rule_economic` | — | Not applied to BAFU or ecoinvent. ecoinvent's price allocation is verified only for the dairy co-products; BAFU's allocation keys are in its LCI reports, not in the data. |
| Circular Footprint Formula | `cff_material`, `cff_energy_disposal`, `add_cff_exchanges` | PET tray in the dairy foreground | A = 0.5 is the PEF value when Part C of Annex II has no application- or material-specific A (not checked for PET trays). R1 = 0 matches the dairy PEFCR default for PET cups. R2, R3 and Qs/Qp are placeholders. Substituted energy is converted from MJ to the unit of the energy dataset (kWh for electricity). On BAFU, recycling is burden-free: the only PET recycling dataset (`xx Recycling PET`) is obsolete and empty. |

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

- BAFU (chromium fix, pesticides, heavy metals, water factors): change in score =
  Σ_j s_j · Σ_f ΔB_fj · CF_f, on the unchanged supply vector s.
- ecoinvent (re-allocation): each co-product column of the technosphere and biosphere matrices is
  multiplied by AF_PEF / AF_economic, and the production entry is kept.

This is exact, fast, and writes no copy of a licensed database into the project. A persistent
"PEF version" database (e.g. for Activity Browser) would use wurst: extract, apply the same rules as
dict transforms, write.

`bw-functional` (`multifunctional`) is an alternative for the foreground: it stores co-products and
offers equal, property-based and manual allocation inside Brightway. An IDF factor could be passed as
a manual factor, but the rule logic and the CFF stay outside, and the package is still beta.

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

**BAFU vs. PEF version (100 food and feed products, against the corrected baseline).** Only toxicity
and water use change. Freshwater ecotoxicity changes for 90 products (median ×1.35, 34 products more
than doubled, asparagus ×40–119 from deltamethrin in the seedlings, soybeans and soybean oil and meal
BR ×20). Human toxicity changes for 45 (cancer) and 54 (non-cancer) products, mainly through the
heavy-metal layer, and all 37 negative toxicity results disappear. The `water` layer changes water use
for 95 products (median ×1.04; green beans CH ×0.035, palm MY ×0.08, sorghum CN ×0.09). Single score,
over the 96 products with positive scores: median +2.1 %, 25 products up more than 10 %, 5 more than
doubled, and down to ×0.56 where water use falls. These are the PEF method's default rules; under the
dairy PEFCR the pesticide layer would not apply.

**BAFU data issue.** Three greenhouse and tunnel products have negative water use in BAFU's own
published results: conventional strawberries from an unheated plastic tunnel, CH (−660 m³ world eq
per kg) and ES (−599,700), and conventional hors-sol tomatoes ES (−252). Water released exceeds water
withdrawn, by an implausible amount for the ES strawberries, and it dominates their single scores
(strawberries ES −2,449 mPt per kg). They are left out of the single-score statistics above. Organic
tunnel tomatoes CH are slightly negative only in our import with the global factor; with the Swiss
factor, as in BAFU's results, they are positive.

**ecoinvent re-allocation.** ecoinvent is licensed, so neither its data nor results derived from it
are published here. `ecoinvent_pef_realloc.py` writes them to `reports/`, which is not in the
repository.

## HESTIA database on BAFU (`hestia_db.py`)

`hestia-2026-04-30` holds the HESTIA aggregated release 2026-04-30 (CC-BY 4.0) as Brightway
datasets: 367 crop and processed-crop cycles (latest period per product and country, 48 countries),
684 datasets after allocation. It links only to `bafu-2026` and `ef-3.1-biosphere`, so the system
contains no ecoinvent.

| Part | Modelling |
|---|---|
| Products | All products with an economic value share; 191 multi-output cycles (grain + straw, oil + meal, rice milling) split by `split_multioutput`, rule `crop` (economic, LEAP) or `economic` for processing, with HESTIA's value shares. A co-product is dropped only if it is under 0.5 % of the main product's mass and under 2 % of the value (intercrop artefacts of the aggregation). |
| Inputs | Fertilisers per nutrient, diesel/petrol/LPG/fuel oil/natural gas as supply only (HESTIA carries the combustion emissions), grid electricity by country ("Electricity mix", GLO fallback), machinery, lime, tap water, processing chemicals → bafu-2026. Seed, and crop or crude-oil inputs of processing cycles → HESTIA datasets. Manure and compost: no burden (residues). |
| Emissions | HESTIA on-farm emissions → EF 3.1 (non-urban air, PM2.5 close to ground, fresh water). Fossil vs. biogenic vs. land-use-change CO₂/CH₄/CO split by emission source; CO₂ from biofuels is biogenic; dLUC, drained organic soils and vegetation burning count as land use change, carbon gains are not credited. Fuel-combustion totals that repeat per-fuel entries are dropped. P lost with eroded soil counts; N lost with eroded soil does not (N to water is nitrate leaching, as in PEF, ecoinvent and Agribalyse; scored as total N it would dominate marine eutrophication). HESTIA background emissions are dropped. |
| Pesticides | PEF default: 90 % of the applied active ingredient to agricultural soil, 9 % to air, 1 % to fresh water, matched by CAS; plus "Pesticide unspecified" production per kg AI. Active ingredients without an EF flow (unspecified, pesticide classes) are emitted as the cycle's identified mix where that carries at least as much mass, else as the mean mix over all cycles. Bulk, inert and processing substances HESTIA types as ingredients (oils, kaolin, chlorine, sulphur, acids) are left out. |
| Land, water | Occupation (site + unused duration) as Arable, Arable irrigated or Permanent crops; irrigation withdrawals as ground, river or lake water, taken as consumed. |

**Gaps** (`reports/hestia_unmapped.csv`): PM10/TSP/black carbon (no EF factor), coal, wood fuel,
sulphur, gypsum, branded NPK products, saplings and seedlings (counted in pieces). Unidentified
pesticide active ingredients are listed there too, marked "proxy": they are emitted through the
proxy above, not dropped.

**Data issues:** HESTIA's almond cycle for Spain reports 4 kg/ha deltamethrin, roughly 200 times a
normal dose, which dominates its freshwater ecotoxicity. In some aggregated cycles the straw amount
and its value share come from different sets of farms: pearl millet, India, lists 283 kg straw
against 2,135 kg dry matter of residue removed, yet gives straw 36 % of the value (similar for
sorghum, barley and wheat in India). The grain follows HESTIA's value shares; the per-kg impacts
of these straw datasets come out too high.

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
`EF-LCIAMethod_CF(EF-v3.1).xlsx`: all 25 indicators match all 276,704 global factors with no
deviation above 1 %. JRC also publishes about 300 country-specific factors per country for
acidification, eutrophication and water use, which the sentier method does not include.

### Normalisation and weighting factors

`ef31_nw.csv` holds the values of the JRC file `Normalisation_Weighting_Factors_EF_3.1.xlsx`
(global population 6,895,889,018). The same values, rounded, are in Tables 61 and 62 of the dairy
PEFCR. The two human toxicity normalisation factors were EF 3.0 values until 9 Nov 2026: cancer
1.69E-05 instead of 1.73E-05 CTUh, non-cancer 2.30E-04 instead of 1.29E-04 CTUh. The non-cancer error
had understated that category in the single score by 44 %.

### Reproducibility

All scripts were re-run from scratch on 9 Nov 2026. `run_demo.py`, `run_bafu.py`,
`ecoinvent_pef_realloc.py` and `bafu_workbook_check.py` reproduced their results exactly.
`bafu_pef_compare.py` had picked one of two duplicate EF flows at random for Mecoprop-P and Paraquat
(up to 0.4 % in human toxicity, non-cancer). The choice is now deterministic, and two runs agree to
machine precision (1e-16).

### Code review (9 Nov 2026)

An independent review, with each finding checked against the data, led to these corrections:

- Product selection: substring filters had dropped 31 food and feed products (oils and meals, sugar,
  molasses and pulps, greenhouse and tunnel strawberries and tomatoes, rapeseed, coconut) and kept a
  field operation. An explicit filter now gives 100 products.
- Pesticides: chloridazon and sodium methylarsonate were excluded by name prefixes, and silicon,
  chlorine, bromine and iodide were counted as active ingredients. BAFU's emissions to
  non-agricultural soil in agriculture datasets now go 90 % to agricultural soil.
- CFF: the electricity credit was in MJ against a per-kWh dataset (3.6 times too large).
- `run_bafu.py`: obsolete "xx" datasets for hay and PET landfill replaced; recycling is burden-free
  because the only PET recycling dataset is empty.
- Robustness: allocation scales uncertainty parameters with the amounts and no longer changes its
  input; the ecoinvent script checks properties and factors, and indexes the production entry by
  product row; heavy metals are netted per dataset and flow; summaries show sign changes separately.
- The ecoinvent results are unchanged by these corrections (identical files).

### Our BAFU baseline vs. BAFU's published EF 3.1 results

`bafu_workbook_check.py` compares the baseline of `bafu_pef_compare.py` with
`BAFU-2026 v1 LCIA Results_corrected.xlsx` for the 100 products (all present in the workbook; 99 with
a non-zero published value). The baseline includes one correction of the Brightway import (`fix`):
unspecified chromium emitted to soil and water is mapped to Cr(III), as BAFU's own EF mapping does
(731 datasets).

| Category | Within 1 % / 10 %, import as is | With the chromium fix | Median ours / published | Remaining cause |
|---|---|---|---|---|
| Climate change, acidification, eutrophication, PM, ozone, radiation, photochemical ozone, resources | 97–99 / 98–99 of 99 | unchanged | 1.000 | — |
| Human toxicity, non-cancer | 83 / 95 | unchanged | 1.000 | small |
| Human toxicity, cancer | 2 / 9 | **5 / 44** | 1.34 → **1.07** | Other flows (sunflower still ×8.6; soybean meal US ×0.12; sorghum CN sign differs). |
| Freshwater ecotoxicity | 43 / 69 | **59 / 75** | 1.00 → 0.99 | Flows characterised in BAFU's mapping but not in the import (rapeseed and rape oil ×0.08–0.12). |
| Land use | 63 / 76 | unchanged | 1.00 | Land transformation "To / From Arable" nearly cancel; ×5–6 for organic and extensive cereals. |
| Water use | 13 / 80 of 94 | unchanged | 1.05 | Turbine water and released water both characterised (±43 m³/m³) and nearly cancel; products BAFU regionalises (green beans CH ×29, strawberries tunnel CH ×32). |

Whey is left out: BAFU gives it no burden, so both values are about zero.

**Water use is not part of the fix.** Country factors for every dataset (the `water` layer) make the
products BAFU regionalises match exactly (asparagus PE and MX, green beans KE, the greenhouse and
tunnel strawberries and tomatoes), but lower agreement overall from 80 to 33 of 94 products within
10 % (palm MY, sorghum CN drop to ×0.08). BAFU regionalises water only in some datasets, which the
Brightway import does not show. Country-level water use is a PEF requirement, so it is a PEF layer
instead.

**What it means.** The baselines behind the main PEF findings match: asparagus (CH, MX, PE) and
soybeans BR within 0.2 % for freshwater ecotoxicity. The BAFU and PEF versions use the same import,
so the PEF ratios hold. Land use for organic and extensive cereals and water use for some vegetables
still differ from BAFU's published numbers.

**Still open in the import** (would also apply to the sentier loader): turbine water, the remaining
ecotoxicity flows (rapeseed, rape oil), and land transformation. In BAFU's data itself: the negative
water use of three greenhouse and tunnel products (see Main results).

## Open points

- Implement the PEF N field model.
- Enter PEF default values: slaughterhouse prices and mass fractions; CFF A, R2, R3, Qs/Qp from
  Part C of Annex II.
- Full heavy-metal rule (Annex I 4.4.1.6 (a)): emissions equal to field inputs without crop uptake. This
  needs the metal inputs, which BAFU's net balances do not show.
- Decide per product which PEFCR applies; it can override the PEF defaults (e.g. pesticides in the
  dairy PEFCR).
- Characterise ecoinvent with the same EF 3.1 factors as BAFU (e.g. a flowmapper bridge) so both
  databases can be compared directly.

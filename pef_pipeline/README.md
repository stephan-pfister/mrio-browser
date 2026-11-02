# PEF pipeline (Brightway 2.5 + wurst)

A small pipeline that applies PEF (EF Recommendation (EU) 2021/2279) multifunctionality rules
and the Circular Footprint Formula to a Brightway foreground database. It uses wurst for the data step.

```
foreground primary data (multi-output, wurst dicts)
   │  pef_allocation.split_multioutput   rule per dataset: dairy_farm | slaughterhouse | crop | economic | physical
   │  pef_allocation.cff_* + add_cff_exchanges   packaging / EoL as plain exchanges
   ▼
wurst.write_brightway2_database  ──►  single-output foreground DB, linked to background
   ▼
bw2calc LCA with EF 3.1 (16 categories)  ──►  single_score() with EF normalisation & weighting
```

## Rules implemented

| rule | used for | formula |
|---|---|---|
| `dairy_farm` | milk vs. live animals at farm gate | AF_milk = 1 − 6.04 · M_meat / M_milk (FPCM, IDF 2015) |
| `slaughterhouse` | meat, offal, hides, bones, Cat 1/2/3 | economic; enter the PEF default mass fractions and prices (demo values are illustrative) |
| `crop` | grain vs. straw / residues | economic (LEAP) |
| `physical` | e.g. dairy plant (dry matter), mass, energy | property named in `pef_property` |
| `economic` | fallback (hierarchy step 3) | `price` property |

CFF: `cff_material` and `cff_energy_disposal` return `(dataset key, coefficient)` pairs. The CFF is
linear in the E-terms, so it becomes ordinary technosphere exchanges. Take A, B, R1–R3 and Qs/Qp
from Annex C.

## Run the demo

```bash
python run_demo.py
```

The demo creates its own project `pef_pipeline_demo` with a toy background and a GWP-only method. All
numbers are illustrative.

## Using real data

- Background: EF 3.1-compliant datasets where you have them; otherwise ecoinvent.
  Ecoinvent *cut-off* is already allocated with ecoinvent rules. Re-allocating background processes
  the PEF way needs the unallocated ("undefined") system model. In practice, apply PEF rules
  to the foreground and to the key agricultural processes you copy from the background.
- Method: EF 3.1 in Brightway, with its biosphere matched to the background database.
  Normalisation and weighting factors come from the EF reference package.
- Alternative to the wurst split: `bw-functional` (`multifunctional`) stores the co-products and does
  property-based allocation inside Brightway. It cannot do the dairy IDF formula or CFF out of the
  box, and it is still beta.

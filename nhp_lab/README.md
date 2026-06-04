# NHP lab customizations to phy 2.0b5

Customizations for the Akers-Campbell lab NHP Neuropixels spike-sorting workflow
(Kilosort4 sorts, QC'd by bombcell, curated in phy). Branch: `nhp-customizations`.

Base: phy 2.0b5. The three modified package files on this branch were verified
identical to the pristine 2.0b5 `master` before editing, so the git diff shows
ONLY the lab customizations.

## 1. Modified phy package files (live in site-packages)
Each maps to the conda env install path
`...\envs\phy2\Lib\site-packages\phy\`. A `.bak` of the pristine file sits next
to each in that install.

| File in this repo            | Customization                                                                 |
|------------------------------|-------------------------------------------------------------------------------|
| `phy/utils/color.py`         | `_make_cluster_group_colormap()`: +index 4 orange (nonsoma_single), +index 5 purple (nonsoma_multi). Feeds the `cluster_group` PLOT color scheme. |
| `phy/apps/base.py`           | `_add_default_color_schemes` `group_colors`: + `nonsoma_single:4`, `nonsoma_multi:5`. |
| `phy/cluster/supervisor.py`  | `_CLUSTER_VIEW_STYLES` (~L260): cluster TABLE row colors — good #86D16D, mua #3366E6, noise #D92626, nonsoma_single #FF8C1A, nonsoma_multi #9933CC. |

## 2. Plugins + config (live in `~/.phy/`)
| File in this repo                     | Installs to                          |
|---------------------------------------|--------------------------------------|
| `nhp_lab/plugins/car_filter.py`       | `~/.phy/plugins/car_filter.py`       |
| `nhp_lab/plugins/pair_correlogram.py` | `~/.phy/plugins/pair_correlogram.py` |
| `nhp_lab/phy_config.py`               | `~/.phy/phy_config.py`               |

- **car_filter.py** (CARFilterPlugin): alt+r cycles raw/high_pass/car/car_highpass.
  NOTE phy CAR = flat global CAR (median across displayed channels), NOT ADC-group
  demux. Demux happens upstream in CatGT (-gbldmx).
- **pair_correlogram.py** (PairCCGPlugin): select 2 clusters -> ACG_A | rawCCG +
  jitter-overlay | ACG_B triptych. Jitter flanks distorted at default +/-25ms
  (widen with `cw 100`).

## 3. The 5-category group scheme (data, lives in the sort dir, not here)
Two TSVs in each sort dir carry the bombcell-derived 5-category scheme:
- `cluster_group.tsv` — drives COLOR via the CSS. Values use UNDERSCORE:
  `good, mua, noise, nonsoma_single, nonsoma_multi` (must match the
  `tr[data-group=...]` selectors in supervisor.py).
- `cluster_bc_unitType.tsv` — the readable TEXT column. Values use HYPHEN:
  `somatic, mua, noise, nonsoma-single, nonsoma-multi`.

GOTCHA (fixed 2026-06-04): phy's `_load_metadata` globs EVERY `*.tsv` and keys
metadata on the HEADER column name, not the filename. A backup file
(`cluster_bc_unitType_orig.tsv`) that kept the header `bc_unitType` silently
overwrote edited values (whichever file globs last wins). Fix: rename the backup's
header column (e.g. to `bc_unitType_orig`) so it never collides. Do NOT leave two
TSVs sharing a column header in a sort dir.

## Adding a 6th category later
Update BOTH sort-dir TSVs + `color.py` colormap + `base.py` group_colors + the
`supervisor.py` CSS. Unknown groups silently render red.

# DPOAE Growth Function fitting Toolkit

Distortion-product otoacoustic emission (DPOAE) growth-function pipeline for DP/GR data exported from **Otodynamics Echoport 292 with ILOv6**. The workflow also works with any dataset that preserves the same column layout (e.g., `Freq (Hz)`, `F1 (dB)`, `F2 (dB)`, `DP (dB)`, `Noise±sd (dB)`, `2F2-F1 (dB)`, etc.).

## What the script does
- Discovers participant folders under `Participants_170725`, matching `####_(DP|GR)_*.csv` per ID/ear.
- Reads DP-gram (`DP`) and growth-rate (`GR`) CSVs into structured dictionaries.
- Normalizes GR frequencies to canonical 1414/4243 Hz, removes duplicate entries, and keeps the best average-DP record.
- Filters growth data for sufficient signal-to-noise, creates focused L2 windows (e.g., 35–55 dB), and stores those slices.
- Runs multiple fittings:
  - Linear regression between 35–55 dB SPL.
  - Cubic fit (with derivative slope between 40–60 dB) and quadratic derivative for curvature analysis.
  - Piecewise linear fits (PWLF) in Pascals with breakpoints at `[min(L2), 25, knee, max(L2)]` where `knee ∈ {45, 50, 55}`.
- Extracts slopes, intercepts, quadratic minima, PWLF slopes, and R² metrics, aggregating everything into `_OaeIO_allIDFitParams.xlsx`.
- Produces multi-panel plots of DP I/O in both dB SPL and Pascals, saving SVGs in `Figures/`.

## Expected input format
Each CSV must include the Otodynamics headers:

| Column            | Description                                   |
|-------------------|-----------------------------------------------|
| `Freq (Hz)`       | Nominal f0 (≈1414 or 4243 Hz for GR files).   |
| `F1 (dB)`, `F2 (dB)` | Primary stimulus levels.                  |
| `DP (dB)`         | DPOAE level.                                  |
| `Noise+1sd/2sd (dB)` | Noise floor estimates.                     |
| `2F2-F1`, `3F1-2F2`, … | Optional distortion products.           |

GR files are expected to include all L2 steps (65 → 20 dB). The script rewrites the `F2 (dB)` vector to `[65, 60, …, 20]` if necessary.

## Key outputs
- `_OaeIO_allIDFitParams.xlsx`: participant-level table of slopes, intercepts, quadratic minima, and PWLF stats.
- `Figures/<ID>_growth_Pa.svg`: DP growth curves plotted in Pascals.
- Console logs highlighting duplicate GR entries that were cleaned.

## Dependencies
Python 3 with `numpy`, `pandas`, `matplotlib`, and `pwlf`. The script assumes folder access to `Participants_170725/` and writes results into `Figures/` plus the Excel file in the project root.

## Usage
Run the script directly (e.g., `python extract_DPOAE_fixed_togithub.py`). Ensure your CSV exports reside in the mirrored folder hierarchy. Review the generated Excel metrics and SVG figures for further analysis.

# Local dashboard

From the repository root, start a loopback server and open [the local dashboard](http://127.0.0.1:8734/dashboard/) in your browser.

```sh
python3 -m http.server 8734 --bind 127.0.0.1
```

The page uses the included `data.json` and local assets. No network services or JavaScript packages are needed. Serving the repository root keeps links to the technical reports available.

## Regenerate the data

From the repository root, create the research environment and regenerate the results:

```sh
python3.13 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-research.txt
python -m pip install --no-build-isolation -e .
python experiments/analyze.py
python experiments/sensitivity.py
python experiments/dashboard.py
python -m unittest discover -s tests -v
```

The exporter reads the published CSVs and calculates the input swaps, the Student-t floor curves, the feasible GEMs scenarios and the concentration comparison. It verifies shared results against the CSVs before writing the JSON. The chart connects precomputed floor values and marks a crossing found by root solving. All selectable values come from evaluated cases. A deterministic subset of the three-year block bootstrap draws is exported for the scatter. Reported intervals and correlation use the full bootstrap. For each draw the exporter computes the excess pair probability over PD squared and solves the crossing at that draw's PD and pair target. Draws at or below independence have no finite crossing and remain in the percentile calculation. The band shows uncertainty in the crossing location, measured in degrees of freedom. It does not represent uncertainty in the pair-probability curve.

The single-name tail share counts defaults by the large name above 99.9% VaR and takes a proportional slice of the probability mass at VaR to complete the worst 0.1% of losses.

The baseline uses the supplied simulation histograms. To regenerate them, run `python experiments/simulate.py --workers 2` before `analyze.py`. The persistence, validation and borrower-mix studies have separate scripts in [experiments](../experiments/).

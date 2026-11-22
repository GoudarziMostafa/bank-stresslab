# What default histories can tell a risk manager about tail losses

Two credit models give every borrower the same one-year default probability and every pair the same joint default probability. They agree on expected loss and loss volatility. Yet for 1,000 equal exposures, the Gaussian model estimates **10.48% of total exposure** as the average loss in the worst 0.1% of years. The matched Student‑t model estimates **9.41%**. The Expected Shortfall (ES) gap is **107 basis points**. The illustrative capital proxy, 99.9% VaR minus expected loss, differs by **54 bp**. [Baseline results](../results/reference.csv)

## Where the difference comes from

Given that two specified borrowers default, a third defaults with probability 9.8% under Gaussian and 8.7% under Student‑t. Matching pairs leaves room for different patterns of large default clusters. Giving one borrower 10% of total exposure raises severe losses under both models while narrowing the ES gap to 83 bp. An independent simulation of 6.4 million scenarios estimated the signed difference at −114 bp, with a 95% Monte Carlo interval of −123 to −105 bp. The headline estimates use finite-portfolio integration. [Calibration](../results/calibration.csv) · [Simulation uncertainty](../results/monte_carlo.csv)

Across the sensitivity grid the matched gap ranges from 2 to 304 bp in magnitude, and 14 of 30 combinations cannot match the default pairs. At GEMs-implied PD, LGD and Gaussian correlation, the ten-degree Student‑t model cannot match the estimated pair-default rate. A feasible comparison using 100 degrees of freedom gives a gap of about 3 bp. The 107 bp synthetic headline depends on a different set of assumptions. [Parameter grid and GEMs calibration](sensitivity.md)

## What changes a decision

| Risk management question | Finding |
|---|---|
| How should a bad year change the forecast? | After 43 or more defaults, next-year 99.9% ES stays at 10.48% with independent years and rises to 16.9% at illustrative annual factor persistence 0.6. Both models have the same unconditional annual loss distribution. [Forecast study](persistence_design.md) |
| Would validation detect the difference? | Over 20 years, a lag-one test missed persistence in 52% of simulated histories. A likelihood test that knew both candidate models missed it in 23%. Their false-alarm rates were 5.0% and 5.6%. [Validation study](validation_power.md) |
| What would trigger a 15% ES limit? | At persistence 0.6, the exact-count threshold is 64 defaults. A binary report of “23 or more” also triggers the limit. In 81% of those alerts the exact-count forecast is below the limit. [Reverse stress test](reverse_stress.md) |
| Could a changing loan book look like a persistent cycle? | In a separate two-grade portfolio where only borrower mix persisted, totals alone falsely flagged a cycle in 680 of 8,000 histories. Observing starting grade sizes and grade defaults reduced that to 393. [Attribution study](borrower_mix.md) |

For a review of clustered defaults, retain the exact counts and reconstruct who was on the book before the year began. Both pieces of information changed model decisions in these experiments. The validation experiments use illustrative portfolios and parameters. Selecting a model or setting a capital limit for a lender requires evidence about its actual exposures, defaults and recoveries. [Reproduction instructions](../README.md#reproduce)

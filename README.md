# Bank StressLab

[Open the live dashboard](https://goudarzimostafa.github.io/bank-stresslab/)

If two credit models agree on how often borrowers default alone and in pairs, how far apart can their worst-loss estimates be?

In a synthetic book of 1,000 equal loans, Gaussian 99.9% expected shortfall is **10.48% of exposure**. Matched Student-t gives **9.41%**, a gap of **107 basis points**. With 30 years of GEMs private-lending data, the ten-degree Student-t model cannot match the estimated clustering. An illustrative 100-degree match gives a gap of **3.3 bp**. The two comparisons use different inputs.

Try **Feasibility** to see where a match becomes possible, then **Concentration** to see what happens when one borrower accounts for 10% of exposure. These are model comparisons on hypothetical portfolios.

[Dashboard screenshot](figures/dashboard-preview.png) · [Methods and results](reports/sensitivity.md) · [Run locally and reproduce](dashboard/README.md)

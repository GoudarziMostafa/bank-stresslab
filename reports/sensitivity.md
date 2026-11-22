# Does the tail gap survive a change of assumptions?

## Holding latent correlation fixed or matching default pairs

Consider 1,000 equal exposures with a 1% annual probability of default (PD) and 45% loss given default (LGD). Expected shortfall (ES) averages losses in the worst 0.1% of years. The choice of calibration changes even the direction of the comparison.

| Model | Latent correlation | Probability two specified borrowers default | 99.9% ES as % of exposure |
|---|---:|---:|---:|
| Gaussian benchmark | 0.25 | 0.044% | 10.48% |
| Student‑t, 10 degrees of freedom, same latent correlation | 0.25 | 0.093% | 18.48% |
| Student‑t, 10 degrees of freedom, matched default pairs | 0.0395 | 0.044% | 9.41% |

A small shared scale draw pushes the Student‑t borrowers toward extreme outcomes together. At equal latent correlation, it more than doubles the probability of a default pair. Refitting the Student‑t correlation removes enough common Gaussian exposure to reverse the ES ordering at this horizon and cutoff. Pair matching fixes the first two moments of the conditional PD. It leaves the distribution of larger default clusters free. For example, after two specified defaults the probability of a third is 9.8% under Gaussian and 8.7% under matched Student‑t. This helps interpret the computed ES ordering without implying a general theorem. [Original design and equations](research_design.md) · [Third-default calculation](../results/calibration.csv)

## Sensitivity across 30 combinations

Each entry gives Student‑t minus Gaussian ES in whole basis points of total exposure after refitting to the same pair-default probability. A dash means the pair target is below the Student‑t floor at zero latent correlation. The model permits nonnegative latent correlations. Infeasibility is a property of that model and target.

| Student‑t degrees of freedom | Gaussian 0.10 | Basel 0.193 | Baseline 0.25 | 0.35 | 0.50 |
|---:|---:|---:|---:|---:|---:|
| 4  | — | — | — | — | −304 |
| 5  | — | — | — | — | −186 |
| 8  | — | — | — | −116 | −68 |
| 10 | — | — | −107 | −71 | −43 |
| 20 | — | −31 | −24 | −17 | −10 |
| 50 | −8 | −5 | −4 | −3 | −2 |

The matched gap is negative in all 16 feasible cells, ranging from about 2 to 304 bp in magnitude. At the same latent correlation the Student‑t ES is higher in all 30 cells. At four degrees of freedom, a Student‑t model with zero latent correlation already gives a 0.095% pair-default probability, above the baseline target of 0.044%. It becomes feasible at the 0.50 Gaussian benchmark. These results cover the stated PD, allocation, confidence level and grid.

The experiment checks total probability, marginal PD and pair-default probability for each integrated model. Successive quadrature orders must agree on the value-at-risk (VaR) cutoff and on ES within 0.05 bp. The independent borrower simulation in the original study checks the baseline. [All targets, floors, fitted correlations and both ES comparisons](../results/dependence_grid.csv) · [Simulation check](matched_dependence.md)

## GEMs private lending

The 2024 GEMs publication gives active counterparties and defaults by calendar year for private lending by multilateral development banks and development finance institutions. The 1994–2023 series contains 1,996 defaults in 55,993 counterparty-years. Their ratio reproduces the published 3.56% rate. For a model of a randomly selected **year**, we give each year's default rate equal weight, yielding a 3.80% PD. Larger, later cohorts get more weight in the published pooled rate. Table 2 gives 72.2% average recovery on 1,561 defaulted private contracts. We use its complement, 27.8%, as fixed LGD throughout. [GEMs report, Tables 2 and 5](https://www.eib.org/files/publications/20240217_default_and_recovery_statistics_private_and_public_lending_1994_2023_en.pdf) · [Annual counts used here](../data/gems_private_1994_2023.csv)

For year $t$ with $d_t$ defaults among $n_t$ active counterparties, the observed rate is $r_t=d_t/n_t$. We average $d_t(d_t-1)/[n_t(n_t-1)]$ across years to estimate the chance of two **distinct** defaults, avoiding the contribution from pairing a borrower with itself. We fit the one-factor Gaussian latent correlation $\rho$ to this pair estimate at the calendar-year PD $p$. To estimate annual factor persistence $\phi$, we use the average adjacent-year product $(r_t-p)(r_{t+1}-p)$. In the stationary model with fresh borrower shocks each year, its population counterpart is the Gaussian pair probability at latent correlation $\rho\phi$, less $p^2$.

| Calendar years | Equal-year PD | Gaussian latent correlation | Factor persistence | Gaussian 99.9% ES | Matched t, 10 df | Matched t, 100 df |
|---|---:|---:|---:|---:|---|---:|
| 1994–2023 | 3.80% | 0.027 | 0.40 | 3.13% | Infeasible | 3.09% |
| 1994–2008 | 4.26% | 0.034 | 0.38 | 3.81% | Infeasible | 3.78% |
| 2009–2023 | 3.34% | 0.007 | No nonnegative fit | 1.90% | Infeasible | Infeasible |

The ES figures apply to a hypothetical book of 1,000 equal exposures with those inputs. The 30-year pair-default estimate is 0.163%, below the ten-degree Student‑t floor of 0.262%. A 100-degree model can match it. At 3.80% PD, 27.8% LGD and Gaussian latent correlation 0.027, the ES gap is about 3 bp. The 107 bp synthetic result uses 1% PD, 45% LGD, Gaussian correlation 0.25 and ten degrees of freedom. The GEMs comparison changes all four inputs, so it describes another calibration without attributing the difference to a single parameter. We chose 100 degrees to illustrate a lighter-tailed alternative. The counts do not identify it as the best fit.

The approximate 95% interval for the annual pair probability is **0.107% to 0.234%**. For the feasibility question the relevant quantity is the ten-degree Student‑t floor minus that pair probability. Its point estimate is **+0.098 percentage points**, with a 95% interval of **+0.089 to +0.111 percentage points**. In each of 12,000 bootstrap draws we resample the 30 years in circular three-year blocks with their counts and denominators together, re-estimate PD, recompute the floor at that PD, and subtract the draw's pair estimate. There are ten block starts per three-year resample. Resamples with more high-default years tend to raise both the pair estimate and the Student‑t floor through the re-estimated PD. This positive covariance reduces the variance of their difference. Percentile intervals using one- and five-year blocks also exclude zero. The intervals describe variation under resampling this historical sequence. Thirty years, changing membership and possible structural breaks limit coverage for a different book or period. [Bootstrap results](../results/gems_uncertainty.csv)

The fitted GEMs pair probability is about **0.019 percentage points above independence**, with an approximate paired 95% bootstrap interval of **0.0023 to 0.0329 percentage points**. Solving separately in each draw puts the zero-correlation Student‑t crossing at **ν ≈ 58** in the full sample, with an approximate 95% bootstrap interval of **38 to 389**. A draw at or below independence would have no finite crossing. All draws in this three-year-block resample had a finite crossing. The wide interval reflects how gradually the floor approaches independence. The floor at ν = 10 exceeds the target in all these draws. [Dashboard setup and calculation](../dashboard/README.md)

The two 15-year fits show how uncertain the interpretation is. The later sample has negative lag covariance, outside the dynamic model's nonnegative persistence range. Annual portfolio growth, changes in borrower mix and repeated counterparties can all affect these moments. The fixed LGD also hides variation across contracts and downturns. These are in-sample, model-implied loss estimates. Thirty yearly observations offer little direct evidence about a one-in-a-thousand-year tail. [Empirical fit, feasibility and unrounded outputs](../results/gems_calibration.csv) · [Composition diagnostic](borrower_mix.md)

## What this means for capital planning

Basel's asymptotic single risk factor (ASRF) model assumes enough small exposures for borrower-specific risk to diversify away. At 1% PD, the unadjusted Basel corporate asset-correlation curve gives 0.1928. The grid uses the curve before any firm-size or financial-institution adjustment. The foundation IRB rules assign 40% LGD to senior unsecured claims on other corporates and 45% to the corresponding claims on sovereigns, banks and other financial institutions specified in CRE32.6. The experiment's 45% is an illustrative common severity. [Basel risk-weight functions, CRE31](https://www.bis.org/basel_framework/chapter/CRE/31.htm?inforce=20230101&published=20200327) · [Basel risk components, CRE32.6](https://www.bis.org/basel_framework/chapter/CRE/32.htm?inforce=20230101&published=20200327)

At the original Gaussian inputs (latent correlation 0.25), the infinitely granular 99.9% loss quantile is about 826 bp. It follows from multiplying LGD by the conditional PD in the 0.1% adverse factor state. Finite-portfolio VaR is about 833 bp for 1,000 equal names and 936 bp when one name carries 10% of exposure. A granularity adjustment approximates the extra capital associated with residual borrower-specific risk. Here the finite-portfolio calculation measures the quantiles directly. This comparison is distinct from the matched-model ES gap. [Finite-portfolio results](../results/reference.csv) · [Basel IRB explanatory note](https://www.bis.org/bcbs/irbriskweight.htm) · [Basel study of credit concentration](https://www.bis.org/publ/bcbs_wp15.htm)

For an MDB single-name limit review, moving 10% of exposure to one obligor increased Gaussian 99.9% VaR by about 104 bp versus equal weights, holding total exposure and borrower risk parameters fixed. A risk committee can use that incremental unexpected loss to compare a proposed large commitment with the available economic-capital buffer or a syndication alternative. The granularity adjustment addresses the borrower-specific risk that survives at large weights. Country and sector shocks, guarantees and project recoveries need separate assumptions before translating this experiment into a limit.

In a capital review, rerun the challenger with the bank's exposure weights, borrower grades and recovery assumptions. Compare one-year VaR less expected loss, a common economic-capital proxy, and use ES to show how severe losses beyond VaR can be. At baseline the two models differ by 54 bp on the proxy and 107 bp on ES. Booking the ES difference as a capital add-on would require a separate judgment about model fit, the institution's capital policy and its risk appetite. Regulatory IRB capital has its own risk weights, maturity treatment and input rules.

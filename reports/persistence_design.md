# Forecasting after a high-default year

## Design fixed before calculation

Use the granular Gaussian benchmark with 1,000 equal exposures, one-year PD 1%, asset correlation 0.25 and fixed LGD 45%. Each year starts with the same portfolio composition. Replace exited borrowers with equivalent exposures and draw new idiosyncratic shocks. Only the common credit environment can transmit information across years.

Let the common risk factor follow a stationary Gaussian process,

$$
Z_0\sim N(0,1),\qquad Z_1=\phi Z_0+\sqrt{1-\phi^2}\,\eta,\qquad \eta\sim N(0,1).
$$

Independently of the past given the contemporaneous factor, the annual count has law $K_t\mid Z_t=z\sim\mathrm{Binomial}(n,u(z))$, with $u(z)=\Phi((\Phi^{-1}(p)-\sqrt{\rho}\,z)/\sqrt{1-\rho})$. Loss as a share of total exposure is $gK_t/n$. Lower factor values mean worse credit conditions. Every choice of persistence has exactly the same stationary one-year count and loss distributions. The forecast starts from the stationary prior and uses one completed year of counts.

**Primary contrast.** Find the smallest count $v$ with $P(K_0\le v)\ge0.95$. Observe the event $A=\{K_0>v\}$ and report its actual probability. Discreteness can make this probability smaller than 5%. Measure the change in next-year 99.9% ES from its unconditional value in basis points of total exposure. Use illustrative annual persistence values 0, 0.3, 0.6 and 0.9. The main comparison is 0.6 versus 0. Report conditional EL and VaR alongside ES.

**Exact-count sensitivity.** Also forecast given $K_0=v+1$, the mildest count qualifying as stress. Use its binomial point probability as the likelihood. The primary forecast pools all qualifying counts before taking ES. The sensitivity forecast uses a single count.

**Numerical checks.** Recover the published granular Gaussian baseline. Cross-check every persistence setting's full unconditional count distribution against that baseline. Integrate event and exact-count forecasts at two resolutions, require ES changes below 0.1 bp and stable VaR cutoffs. Combine forecasts after stress and its complement to recover the unconditional count probabilities. Independently simulate pairs of years for the primary comparison with 32 independent batches of 400,000 paths, fixed seed 20261121, and report a paired 95% Monte Carlo interval for the ES revision. That interval measures numerical uncertainty at fixed parameters.

**Interpretation.** Duffie, Eckner, Horel and Saita motivate inference about an unobserved common credit condition. Their estimated model uses firm-specific default intensities, covariates and a filtered frailty process. Here the annual factor persistence is a sensitivity parameter. Stationary portfolio composition, constant PD and LGD, and fresh borrower shocks isolate the forecast effect of persistence. [Duffie et al. (2009), Figures 4–5](https://web.stanford.edu/~duffie/fcdt.pdf)

## Results

The stationary 95th percentile is 42 defaults. A report of **43 or more** occurs with probability **4.8611%**. A report of **exactly 43** occurs with probability **0.1884%**. Before receiving either observation, 99.9% ES is **10.4825% of total exposure** in every model.

| Annual persistence | EL after 43 or more | 99.9% VaR after 43 or more | 99.9% ES after 43 or more | ES revision (bp) | 99.9% ES after exactly 43 |
|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.45% | 8.33% | 10.48% | 0.00 | 10.48% |
| 0.3 | 0.94% | 12.24% | 14.69% | +420.41 | 13.34% |
| 0.6 | 1.70% | 14.58% | 16.94% | +645.56 | 13.21% |
| 0.9 | 2.80% | 14.99% | 17.27% | +679.02 | 7.82% |

All percentages in the table are shares of total EAD. At persistence 0.6, the expected number of defaults next year increases from 10 to about 38 after the stress trigger. The 645.56 bp ES revision measures the effect of the observation within that model. It also equals the difference from the independent-year forecast after the same observation.

The exact-count column illustrates how a loss forecast can become more pessimistic about the average and less pessimistic about the extreme tail. At persistence 0.9, observing exactly 43 defaults raises expected loss from 0.45% to **1.72%**, but lowers 99.9% ES from 10.48% to **7.82%**. Under the forecast the next factor has mean −1.48 and standard deviation 0.45, versus zero and one before observing the count. The exact observation provides evidence against extremely adverse current states. High persistence carries that information into next year. A report giving only a lower bound includes years with many more than 43 defaults and produces a different predictive distribution. [Full forecast results](../results/persistence_forecasts.csv)

## How the results were checked

For an observation with likelihood $h(z)$, the filtered density of the current factor is proportional to $h(z)\varphi(z)$, where $\varphi$ is the standard normal density. The calculation integrates tomorrow's normal transition density against this posterior and mixes the finite binomial count law over the resulting predictive density. For the stress observation, $h(z)=P(K_0>42\mid Z_0=z)$. For the exact count, $h(z)=P(K_0=43\mid Z_0=z)$. Discrete ES assigns only the requisite probability mass to the VaR cutoff.

Composite Gauss–Legendre integration uses 48 panels on $[-12,12]$ and compares orders 16 and 32 per panel. All VaR cutoffs were stable and all ES changes met the 0.1 bp threshold. For each persistence value, the full unconditional count distribution and the mixture of the stress and complementary forecasts matched the original Gaussian benchmark to within $10^{-11}$ absolute probability per count. A separate one-borrower test compares the forecast conditional on default to a bivariate normal calculation with cross-year latent correlation $\rho\phi$.

In an independent simulation, **621,859 of 12.8 million** two-year paths met the stress trigger. The paired ES revision for persistence 0.6 versus zero was **642.74 bp**. Its approximate 95% Monte Carlo interval was **618.26 to 667.23 bp**, containing the quadrature estimate of 645.56 bp. The simulation draws both years' counts directly conditional on their factors. The persistence settings share the initial histories and future common-factor innovations. Interval estimates use a delete-one-batch jackknife over 32 independent batches, centered on the pooled estimate.

After seeing the exact-count reversal, we added a simulation check for persistence 0.9. Among **24,029** qualifying paths, the estimate was **7.98%** with an approximate 95% interval of **7.55% to 8.40%**, compared with **7.82%** by integration. Only about 24 simulated observations contribute to that 0.1% tail, so the interval is a rough check. The central results use quadrature. [Simulation output](../results/persistence_simulation.csv) · [Reproduction instructions](../README.md#reproduce)

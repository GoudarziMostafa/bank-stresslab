# Bank StressLab research design

The design and calibration were recorded on 20 November 2026 before estimating portfolio tail losses. Realized precision and outcomes appear in the [executive summary](executive_summary.md).

## Question and prior work

How far apart are two models' one-year credit loss tails when every borrower's default probability and every pair's joint default probability agree? Does reallocating exposure to one large borrower widen or narrow that gap?

Frey and McNeil already study credit loss models calibrated to common default probabilities and default correlations. This project builds on that experiment with a finite portfolio, a change in single-name concentration, and two numerical methods for evaluating the same tail quantities. Synthetic data permit a controlled model comparison. They cannot tell us which dependence law describes real borrowers. [Frey and McNeil (2002)](https://statmath.wu.ac.at/~frey/publications/rome.pdf) and [Frey and McNeil (2003), Example 4.8](https://www.wu.ac.at/fileadmin/wu/d/i/statmath/Publications/Publications_Frey/Journal_Papers/frey-mcneil-JoR.pdf)

## Portfolio and estimands

Use 1,000 distinct borrowers, a one-year marginal default probability $p=0.01$, constant LGD $g=0.45$, and total EAD normalized to one. All inputs are illustrative. For weights $w_i$ summing to one and default indicators $D_i$, the loss as a fraction of total EAD is $L=g\sum_i w_iD_i$.

| Allocation | Weight of designated borrower | Weight of each other borrower |
|---|---:|---:|
| Granular | 1/1,000 | 1/1,000 |
| Concentrated | 0.10 | 0.90/999 |

The primary estimand is $10^4(\mathrm{ES}^{t}_{0.999}-\mathrm{ES}^{G}_{0.999})$ for the granular allocation, in basis points of total EAD. The concentration comparison is the same signed gap under the concentrated allocation minus the granular gap. Either sign is admissible. Report both unrounded gaps before calculating their difference. Secondary outputs are VaR and ES at 99% and 99.9%, expected loss, the illustrative capital proxy VaR minus expected loss, and $P(L>\ell)$ versus a common loss threshold $\ell$. One basis point of total EAD is 0.0001 in these units.

## Dependence and calibration

Let $Z,\epsilon_1,\ldots,\epsilon_N$ be independent standard normals and $W\sim\chi^2_{\nu}$ independent of them. Define $Y_i(r)=\sqrt{r}Z+\sqrt{1-r}\epsilon_i$. The Gaussian model defaults when $Y_i(r_G)\leq\Phi^{-1}(p)$. The Student model defaults when $Y_i(r_T)/\sqrt{W/\nu}\leq t_{\nu}^{-1}(p)$. The scale draw is common to all borrowers. For the chosen $\nu>2$, $r_T$ is the correlation of the latent t variables. The default indicators have a different correlation.

Choose $r_G=0.25$ and $\nu=10$. Write $C^G_r$ and $C^t_{r,\nu}$ for the corresponding bivariate copulas. Set $q=C^G_{r_G}(p,p)$ and solve $C^t_{r_T,\nu}(p,p)=q$ for $r_T\geq0$. At these inputs the Student pair probability ranges from 0.000372870040 at $r_T=0$ to 0.000930306883 at $r_T=0.25$, bracketing the target. Conditional on the common scale, the bivariate normal lower-tail probability increases strictly with correlation. Averaging over that scale preserves monotonicity, so the solution is unique.

| Quantity | Calibration target or solution |
|---|---:|
| Joint probability of two specified defaults, $q$ | 0.000437515126 |
| Common default-indicator correlation, $(q-p^2)/[p(1-p)]$ | 0.034092437 |
| Student latent correlation, $r_T$ | 0.039492988 |

Compute the root at full precision in the experiment. A t model with five degrees of freedom has a floor of 0.000746693607 at this PD, so cannot match this benchmark with a nonnegative one-factor loading parameter. Degrees of freedom were selected for feasibility before estimating any portfolio tail. There is no empirical estimate of degrees of freedom in this design.

Two deterministic calculations agreed on $q$ to within $6\times10^{-14}$ and on $r_T$ to within $4\times10^{-11}$. The first integrates the bivariate law by conditioning on one latent variable. The second exploits radial symmetry. For a common threshold $h<0$ on the latent scale and latent correlation $r$, integrate $S_R(|h|/\cos\theta)/\pi$ over $[\arccos(r)/2,\pi/2]$. The radial survival function is $S_R(x)=\exp(-x^2/2)$ for the Gaussian and $(1+x^2/\nu)^{-\nu/2}$ for the Student model. The lower limit arises from the angle between the two default half-planes after whitening.

## Mathematical checks

Write $H=\sum_i w_i^2$. In either model,

$$
E[L]=gp,\qquad \operatorname{Var}(L)=g^2\{q-p^2+(p-q)H\}.
$$

These identities follow by splitting the second moment into single-borrower and distinct-pair terms. Expected loss equals 45 bp in both allocations. Loss standard deviation equals 83.8351 bp for granular weights and 94.4886 bp for concentrated weights, in each model. These are analytical test targets.

With $h_G=\Phi^{-1}(p)$ and $h_T=t_\nu^{-1}(p)$, the conditional default probabilities are

$$
u_G(z)=\Phi\!\left(\frac{h_G-\sqrt{r_G}z}{\sqrt{1-r_G}}\right),\qquad
u_T(z,w)=\Phi\!\left(\frac{h_T\sqrt{w/\nu}-\sqrt{r_T}z}{\sqrt{1-r_T}}\right).
$$

For a model's conditional probability $U$, check $E[U]=p$ and $E[U^2]=q$. As an interpretation aid, report $E[U^3]/q$, the probability that a third specified borrower defaults given that two others defaulted. This is the first order of joint default left free by the calibration. The entire distribution of $U$ matters for ES.

## Numerical protocol

1. **Finite-portfolio reference.** Given $U=u$, let $B\sim\mathrm{Bernoulli}(u)$ denote the designated borrower's default and $K\sim\mathrm{Binomial}(999,u)$ count defaults among the others. They are conditionally independent. For designated weight $s$, tabulate the conditional distribution of $g[sB+(1-s)K/999]$ and integrate over the common factors. Combine coincident loss outcomes using integer exposure units before computing quantiles. As a cross-check, the granular count conditional on $U$ is exactly binomial with size 1,000.
2. **Risk measures.** Set $v=\inf\{x:F_L(x)\geq\alpha\}$. Use $\mathrm{ES}_\alpha=v+E[(L-v)_+]/(1-\alpha)$. This handles probability mass at the cutoff. Evaluate the upper tail directly when possible. [Acerbi and Tasche (2002)](https://arxiv.org/abs/cond-mat/0104295)
3. **Independent simulation.** In a separate path, draw shared factors, generate latent borrower values and threshold defaults directly. Reuse Gaussian draws between models and the same default indicators across allocations. Process scenarios in chunks. Use fixed, recorded seeds and independent paired batches. Calculate the headline estimate from pooled scenarios. Use paired batch differences for Monte Carlo uncertainty and record realized precision and scenario count.
4. **Acceptance checks.** Recover $p$ and $q$ by quadrature to within $10^{-9}$ absolute. Verify loss probability sums to one, support lies in $[0,g]$, and moments recover the analytical targets. Check the independent Gaussian binomial benchmark at $r_G=0$ and a one-borrower portfolio where ES can be calculated directly. Refine quadrature until successive ES estimates move by less than 0.1 bp and the reported VaR cutoffs are stable. Check $F(v^-)<\alpha\leq F(v)$ against estimated numerical error. Flag unresolved cutoffs.
5. **Precision and scope.** Start with 32 independent batches of 50,000 scenarios (about 1,600 draws in the pooled worst 0.1%). Budget up to 6.4 million scenarios with the same batch count in the precision study. Target a 95% Monte Carlo interval half-width of 2 bp for the primary ES gap. Report actual precision even if the target is missed. Investigate material discrepancies against quadrature before interpreting the gap. Report any changes to the design alongside the results.

An optional secondary stress doubles marginal PD to 2% while freezing both models' latent parameters and changing their default thresholds. Recompute and disclose each pair-default probability in that scenario. It answers a separate question about extrapolating the baseline calibration. Macro transmission, regulatory capital and empirical model ranking are outside this experiment.

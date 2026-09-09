# Proof Memo: Conditional Monte Carlo for Gaussian Quadratic First Passage

## Claim

This memo proves an exact conditional Monte Carlo identity for a centered Gaussian first-passage event. It also states when the codebase evaluates that identity: the angular measure, chi-square degrees of freedom, and directional maximum must match the stated model. Consequently, for independent Monte Carlo samples, the estimator is unbiased and has variance no greater than the direct indicator estimator. This operation is commonly called Rao-Blackwellization in Monte Carlo; more precisely, it is conditional Monte Carlo or analytic marginalization of the radial variable.

## Assumptions and definitions

Let $r\geq1$, let $J\subseteq\mathbb{R}$ be an interval with left endpoint $x_0$, and let $G:J\to\mathbb{R}^{r\times r}$ be continuous and symmetric. Fix $x\in J$ and set $I_x=[x_0,x]$. Then $I_x$ is nonempty and compact. Let

$$
\boldsymbol{w}\sim N(\boldsymbol{0},I_r).
$$

Define the running maximum and its threshold indicator by

$$
M_x(\boldsymbol{w})=
\max_{\xi\in I_x}\boldsymbol{w}^TG(\xi)\boldsymbol{w},
\qquad
Y_x=\mathbf{1}\{M_x(\boldsymbol{w})\geq1\}.
$$

Define $h(\xi,\boldsymbol{u})=\boldsymbol{u}^TG(\xi)\boldsymbol{u}$. Joint continuity on the compact set $I_x\times S^{r-1}$ implies that $m_x$ is continuous, hence measurable, by the maximum theorem. The same argument makes $M_x$ continuous in $\boldsymbol{w}$. The extreme-value theorem makes both maxima finite and attained. Define the first-passage event at $x$ by

$$
\{X_{\mathrm{fp}}\leq x\}=\{M_x(\boldsymbol{w})\geq1\}.
$$

Therefore, the target probability is

$$
F_{\mathrm{fp}}(x)=\Pr(X_{\mathrm{fp}}\leq x)=E[Y_x].
$$

This event-based definition avoids any separate infimum-attainment condition. Equivalently, define $X_{\mathrm{fp}}=\inf\{\xi\in J:\boldsymbol{w}^TG(\xi)\boldsymbol{w}\geq1\}$ with $\inf\varnothing=+\infty$. Since $J$ starts at $x_0$ and the exceedance set is closed by continuity, the same event identity holds. Here, first passage means threshold attainment or exceedance, not equality crossing alone.

## Polar decomposition

Write the standard Gaussian vector as

$$
\boldsymbol{w}=R\boldsymbol{u},
$$

where

$$
R^2\sim\chi_r^2,
\qquad
\boldsymbol{u}\sim\operatorname{Unif}(S^{r-1}),
\qquad
R\perp\boldsymbol{u}.
$$

For each direction, define

$$
m_x(\boldsymbol{u})=
\max_{\xi\in I_x}\boldsymbol{u}^TG(\xi)\boldsymbol{u}.
$$

Quadratic homogeneity gives, almost surely,

$$
\begin{aligned}
M_x(R\boldsymbol{u})
&=\max_{\xi\in I_x}(R\boldsymbol{u})^TG(\xi)(R\boldsymbol{u})\\
&=R^2\max_{\xi\in I_x}\boldsymbol{u}^TG(\xi)\boldsymbol{u}\\
&=R^2m_x(\boldsymbol{u}).
\end{aligned}
$$

Thus

$$
Y_x=\mathbf{1}\{R^2m_x(\boldsymbol{u})\geq1\}.
$$

## Conditional Monte Carlo identity

Let $\mathcal{U}=\sigma(\boldsymbol{u})$. Since $R$ and $\boldsymbol{u}$ are independent,

$$
\begin{aligned}
Z_x
&=E[Y_x\mid\mathcal{U}]\\
&=\Pr\left(R^2m_x(\boldsymbol{u})\geq1\mid\boldsymbol{u}\right)\\
&=\mathbf{1}_{\{m_x(\boldsymbol{u})>0\}}
\overline F_{\chi_r^2}\left(\frac{1}{m_x(\boldsymbol{u})}\right).
\end{aligned}
$$

The stated assumptions make $m_x(\boldsymbol{u})<\infty$. When $m_x(\boldsymbol{u})\leq0$, the event $R^2m_x(\boldsymbol{u})\geq1$ is impossible and the indicator assigns $Z_x=0$. For positive finite $m_x$, continuity of the chi-square law makes $\Pr(R^2\geq1/m_x)=\Pr(R^2>1/m_x)$. Hence the angular integrand is exactly $E[Y_x\mid\mathcal{U}]$.

By the tower property,

$$
E[Z_x]=E[E[Y_x\mid\mathcal{U}]]=E[Y_x]=F_{\mathrm{fp}}(x).
$$

Thus analytic integration of the radius does not change the first-passage probability.

## Variance identity

Conditional on $\mathcal{U}$, $Y_x$ is Bernoulli with success probability $Z_x$. Therefore,

$$
\operatorname{Var}(Y_x\mid\mathcal{U})=Z_x(1-Z_x).
$$

The law of total variance gives

$$
\operatorname{Var}(Y_x)
=E[Z_x(1-Z_x)]+\operatorname{Var}(Z_x),
$$

or equivalently,

$$
\boxed{
\operatorname{Var}(Y_x)-\operatorname{Var}(Z_x)
=E[Z_x(1-Z_x)]\geq0.
}
$$

The inequality is strict exactly when

$$
\Pr\{0<Z_x<1\}>0.
$$

Under the finiteness assumption, $0<m_x(\boldsymbol{u})<\infty$ implies $0<Z_x<1$ because the chi-square law has positive density on $(0,\infty)$. Hence strict reduction occurs exactly when

$$
\Pr\{m_x(\boldsymbol{u})>0\}>0.
$$

For independent uniform directions $\boldsymbol{u}_1,\ldots,\boldsymbol{u}_K$,

$$
\widehat F_{\mathrm{CMC}}(x)=\frac{1}{K}\sum_{k=1}^{K}Z_x(\boldsymbol{u}_k)
$$

is unbiased, and

$$
\operatorname{Var}(\widehat F_{\mathrm{CMC}}(x))
=\frac{\operatorname{Var}(Z_x)}{K}
\leq\frac{\operatorname{Var}(Y_x)}{K}.
$$

For comparison, independent direct Gaussian trajectories $\boldsymbol{w}_1,\ldots,\boldsymbol{w}_K$ give $\widehat F_{\mathrm{direct}}(x)=K^{-1}\sum_kY_x(\boldsymbol{w}_k)$, with variance $F_{\mathrm{fp}}(x)[1-F_{\mathrm{fp}}(x)]/K$. This is a per-sample variance result. It does not itself compare computational cost. Deterministic quadrature and RQMC use the same exact radial integration, but this independent-sample variance formula does not apply directly to them.

## Code-to-model correspondence

`generateAngularRule` samples random directions by normalizing independent standard Gaussian vectors. Its RQMC rule maps each scrambled Sobol' point coordinatewise to a vector with standard-normal marginal coordinates, then normalizes it. Each unmodified point therefore targets the uniform spherical marginal law, but RQMC directions are dependent across rule points. Probability clipping before the inverse-normal map is a controlled numerical approximation. The circle parameterization applies only for $r=2$; its half-circle reduction is valid because $m_x(-\boldsymbol{u})=m_x(\boldsymbol{u})$. User-supplied rules require separate verification of their spherical measure and weights.

`computeTransitionCDF` evaluates the required tail as `gammainc(1 ./ (2 .* m), r / 2, 'upper')` for positive directional maxima and assigns zero transition probability otherwise. Thus, subject to floating-point evaluation, it computes $Z_x$ for the directional maxima supplied to it.

`computeRunningMaxGrid` computes the maximum over saved stations. At the mathematical-algorithm level, it therefore applies the theorem to the discrete first-passage event defined on that station set, not necessarily to a continuum event. `computeRunningMaxHermite` targets the maximum of a cubic-Hermite interpolant by evaluating endpoints and candidate stationary points on each interval. Cubic interpolation need not preserve positive semidefiniteness, but the theorem requires only symmetry. Subject to successful global critical-point evaluation and floating-point error, it therefore evaluates the conditional Monte Carlo integrand for the interpolated model. Its metadata records that this interpolation is not a certified continuum bound. Convergence or interpolation analysis is required before either numerical maximum can be identified with the exact continuum maximum.

## Scope

The proof applies to centered Gaussian vectors and homogeneous quadratic forms. Let $\boldsymbol{w}\sim N(\boldsymbol{0},\Sigma)$ in $\mathbb{R}^d$, let $q=\operatorname{rank}(\Sigma)\geq1$, and choose $L\in\mathbb{R}^{d\times q}$ with full column rank such that $\Sigma=LL^T$. Writing $\boldsymbol{w}=L\boldsymbol{z}$, where $\boldsymbol{z}\sim N(\boldsymbol{0},I_q)$, replaces $G(\xi)$ by $\widetilde G(\xi)=L^TG(\xi)L$ and the radial law by $\chi_q^2$. If $q=0$, then $\boldsymbol{w}=\boldsymbol{0}$ almost surely and the positive-threshold first-passage probability is zero.

Noncentered Gaussian inputs are a possible extension, but they do not have this simple radial factorization.

The conditional-expectation variance result is a specialization of Blackwell, "Conditional Expectation and Unbiased Sequential Estimation," *The Annals of Mathematical Statistics*, 18(1), 105-110 (1947), DOI: [10.1214/aoms/1177730497](https://doi.org/10.1214/aoms/1177730497). No sufficiency claim is required because conditioning on any sub-sigma-algebra preserves expectation and weakly reduces variance.

The real-valued implementation is checked against direct trajectory Monte Carlo in [tests/runSyntheticRegressionTests.m](tests/runSyntheticRegressionTests.m).

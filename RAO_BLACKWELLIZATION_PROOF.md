# Proof Memo: Conditional Monte Carlo for Gaussian Quadratic First Passage

## Claim

This memo proves an exact conditional Monte Carlo identity for a centered Gaussian first-passage event. It also states when the codebase evaluates that identity: the angular measure, chi-square degree of freedom, and running maximum all match the implemented quadratic-form tail probability.

## Assumptions and definitions

Let $r\geq1$, let $J\subseteq\mathbb{R}$ be an interval with left endpoint $x_0$, and let $G:J\to\mathbb{R}^{r\times r}$ be continuous and symmetric. Fix $x\in J$ and set $I_x=[x_0,x]$. Then $I_x$ is compact.

Let

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

Define $h(\xi,\boldsymbol{u})=\boldsymbol{u}^TG(\xi)\boldsymbol{u}$. Joint continuity on the compact set $I_x\times S^{r-1}$ implies that $m_x$ is continuous, hence measurable, by the maximum theorem.

The first-passage event is

$$
\{X_{\mathrm{fp}}\leq x\}=\{M_x(\boldsymbol{w})\geq1\}.
$$

Therefore, the target probability is

$$
F_{\mathrm{fp}}(x)=\Pr(X_{\mathrm{fp}}\leq x)=E[Y_x].
$$

This event-based definition avoids any separate infimum-attainment condition. Equivalently, define $X_{\mathrm{fp}}=\inf\{\xi\in J:\boldsymbol{w}^TG(\xi)\boldsymbol{w}\geq1\}$ with $\inf\varnothing=+\infty$.

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

The stated assumptions make $m_x(\boldsymbol{u})<\infty$. When $m_x(\boldsymbol{u})\leq0$, the event $R^2m_x(\boldsymbol{u})\geq1$ is impossible and the indicator assigns $Z_x=0$. For positive finite $m_x(\boldsymbol{u})$, the chi-square survival function is evaluated at the reciprocal threshold.

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

For comparison, independent direct Gaussian trajectories $\boldsymbol{w}_1,\ldots,\boldsymbol{w}_K$ give $\widehat F_{\mathrm{direct}}(x)=K^{-1}\sum_kY_x(\boldsymbol{w}_k)$, with variance $F_{\mathrm{fp}}(x)(1-F_{\mathrm{fp}}(x))/K$.

## Code-to-model correspondence

`generateAngularRule` samples random directions by normalizing independent standard Gaussian vectors. Its RQMC rule maps each scrambled Sobol' point coordinatewise to a vector with standard-normal marginals and then normalizes.

`computeTransitionCDF` evaluates the required tail as `gammainc(1 ./ (2 .* m), r / 2, 'upper')` for positive directional maxima and assigns zero transition probability otherwise. Thus, subject to floating-point error, it computes the same conditional expectation as above.

`computeRunningMaxGrid` computes the maximum over saved stations. At the mathematical-algorithm level, it therefore applies the theorem to the discrete first-passage event defined on that station set, not to the continuous-time event.

## Scope

The proof applies to centered Gaussian vectors and homogeneous quadratic forms. Let $\boldsymbol{w}\sim N(\boldsymbol{0},\Sigma)$ in $\mathbb{R}^d$, let $q=\operatorname{rank}(\Sigma)\geq1$, and choose any factorization $\Sigma=LL^T$. Then $L\boldsymbol{z}$ with $\boldsymbol{z}\sim N(\boldsymbol{0},I_q)$ reduces the problem to the identity above on the $q$-dimensional support.

Noncentered Gaussian inputs are a possible extension, but they do not have this simple radial factorization.

The conditional-expectation variance result is a specialization of Blackwell, "Conditional Expectation and Unbiased Sequential Estimation," *The Annals of Mathematical Statistics*, 18(1), 105-110 (1947).

The real-valued implementation is checked against direct trajectory Monte Carlo in [tests/runSyntheticRegressionTests.m](tests/runSyntheticRegressionTests.m).

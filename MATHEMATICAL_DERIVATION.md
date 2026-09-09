# Mathematical Derivation: First-Passage Probability for a Linear Spatial System

## 1. Linear spatial model

Let $\boldsymbol{q}(x) \in \mathbb{C}^{n}$ be a state that evolves from inlet location $x_0$ under the linear spatial operator

$$
\boldsymbol{q}(x)=\mathcal{P}(x,x_0)\boldsymbol{q}_0.
$$

Represent the random inlet state with $r$ independent real Gaussian coordinates:

$$
\boldsymbol{q}_0=B_0\boldsymbol{w},
\qquad
\boldsymbol{w}\sim N(\boldsymbol{0},I_r).
$$

Define the propagated factor

$$
B(x)=\mathcal{P}(x,x_0)B_0,
$$

so that

$$
\boldsymbol{q}(x)=B(x)\boldsymbol{w}.
$$

Let $H(x)\in\mathbb{C}^{n\times n}$ be Hermitian positive semidefinite. Define the scalar state measure

$$
e(x)=\boldsymbol{q}(x)^*H(x)\boldsymbol{q}(x).
$$

Substitution gives

$$
\begin{aligned}
e(x)
&=\boldsymbol{w}^T B(x)^*H(x)B(x)\boldsymbol{w}\\
&=\boldsymbol{w}^T A(x)\boldsymbol{w},
\end{aligned}
$$

where

$$
A(x)=\operatorname{Re}\!\left(B(x)^*H(x)B(x)\right)
\in\mathbb{R}^{r\times r}.
$$

For real $\boldsymbol{w}$, the imaginary skew-symmetric part of $B^*HB$ has zero quadratic form. Thus $A(x)$ is real symmetric positive semidefinite.

## 2. Normalized threshold event

Let $e_{\mathrm{thres}}(x)>0$ be a prescribed threshold. Define

$$
G(x)=\frac{A(x)}{e_{\mathrm{thres}}(x)}.
$$

Then $G(x)$ is real symmetric positive semidefinite, and the normalized exceedance event is

$$
e(x)\geq e_{\mathrm{thres}}(x)
\quad\Longleftrightarrow\quad
\boldsymbol{w}^T G(x)\boldsymbol{w}\geq 1.
$$

On an ordered spatial domain $J$, define the first-passage location by

$$
X_{\mathrm{fp}}
=\inf\left\{x\in J:
\boldsymbol{w}^TG(x)\boldsymbol{w}\geq1\right\},
$$

with $X_{\mathrm{fp}}=\infty$ when no threshold attainment or exceedance occurs. If $G$ is continuous on $[x_0,x]\subseteq J$, the exceedance set is closed.

For any $x$, the first-passage event is

$$
\{X_{\mathrm{fp}}\leq x\}
=\left\{\sup_{\xi\in[x_0,x]}
\boldsymbol{w}^TG(\xi)\boldsymbol{w}\geq1\right\}.
$$

Therefore, the first-passage CDF is

$$
F_{\mathrm{fp}}(x)=\Pr(X_{\mathrm{fp}}\leq x).
$$

## 3. Radial-angular decomposition

For $\boldsymbol{w}\sim N(\boldsymbol{0},I_r)$, write

$$
\boldsymbol{w}=R\boldsymbol{u},
$$

where

$$
R=\|\boldsymbol{w}\|_2,
\qquad
\boldsymbol{u}=\frac{\boldsymbol{w}}{\|\boldsymbol{w}\|_2}.
$$

The standard normal invariance gives

$$
R^2\sim\chi_r^2,
\qquad
\boldsymbol{u}\sim\operatorname{Unif}(S^{r-1}),
\qquad
R\perp\boldsymbol{u}.
$$

For a fixed direction $\boldsymbol{u}$, define the directional gain

$$
a_{\boldsymbol{u}}(x)=\boldsymbol{u}^TG(x)\boldsymbol{u}\geq0
$$

and its running maximum

$$
m_{\boldsymbol{u}}(x)=
\max_{\xi\in[x_0,x]}a_{\boldsymbol{u}}(\xi).
$$

Since

$$
\boldsymbol{w}^TG(\xi)\boldsymbol{w}
=R^2a_{\boldsymbol{u}}(\xi),
$$

the conditional first-passage event is

$$
\{X_{\mathrm{fp}}\leq x\mid\boldsymbol{u}\}
=\{R^2m_{\boldsymbol{u}}(x)\geq1\}.
$$

When $m_{\boldsymbol{u}}(x)>0$,

$$
\Pr(X_{\mathrm{fp}}\leq x\mid\boldsymbol{u})
=\Pr\left(R^2\geq\frac{1}{m_{\boldsymbol{u}}(x)}\right)
=\overline F_{\chi_r^2}\left(\frac{1}{m_{\boldsymbol{u}}(x)}\right).
$$

This conditional probability is zero when $m_{\boldsymbol{u}}(x)=0$. Taking expectation over the sphere yields

$$
\boxed{
F_{\mathrm{fp}}(x)=
E_{\boldsymbol{u}}\left[
\overline F_{\chi_r^2}\left(\frac{1}{m_{\boldsymbol{u}}(x)}\right)
\right].
}
$$

In the implementation, the chi-square survival function is evaluated as

$$
\overline F_{\chi_r^2}(z)=
\operatorname{gammainc}\left(\frac{z}{2},\frac{r}{2},\mathrm{upper}\right).
$$

## 4. Per-mode to composite construction

Let $M$ independent modes have centered Gaussian coordinate vectors

$$
\boldsymbol{w}_i\sim N(\boldsymbol{0},I_{r_i}),
\qquad i=1,\ldots,M,
$$

and mode energy matrices $A_i(x)$. The total energy is

$$
e_{\mathrm{total}}(x)=
\sum_{i=1}^{M}\boldsymbol{w}_i^TA_i(x)\boldsymbol{w}_i.
$$

Let $e_{\mathrm{thres}}(x)>0$ be one common total-energy threshold and define the mode-normalized matrices

$$
G_i(x)=\frac{A_i(x)}{e_{\mathrm{thres}}(x)}.
$$

Each mode may use an orthonormal reduction basis $V_i\in\mathbb{R}^{r_i\times k_i}$. With reduced coordinates $\boldsymbol{z}_i\sim N(\boldsymbol{0},I_{k_i})$, the reduced matrix is

$$
G_{i,\mathrm{red}}(x)=V_i^TG_i(x)V_i.
$$

Concatenate the independent reduced coordinates and form the block-diagonal composite matrix:

$$
\boldsymbol{z}=
\begin{bmatrix}
\boldsymbol{z}_1\\
\vdots\\
\boldsymbol{z}_M
\end{bmatrix}
\sim N(\boldsymbol{0},I_{r_{\mathrm{total}}}),
\qquad
r_{\mathrm{total}}=\sum_{i=1}^{M}k_i,
$$

$$
G_{\mathrm{total}}(x)=
\operatorname{blkdiag}\left(
G_{1,\mathrm{red}}(x),\ldots,G_{M,\mathrm{red}}(x)
\right).
$$

Then

$$
\boldsymbol{z}^TG_{\mathrm{total}}(x)\boldsymbol{z}
=\sum_{i=1}^{M}
\boldsymbol{z}_i^TG_{i,\mathrm{red}}(x)\boldsymbol{z}_i.
$$

Thus the composite system has the same normalized total-energy threshold event as the sum of reduced mode energies. The radial-angular derivation applies directly in dimension $r_{\mathrm{total}}$. Reduction changes the model through projection; block-diagonal assembly itself preserves the reduced independent-mode quadratic form exactly.

## 5. Related probability quantities

The survival probability is

$$
S_{\mathrm{fp}}(x)=1-F_{\mathrm{fp}}(x).
$$

For ordered reporting locations $x_1<\cdots<x_N$, the discrete interval mass is

$$
p_1=F_{\mathrm{fp}}(x_1),
\qquad
p_n=F_{\mathrm{fp}}(x_n)-F_{\mathrm{fp}}(x_{n-1}),\quad n\geq2.
$$

The right-censoring probability at $x_N$ is

$$
p_{\mathrm{cens}}=S_{\mathrm{fp}}(x_N),
$$

and probability conservation is

$$
\sum_{n=1}^{N}p_n+p_{\mathrm{cens}}=1.
$$

The local exceedance probability is

$$
p_{\mathrm{local}}(x)=
E_{\boldsymbol{u}}\left[
\overline F_{\chi_r^2}\left(\frac{1}{a_{\boldsymbol{u}}(x)}\right)
\right].
$$

Since $a_{\boldsymbol{u}}(x)\leq m_{\boldsymbol{u}}(x)$ and the survival function is decreasing,

$$
p_{\mathrm{local}}(x)\leq F_{\mathrm{fp}}(x).
$$

Equality holds at the inlet because $m_{\boldsymbol{u}}(x_0)=a_{\boldsymbol{u}}(x_0)$.

## 6. Angular quadrature estimator

Let $\{(\boldsymbol{u}_k,\omega_k)\}_{k=1}^K$ be directions and nonnegative normalized weights, with

$$
\sum_{k=1}^{K}\omega_k=1.
$$

The numerical estimator is

$$
\widehat F_{\mathrm{fp}}(x)=
\sum_{k=1}^{K}\omega_k
\overline F_{\chi_r^2}\left(\frac{1}{m_k(x)}\right),
$$

where

$$
m_k(x)=\max_{\xi\in[x_0,x]}\boldsymbol{u}_k^TG(\xi)\boldsymbol{u}_k.
$$

The same formula with $m_k(x)$ replaced by $a_k(x)$ estimates local exceedance. Using the same quadrature directions for both estimates preserves $\widehat p_{\mathrm{local}}(x)\leq\widehat F_{\mathrm{fp}}(x)$ up to floating-point error.

For a saved spatial grid, the maximum is over the grid-defined event. A cubic-Hermite maximum instead defines an interpolated event. Either numerical maximum agrees with the continuum model only after a separate grid-convergence or interpolation-error analysis.

For RQMC, let $\widehat F_j(x)$ be the estimate from scramble $j$, $j=1,\ldots,J$. Then

$$
\widehat F(x)=\frac{1}{J}\sum_{j=1}^{J}\widehat F_j(x),
$$

with replicate standard-error estimate

$$
\widehat{\operatorname{SE}}[\widehat F(x)]
=\frac{\operatorname{sd}\{\widehat F_j(x)\}_{j=1}^{J}}{\sqrt{J}}.
$$

## 7. Proper-complex variant

For a proper complex standard Gaussian vector $\boldsymbol{z}\sim\mathcal{CN}(\boldsymbol{0},I_r)$, use

$$
\boldsymbol{z}=R\boldsymbol{u},
\qquad
R^2\sim\operatorname{Gamma}(r,1).
$$

The directional maximum remains unchanged, but the conditional radial tail becomes

$$
\Pr(X_{\mathrm{fp}}\leq x\mid\boldsymbol{u})
=\operatorname{gammainc}\left(\frac{1}{m_{\boldsymbol{u}}(x)},r,\mathrm{upper}\right).
$$

Thus the derivation is identical after replacement of the real chi-square radial law by the proper-complex gamma radial law.

## 8. Validation references

The real-valued derivation is numerically checked by the analytic synthetic matrix reconstruction, angular refinement, direct trajectory Monte Carlo comparison, local-versus-first-passage inequality, and probability-conservation checks in [tests/runSyntheticRegressionTests.m](tests/runSyntheticRegressionTests.m).

The proper-complex radial law and its first-passage implementation are checked against a scalar analytic result, quadratic-form moments, and direct proper-complex trajectory Monte Carlo in [tests/runProperComplexRegressionTests.m](tests/runProperComplexRegressionTests.m).
# Conceptual Guide: First-Transition Probability

## The question

Suppose a linear spatial system carries a small uncertain input downstream. At each location, the state has an energy. We want the probability that the energy has crossed a threshold by that location.

This is a first-passage question. It is not the same as asking whether energy is above the threshold at one location. A realization can cross early, later decrease, and still count as transitioned.

## A two-coordinate picture

Consider two independent random inlet amplitudes,

$$
\boldsymbol{w}=
\begin{bmatrix}
w_1\\
w_2
\end{bmatrix},
\qquad
w_1,w_2\sim N(0,1).
$$

Each point $(w_1,w_2)$ is one possible inlet disturbance. The spatial operator propagates that disturbance. Because propagation is linear, the state at location $x$ is a linear combination of the same two inlet amplitudes:

$$
\boldsymbol{q}(x)=B(x)\boldsymbol{w}.
$$

The matrix $B(x)$ contains the spatial amplification and rotation of the two inlet patterns.

## Energy becomes an ellipse

Energy is a quadratic measure of the state. In the inlet-coordinate plane, it has the form

$$
e(x)=\boldsymbol{w}^TA(x)\boldsymbol{w}.
$$

At a fixed location, a constant-energy contour is usually an ellipse. The ellipse can change size and rotate as $x$ changes. A threshold defines a boundary in the $(w_1,w_2)$ plane. Points outside that boundary have crossed at that location.

After division by the threshold, the crossing condition is

$$
\boldsymbol{w}^TG(x)\boldsymbol{w}\geq1.
$$

## Why direct simulation is expensive

A direct Monte Carlo calculation draws many vectors $\boldsymbol{w}$. For each draw, it propagates the state or evaluates its energy at every saved location. It then records whether the energy ever exceeds the threshold.

This method uses a binary answer for each draw: crossed or did not cross. Binary indicators have high statistical noise, especially for rare events. Many draws are needed for a stable probability estimate.

## Separate size from direction

Every Gaussian vector can be written as

$$
\boldsymbol{w}=R\boldsymbol{u},
$$

where $R$ is its length and $\boldsymbol{u}$ is a unit direction. In two dimensions, $\boldsymbol{u}$ is an angle on a circle. The Gaussian distribution separates into:

- a direction that is uniform around the circle;
- a nonnegative radius with $R^2\sim\chi_2^2$.

The chi-square distribution is the distribution of the sum of squared independent standard normal variables. In $r$ coordinates, $R^2=w_1^2+\cdots+w_r^2\sim\chi_r^2$.

For a fixed direction, the normalized energy is

$$
\boldsymbol{w}^TG(x)\boldsymbol{w}
=R^2\boldsymbol{u}^TG(x)\boldsymbol{u}.
$$

The direction selects a ray from the origin. The radius selects how far out the realization lies along that ray.

## First passage along one direction

For one direction, calculate the directional gain at each location:

$$
a_{\boldsymbol{u}}(x)=\boldsymbol{u}^TG(x)\boldsymbol{u}.
$$

Then retain its running maximum:

$$
m_{\boldsymbol{u}}(x)=\max_{\xi\leq x}a_{\boldsymbol{u}}(\xi).
$$

This is the memory mechanism. If the directional gain was large at an earlier location, that earlier value remains active downstream.

For the selected direction, crossing by $x$ occurs when

$$
R^2m_{\boldsymbol{u}}(x)\geq1.
$$

Instead of randomly drawing $R$, the method calculates this probability exactly from the chi-square tail:

$$
\Pr\{\text{crossed by }x\mid\boldsymbol{u}\}
=\overline F_{\chi_r^2}\left(\frac{1}{m_{\boldsymbol{u}}(x)}\right).
$$

Here $\overline F$ is the survival function: the probability that a random variable is greater than a stated value.

## Average over directions

The final first-transition probability is the average of the conditional probabilities over all directions:

$$
F_{\mathrm{fp}}(x)=
E_{\boldsymbol{u}}\left[
\overline F_{\chi_r^2}\left(\frac{1}{m_{\boldsymbol{u}}(x)}\right)
\right].
$$

In two dimensions, this is an average around a circle. In more dimensions, directions lie on a sphere. The code uses circle quadrature, random directions, or scrambled Sobol' directions for this average.

Scrambled Sobol' sampling is a randomized quasi-Monte Carlo method. It places points more evenly than ordinary random sampling while retaining repeated randomized estimates for uncertainty assessment.

## Why the method is more efficient

The method removes one full random dimension: the Gaussian radius. Direct simulation repeatedly draws both direction and radius, then converts them into a binary crossing result. This method samples or integrates only directions and replaces the random radial decision with an exact probability.

This is conditional Monte Carlo, also called Rao-Blackwellization. It does not change the target probability. For independent random directions, it cannot increase estimator variance.

## Multiple independent modes

For several independent modes, each mode has its own Gaussian coordinates and energy contribution. Total energy is their sum.

The code first reduces each mode to its important stochastic directions. It then concatenates the reduced coordinates and creates one block-diagonal composite matrix. Block diagonal means that each mode has its own matrix block and there are no cross terms between independent modes.

The composite problem has the same form as the two-coordinate example, but with more coordinates. The radial-angular method then operates on the combined coordinate vector.

## Practical accelerations

- **Dimension reduction:** retain a small basis that captures transition-relevant energy growth.
- **Station downsampling:** construct matrices only at a selected subset of spatial locations.
- **Reduced cache:** store the smaller coordinate, matrix, and metadata representation instead of re-reading large raw solution files.
- **Parallel cache build:** build independent missing cache files concurrently within memory and CPU limits.
- **Hermite maxima:** use endpoint values and derivatives to resolve a cubic interpolant inside each saved spatial interval.

The grid maximum defines first passage on saved stations. The Hermite maximum defines first passage for the interpolated curve. Both require convergence checks before they represent the exact continuum model.

## A concise explanation

The code maps uncertain inlet disturbances into a low-dimensional Gaussian coordinate space. At each spatial location, energy becomes an ellipse in that space. First transition asks whether a point has ever left the moving threshold ellipse. The method examines one direction at a time, keeps the largest amplification seen so far, and calculates the probability of having sufficient Gaussian radius to cross. Averaging those exact conditional probabilities over directions gives the first-transition curve with less noise than direct trajectory sampling.
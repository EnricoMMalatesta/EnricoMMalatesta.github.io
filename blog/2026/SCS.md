[**Could be worse**](https://www.youtube.com/watch?v=GPfKEhFh9ww) - Random topics in statistical physics, spin glasses, optimization and average case hardness

# Chaos in Random Neural Networks


@def title = "Chaos in Random Neural Networks"
@def authors = "E. Malatesta"
@def published = "9 August 2026"
@def pt_lang = false
@def rss_pubdate = Date(2026, 8, 9)
@def rss = "random_neural_networks_scs_dmft"
@def rss_description = "Sompolinsky-Crisanti-Sommers random neural networks, the large-N dynamical mean-field reduction, and the Newtonian equation for the autocorrelation"


{{ published }} | **{{ authors }}**


@@post-top-nav
[← All blog posts](/pages/blog/)
@@


I recently went back to the classic paper of Sompolinsky, Crisanti and Sommers on chaos in random neural networks [\cite{SCS1988}], which was revisited in detail in [\cite{CS2018}].

The model they considered is a recurrent neural network with random asymmetric couplings. Its remarkable feature is that, using Dynamical Mean Field Theory (DMFT), one can show that in the large-$N$ limit the collective dynamics reduces to an effective one-neuron stochastic problem. The noise in this effective problem is generated self-consistently by the activity of the rest of the network.

In this first post I want to derive this reduction and the resulting closed equation for the autocorrelation function. I will then show how, in this model, the effective one-neuron problem can be recast as Newton's equation for a fictitious particle moving in a self-consistent potential. This mechanical analogy gives a very transparent way of organizing the formal DMFT solutions: static solutions, periodic solutions, and the special decaying solution associated with the chaotic phase.


<!-- 
I recently went back to the classic paper of Sompolinsky, Crisanti and Sommers on chaos in random neural networks [\cite{SCS1988}] which has been recently reviewed back in [\cite{CS2018}]. 

The model they considered is a recurrent neural network with random asymmetric couplings. The surprising result is that, using Dynamical Mean Field Theory (DMFT) it is possible to show that when the number of neurons is large, the collective dynamics can be reduced to a one-neuron effective stochastic problem whose noise is generated self-consistently by the rest of the network. 

In this first post I want to derive this reduction and the resulting equation for the autocorrelation function. I will also show how in this model this effective one neuron problem can be explicitly recasted as a Newton's equation, and the resulting physics (TODO: mention that the DMFT solutions corresponds to limit cycles, fixed points and chaos) can be unveiled easily. In the next post TODO: 1) which of the DMFT solutions corresponds to a  stable attractor in the dynamics (the chaotic one) 2) Lyapunov exponents of the chaotic solution -->


@@toc-title
Contents
@@
\toc


## The model

Consider $N$ continuous variables $h_i(t)$ $i=1,\dots,N$, representing the *post-synaptic potentials* of $N$ neurons. The activity, or output, of neuron $i$ is

$$
S_i(t)=\phi(g h_i(t)).
$$
Here $g$ is the gain and $\phi$ is an odd sigmoidal-shaped function. I will mostly keep the notation general, but the canonical example I will consider is $\phi(x)=\tanh x$. The dynamics of the post-synaptic potentials are given by the following $N$ coupled differential equations
$$
\label{eq::dynamical_equations}
\frac{d h_i}{dt} = -h_i(t)+\frac{1}{\sqrt{N}}\sum_{j=1}^N J_{ij}S_j(t).
$$

The first term is the *leak*: without any input, $h_i(t)$ decays exponentially to zero. The second term is the recurrent input from the rest of the network, which depends on the coupling matrix $J_{ij}$. Here we will consider random, independent, asymmetric couplings, having mean 0 and variance 1:
$$
\mathbb E_J[J_{ij}]=0,
\qquad
\mathbb E_J[J_{ij}^2]=1,
\qquad
\forall i \ne j = 1, \dots, N
%J_{ij}\ \text{independent of}\ J_{ji}.
$$
We also set $J_{ii}=0$, which is irrelevant at large $N$. The scaling $1/\sqrt{N}$ in equation \eqref{eq::dynamical_equations} maintains each recurrent input of order one as $N\to\infty$.

The configuration where all the neurons are silent
$$
h_i=0 \qquad \forall i
$$
is always a fixed point. Linearizing around it gives
$$
\partial_t h_i
=
-h_i+\frac{g}{\sqrt{N}}\sum_j J_{ij}h_j.
$$
Since the eigenvalues of $J/\sqrt{N}$ fill the unit disk in the complex plane at large $N$, the rightmost real part of the spectrum of $gJ/\sqrt{N}$ is $g$. Hence the zero fixed point is linearly stable for $g<1$ and loses stability at $g=1$. This simple calculation identifies the transition point, but it does not tell us what happens for $g>1$. For that we need the nonlinear mean-field theory.

## Dynamical Mean Field Theory

Define the recurrent input to neuron $i$:
$$
\eta_i(t)
\equiv
\frac{1}{\sqrt{N}}\sum_{j=1}^N J_{ij} \, \phi(g h_j(t)).
$$

The equation of motion becomes

$$
\partial_t h_i(t)=-h_i(t)+\eta_i(t).
$$


The idea of dynamical mean-field theory (DMFT) is that, at large $N$, $\eta_i(t)$ becomes a Gaussian process. In the following we will perform a simple, naive computation, which gives the right equations.


At fixed trajectories $\{h_j(t)\}$, $\eta_i(t)$ is a sum of many independent random variables. Its mean and covariance are:
$$
\begin{split}
\mathbb E_J[\eta_i(t)]&=0 \,, \\
\mathbb E_J[\eta_i(t)\eta_i(t')]
%&=
%\sum_{j,k} \mathbb E_J[J_{ij}J_{ik}] S_j(t)S_k(t') \\
&=
\frac{1}{N}
\sum_{j=1}^N
\phi(g h_j(t)) \phi(g h_j(t'))
\equiv C(t, t').
\end{split}
$$

At large $N$, the empirical average appearing in the covariance self-averages and converges to a deterministic two-time autocorrelation function $C(t,t')$. The original $N$-dimensional network can therefore be reduced to a single effective neuron driven by a Gaussian input $\eta(t)$ with covariance

$$
\boxed{
\langle \eta(t)\eta(t')\rangle=C(t,t')\,.} 
$$

Here and in the following, $\langle\bullet\rangle$ denotes the average over the effective Gaussian process. The dynamics of the effective neuron is therefore

$$
\label{eq::post_syn_pot}
\boxed{
\partial_t h(t)=-h(t)+\eta(t)
}
$$

The covariance of the effective noise must be determined self-consistently from the activity of the effective neuron via

$$
\boxed{
C(t,t') =
\left\langle
\phi(g h(t))\phi(g h(t'))
\right\rangle.
}
$$

Thus, the original dynamics of the $N$ dimensional network has been replaced by a by a stochastic single-neuron problem whose noise statistics are generated self-consistently by the neuron itself. This is the central dynamical mean-field equation.

A small subtlety is hidden in this argument. The trajectories $h_j(t)$ are generated by the same coupling matrix $J$, so the factors $\phi(g h_j(t))$ cannot literally be regarded as fixed independently of the couplings $J_{ij}$. The reason the conclusion is nevertheless correct is the full asymmetry and mean-field scaling of the couplings and can be understood with a cavity argument, see the footnote[^1]. The same result can be derived more systematically using the Martin–Siggia–Rose–Janssen–De Dominicis (MSRJD) generating-functional formalism; see [\cite{CS2018}]. For a pedagogical introduction to this method in the context of mean-field spin-glass dynamics, see [\cite{CavagnaCastellani2005}].




## The autocorrelation equation

We now focus on a stationary state where two-time quantities like the activity autocorrelation $C(t, t')$ depend only on the time difference $\tau=t-t'$
$$
C(\tau) = C(t - t')
$$
Similarly the post-synaptic field autocorrelation
$$
\Delta(\tau) = 
\left\langle h(t')h(t'+\tau)\right\rangle,
\qquad
\Delta_0 \equiv \Delta(0).
$$

Next, let's take the effective equation $(1+\partial_t)h(t)=\eta(t)$, multiply the left and the right hand side by $\eta(t')$ and take the average with respect to $\eta$:
$$
(1+\partial_t)(1+\partial_{t'})
\left\langle h(t)h(t')\right\rangle
=
\left\langle \eta(t)\eta(t')\right\rangle.
$$
The right and left hand side contain by definition repsectively $C(\tau)$ and $\Delta(\tau)$. Using the fact that $\partial_t=\partial_\tau$ and $\partial_{t'}=-\partial_\tau$ one finds that the first order derivative cancel, obtaining
$$
\label{eq::ODE_Delta}
\boxed{
  \Delta(\tau)-\partial_\tau^2\Delta(\tau)=C(\tau)
%\partial_{\tau}^2 \Delta(\tau) = -C(\tau) + \Delta(\tau).
}
$$

This is yet not enough, as we need to express $C$ self-consistently in terms of $\Delta$. The key simplification is that $h(t)$ is Gaussian. Indeed since the effective equation \eqref{eq::post_syn_pot} is linear in $h$ the solution can be written as
$$
 h(t)=\int_{-\infty}^t ds\, e^{-(t-s)}\eta(s).
$$
This is a linear functional of a Gaussian process, hence it is itself a Gaussian process. Therefore the pair $ (h(t),h(t+\tau))$ on which $\Delta(\tau)$ depends is a two-dimensional Gaussian vector. Its covariance matrix is
$$
\Sigma(\tau)
=
\begin{pmatrix}
\Delta_0 & \Delta(\tau) \\
\Delta(\tau) & \Delta_0
\end{pmatrix}.
$$
A two-dimensional Gaussian is completely fixed by this covariance matrix. Therefore $C(\tau)$ is a function of $\Delta(\tau)$ and $\Delta_0$:
$$
\boxed{
C(\tau)=F(\Delta(\tau);\Delta_0).
}
$$
We can represent the two correlated Gaussian variables using three independent standard Gaussians $x,y,z$:
$$
\begin{split}
%h&=\sqrt{\Delta_0-\Delta}\,x+\sqrt{\Delta}\,z, \\
%h'&=\sqrt{\Delta_0-\Delta}\,y+\sqrt{\Delta}\,z.
h&=\sqrt{\Delta_0-|\Delta|}\,x+\sqrt{|\Delta|}\,z \,, \\
h'&=\sqrt{\Delta_0-|\Delta|}\,y+\text{sign}(\Delta)\sqrt{|\Delta|}\,z.
\end{split}
$$
Then
$$
\langle h^2\rangle=\langle h'^2\rangle=\Delta_0,
\qquad
\langle hh'\rangle=\Delta.
$$
So we find, using the fact that $\phi$ is an odd function
$$
\label{eq::F}
F(\Delta;\Delta_0) = \mathrm{sign}({\Delta})\int Dz \left[ \int Dx\, \phi\left( g\sqrt{\Delta_0-|\Delta|}\,x + g\sqrt{ |\Delta|}\,z \right) \right]^2.
$$
where we have introduced the notation $Dx\equiv \frac{dx}{\sqrt{2\pi}}e^{-x^2/2}$. The closed DMFT equation is therefore
$$
\label{eq::final_DMFT_equation}
\boxed{
\Delta(\tau)-\ddot\Delta(\tau)
=
F(\Delta(\tau);\Delta_0).
}
$$
The unknown is a single function $\Delta(\tau)$, and all dependence on the nonlinear transfer function is hidden in the scalar function $F$.



### Properties of $\Delta(\tau)$

The equation \eqref{eq::ODE_Delta} can be inverted explicitly, using the face that the Green function of the operator $1-\partial_\tau^2$ is
$$
G(\tau)=\frac{1}{2}e^{-|\tau|},
$$
as $(1-\partial_\tau^2)G(\tau)=\delta(\tau)$. Therefore we find
$$
\label{eq::Sol_ODE_Delta}
\boxed{
\Delta(\tau) = \frac{1}{2} \int_{-\infty}^{+\infty} d\tau'\, e^{-|\tau-\tau'|} C(\tau'). 
}
$$

This formula is useful because it gives the properties of $\Delta$ directly. First, $\Delta$ is an even function $\Delta(-\tau)=\Delta(\tau)$ (the same is of course true for $C$). Secondly, $\Delta$ is differentiable function with zero derivative in the origin[^2]
$$
\boxed{
\dot\Delta(0)=0.
}
$$
Third, since $\Delta$ is an autocorrelation, it obeys the Cauchy-Schwarz bound
$$
|\Delta(\tau)|
=
|\langle h(t)h(t+\tau)\rangle|
\le
\sqrt{\langle h(t)^2\rangle\langle h(t+\tau)^2\rangle}
=
\Delta_0.
$$
Thus acceptable solutions must satisfy
$$
\boxed{
|\Delta(\tau)|\le \Delta_0.
}
$$
Collecting everything, the steady-state DMFT solutions are curves $\Delta(\tau)$ satisfying
$$
\Delta(0)=\Delta_0,
\qquad
\dot\Delta(0)=0,
\qquad
|\Delta(\tau)|\le \Delta_0.
$$
The initial condition $\Delta_0$ must be chosen so that the resulting solution is self-consistent.


## The Newtonian particle analogy


Much insight can be gained into the problem by rewriting the DMFT equation \eqref{eq::final_DMFT_equation} as

$$
\label{eq::Newton_equation}
\frac{d^2 \Delta}{d \tau^2} = \Delta-F(\Delta;\Delta_0) = -\frac{\partial V}{\partial \Delta}
$$

where we have defined the potential[^3]

$$
\label{eq::potential}
\boxed{
V(\Delta;\Delta_0)
=
-\frac{\Delta^2}{2}
+
\int_0^\Delta du\,F(u;\Delta_0).
}
$$

The DMFT equation therefore turns out to be a Newton's equation for a particle with position $\Delta$, moving in the time variable $\tau$, in the potential $V$ starting at position $\Delta_0$ with zero velocity. The energy

$$
\boxed{
E
=
\frac{1}{2}\dot\Delta^2+V(\Delta;\Delta_0)
}
$$

is conserved. The important point is that the potential itself depends on $\Delta_0$. By changing $\Delta_0$ we are not only changing the initial point of the Newtonian particle; we are also changing the shape of the potential in which it moves. We will see how the potential changes very concretely in the next section. 

Therefore the problem is self-consistent in a very concrete mechanical sense. We must choose a starting point $\Delta_0$, build the potential $V(\cdot;\Delta_0)$, release the particle from rest at $\Delta=\Delta_0$, and then check whether the orbit is an admissible autocorrelation.


### The shape of the potential

In order to classify the possible solutions to the DMFT equation \eqref{eq::Newton_equation},  it is useful to understand how the potential can look like vs $\Delta$, for each chosen initial condition $\Delta_0$.

Start noticing that for an odd transfer function, such as $\phi(x)=\tanh x$, the function $F(\Delta;\Delta_0)$ is odd in $\Delta$. Therefore the potential $V(\Delta; \Delta_0)$ in \eqref{eq::potential} is an even function of $\Delta$. Hence it is enough to study its shape when $\Delta>0$.


In the next discussion we will argue that the potential cannot have an arbitrary shape. We can see this by computing the first three derivatives of the potential. Using the fact illustrated in the footnote[^3], it is easy to show that, for $\Delta>0$

$$
\frac{\partial V}{\partial \Delta}
=
-\Delta+F(\Delta;\Delta_0),
$$

$$
\frac{\partial^2 V}{\partial \Delta^2}
=
-1+
\frac{\partial F}{\partial \Delta} = g^2
\int Dz
\left[
\int Dx\,
\phi'\left(
g\sqrt{\Delta_0-\Delta}\,x
+
g\sqrt{\Delta}\,z
\right)
\right]^2 - 1
$$


$$
\label{eq::third_derivative_potential}
\frac{\partial^3 V}{\partial \Delta^3} = \frac{\partial^2 F}{\partial \Delta^2}
=
g^4
\int Dz
\left[
\int Dx\,
\phi''\left(
g\sqrt{\Delta_0-\Delta}\,x
+
g\sqrt{\Delta}\,z
\right)
\right]^2
>0.
$$

Therefore $\partial^2_\Delta V$ is a monotonically increasing function of $\Delta$ for $\Delta>0$. This is a rather strong constraint: the curvature can change sign at most once. As a consequence the potential can only have few qualitative shapes. We can distinguish between them by looking at the sign of the curvature of the potential at the origin. The first possibility is 
$$ 
\frac{\partial^2 V}{\partial \Delta^2}(0;\Delta_0)\ge0. 
$$ 

Since the curvature is increasing for $\Delta>0$, the potential is convex. Because $V$ is even, $\Delta=0$ is then the unique minimum. The potential is a *single well*. The second possibility is

$$
\frac{\partial^2 V}{\partial \Delta^2}(0;\Delta_0)<0.
$$

Then $\Delta=0$ is a local maximum. Since the curvature can change sign at most once, the potential can bend upward at only once in the physical interval $\Delta_0 > \Delta > \Delta$. If $\partial_{\Delta} V$ crosses zero before the endpoint $\Delta_0$, this crossing gives a unique minimum at positive $\Delta$, and by evenness another minimum at negative $\Delta$. In that case the potential has the shape of a *double well*. 

If $\partial_{\Delta} V$ does not cross zero before the endpoint $\Delta_0$, then we have a *downhill potential* as it is monotonically decreasing in the whole interval $\Delta_0 > \Delta > \Delta$. As shown in the footnote[^4] this case only happens when $g\le 1$, which implies that in this regime the DMFT equation admits only the solution $\Delta=0$. 

Thus the monotonicity of $\partial_\Delta^2 V$ tells us that only very simple shapes are possible: a single well centered at the origin, a double well, or a *downhill potential* on the admissible interval.


The boundary between the single-well and double-well shapes given by the condition

$$
\label{eq::single_to_double_well_transition}
\boxed{g \int Dx\, \phi'\left(g\sqrt{\Delta_0}\,x\right) = 1}
$$


In the case of the $\tanh$ function, for small $\Delta_0$, the integral is close to one, so if $g>1$ the curvature at the origin is positive and the potential is a single well. As $\Delta_0$ increases, the Gaussian variable $g\sqrt{\Delta_0}x$ explores the saturated part of $\tanh$, where $\phi' = \operatorname{sech}^2$ is small. The integral decreases, the curvature at the origin eventually becomes negative, and the potential turns into a double well. 


### Static solutions

The simplest possible solution of the DMFT equation corresponds to an autocorrelation function independent of time

$$
\Delta(\tau)=\Delta_0.
$$

This corresponds to the particle not moving and sitting at the stationary points of the potential. This gives the condition

$$
\label{eq::static_solutions_equation}
\boxed{
\Delta_0=F(\Delta_0;\Delta_0) = \int Dz\,\phi(g\sqrt{\Delta_0}\,z)^2.}
$$

For $\phi(x)=\tanh x$, $\Delta_0=0$ is always a solution. This is the silent fixed point of the original network. The argument of the previous subsection show in the footnote[^4] shows that for $g\le1$ this is the only admissible stationary DMFT solution. For $g>1$ equation \eqref{eq::static_solutions_equation} has also a nonzero solution $\Delta_0>0$. 


### Time-dependent solutions

We here analyze other solutions of the DMFT equation \eqref{eq::Newton_equation} for $g>1$ which are time-dependent. Since the motion conserves energy, at each time the energy stays equal to the initial one corresponding to the one of the particle being at rest at $\Delta_0$ 
$$
E=V(\Delta_0;\Delta_0).
$$

Different shapes of the potential lead to different types of formal DMFT solutions. If the initial condition $\Delta_0$ is such that $V(\Delta_0; \Delta_0) \ne 0$ then the autocorrelation diplays an oscillatory, periodic profile. This corresponds to a *limit cycle* in the original network. Moreover if $V(\Delta_0; \Delta_0) > 0$, then the $\Delta(\tau)$ changes sign (when the particle reaches $\Delta = 0$). If instead $V(\Delta_0; \Delta_0) < 0$, then the particle is confined into the positive well of the potential, and the oscillations do not change sign. 


There is also a special solution that separates those two regimes when

$$
\boxed{
V(\Delta_0;\Delta_0)=0.
}
$$

This solution decays to zero when for $\tau$ large

$$
\Delta(\tau)\to0
\qquad
\text{as}\qquad
\tau\to\infty.
$$

This suggests that the underlying neural dynamics is *chaotic*, as the networks tends to forget the initial condition. 


### Numerical illustration

We can clarify the previous classificationn of DMFT solutions, by simply solving \eqref{eq::Newton_equation} numerically for each given initial condition $\Delta_0$ and showing the corresponding potential $V(\Delta; \Delta_0)$.  


@@codefile
```julia

module SCS 

using OrdinaryDiffEq, QuadGK, Roots

# ------------------------------------------------------------
# Gaussian integration
# ------------------------------------------------------------

const ∞ = 20.0
const dx = 0.5
const interval = map(x->sign(x)*abs(x)^2, -1:dx:1) .* ∞

G(x) = exp(-x^2/2) / √(2π)

∫D(f, int=interval) = quadgk(z->begin
        r = G(z) .* f(z)
        isfinite(r) ? r : 0.0
    end, int..., atol=1e-5, maxevals=10^5)[1]

# ------------------------------------------------------------
# Transfer function, its derivative and primitive
# ------------------------------------------------------------

ϕ(x) = tanh(x)
∂ϕ(x) = 1 - tanh(x)^2
Φ(x) = abs(x) + log1p(exp(-2abs(x))) - log(2) # Stable version of log(cosh(x))

# ------------------------------------------------------------
# F(Δ; Δ0)
# ------------------------------------------------------------

function F(Δ, Δ0, g)
    abs(Δ) < 1e-12 && return 0.0

    if Δ == Δ0
         return ∫D(z -> ϕ(g * √Δ0 * z)^2)
    end

    A = min(abs(Δ), Δ0)
    a = √(max(Δ0 - A, 0.0))
    b = √A

    m(z) = ∫D(x -> ϕ(g * (a * x + b * z)))

    return sign(Δ) * ∫D(z -> m(z)^2)
end

# ------------------------------------------------------------
# Potential V(Δ; Δ0)
# ------------------------------------------------------------

function V(Δ, Δ0, g)

    m1 = ∫D(x -> Φ(g * √Δ0 * x))

    if Δ == Δ0
        m1 = ∫D(z -> Φ(g * √Δ0 * z))
        m2 = ∫D(z -> Φ(g * √Δ0 * z)^2)
        
        return -0.5 * Δ0^2 + (m2 - m1^2) / g^2
    end

    A = min(abs(Δ), Δ0)
    a = √max(Δ0 - A, 0.0)
    b = √A

    m(z) = ∫D(x -> Φ(g * (a * x + b * z)))

    return - 0.5 * Δ^2 + (∫D(z -> m(z)^2) - m1^2) / g^2
end

# Curvature of the potential at the origin
function Vcurv0(Δ0, g)
    m = ∫D(x -> ∂ϕ(g * √Δ0 * x))
    return g^2 * m^2 - 1
end

# ------------------------------------------------------------
# The three phase-diagram curves
# ------------------------------------------------------------

# Boundary between single-well and double-well potential
function Δ0_boundary(g; Δmin = 1e-9, Δmax = 1.2)
    g < 1 && return NaN

    return find_zero(Δ0 -> Vcurv0(Δ0, g), (Δmin, Δmax), Bisection())
end

# Decaying solution: selects the initial condition such that Δ(τ) → 0 as τ → ∞
function Δ0_decay(g; Δmin = 1e-9, Δmax = 1.2)
    g < 1 && return NaN

    return find_zero(Δ0 -> V(Δ0, Δ0, g), (Δmin, Δmax), Bisection())
end

# Static, fixed point solution
function Δ0_static(g; Δmin=1e-9, Δmax = 1.2)
    g < 1 && return NaN

    return find_zero(Δ0 -> F(Δ0, Δ0, g) - Δ0, (Δmin, Δmax), Bisection())
end

# ------------------------------------------------------------
# Numerical solution of Newton's equation:
#
# Δ'' = Δ - F(Δ; Δ0)
# ------------------------------------------------------------

function rhs!(du, u, p, t)
    Δ, v = u
    Δ0, g = p

    du[1] = v
    du[2] = Δ - F(Δ, Δ0, g)
end

function orbit(Δ0, g; T = 60.0, dt = 0.05)
    init_cond = [Δ0, 0.0]
    tspan = (0.0, T)
    params = (Δ0, g)
    prob = ODEProblem(rhs!, init_cond, tspan, params)
    return solve(prob, Tsit5(); saveat = dt, abstol = 1e-8, reltol = 1e-8)
end


end

```
@@

We recap here also the three special values of $\Delta_0$ given by the conditions:
$$
\label{eq::single_to_double}
V''(0;\Delta_0)=0,
$$

$$
\label{eq::chaos}
V(\Delta_0;\Delta_0)=0,
$$

$$
\label{eq::static}
V'(\Delta_0; \Delta_0) = 0
%F(\Delta_0;\Delta_0)=\Delta_0.
$$

The first equation defines boundary between the single-well and double-well shapes;
the second one selects the decaying solution $\Delta(\tau) \to 0$ for $\tau \to \infty$ and finally the third one gives the static solution. Those three conditions are implemented respectively in the functions `Δ0_boundary`, `Δ0_decay` and `Δ0_static` in the julia code above.



In [Figure 1](#fig-scs-potentials) we show an example of the shapes of the potential for $g=2$ for several values of $\Delta_0$. In [Figure 2](#fig-scs-orbits) I show the corresponding trajectories $\Delta(\tau)$. The script below was used to produce those two figures.



@@codefile
```julia
using Plots, Plots.PlotMeasures, LaTeXStrings, Printf
# ------------------------------------------------------------
# Example: g = 2
# ------------------------------------------------------------

g = 2

Δ_boundary = SCS.Δ0_boundary(g)
Δ_decay = SCS.Δ0_decay(g)
Δ_static = SCS.Δ0_static(g)

Δ0s = [Δ_boundary, (Δ_boundary+Δ_decay)/2, Δ_decay, (Δ_decay + Δ_static)/2, Δ_static]

# ------------------------------------------------------------
# Plot potentials
# ------------------------------------------------------------

pV = plot(xlabel = L"\Delta", ylabel = L"V(\Delta;\Delta_0)", legend = (0.45, 0.3), palette = :Set1)

for Δ0 in Δ0s
    Ds = range(-Δ0, Δ0; length = 400)
    lab = latexstring(@sprintf("\\;\\Delta_0=%.3f", Δ0))
    plot!(pV, Ds, [SCS.V(D, Δ0, g) for D in Ds], label = lab)
    scatter!(pV, [Δ0], [SCS.V(Δ0, Δ0, g)], label = false,  primary = false, markersize = 3)
end


hline!(pV, [0.0], color = :black, linestyle = :dash, label = false)

display(pV)


# ------------------------------------------------------------
# Plot Δ(τ)
# ------------------------------------------------------------

Δ0s = [Δ_boundary, (Δ_boundary+Δ_decay)/2, Δ_decay, (Δ_decay + Δ_static)/2, Δ_static]

pD = plot(xlabel = L"\tau", ylabel = L"\Delta(\tau)", legend = false, palette = :Set1)


for Δ0 in Δ0s
    sol = SCS.orbit(Δ0, g; T = 60.0, dt = 0.05)
    plot!(pD, sol.t, [u[1] for u in sol.u])
end

hline!(pD, [0.0], color = :black, linestyle = :dash)

display(pD)

```
@@


\label{fig-scs-potentials}
@@center ![Potentials for several values of Delta0](/assets/images/blog/scs_dmft_potentials.png) @@
@@center *Figure 1: Self-consistent potentials $V(\Delta;\Delta_0)$ for $g=2$ and several values of $\Delta_0$. The dots mark the release points $\Delta=\Delta_0$. For small $\Delta_0$ the potential is a single well centered at the origin. Increasing $\Delta_0$ changes the curvature at the origin and produces a double-well shape. The value $\Delta_0\simeq0.481$ corresponds to the zero-energy chaotic solution, while $\Delta_0\simeq0.530$ is the static solution.* @@


\label{fig-scs-orbits}
@@center ![Autocorrelation trajectories for several values of Delta0](/assets/images/blog/scs_dmft_orbits.png) @@
@@center *Figure 2: Solutions $\Delta(\tau)$ of the Newtonian DMFT equation \eqref{eq::Newton_equation} for $g=2$. I have used the same initial conditions as in [Figure 1](#fig-scs-potentials).* 
@@


### Phase Diagram

Finally we can compute the three boundaries given in \eqref{eq::single_to_double}, \eqref{eq::chaos} and \eqref{eq::static} for different values of $g$ and draw a phase diagram in the $(\Delta_0,1/g)$ plane. This is done in the julia script below, and [Figure 3](#fig-scs-potentials) summarizes the result. Notice that all three curves meet at the transition point $g=1$, $\Delta_0=0$. 

@@codefile
```julia
# ------------------------------------------------------------
# Plot phase diagram
# ------------------------------------------------------------

gs = range(1.001, 100.0; length = 500)

Δ_boundary_curve = Float64[]
Δ_decay_curve = Float64[]
Δ_static_curve = Float64[]

for gval in gs
    push!(Δ_boundary_curve, SCS.Δ0_boundary(gval; Δmin = 1e-8, Δmax = 1.2))
    push!(Δ_decay_curve,   SCS.Δ0_decay(gval;   Δmin = 1e-8, Δmax = 1.2))
    push!(Δ_static_curve,  SCS.Δ0_static(gval;  Δmin = 1e-8, Δmax = 1.2))
end

push!(Δ_boundary_curve, 2/π)
push!(Δ_decay_curve,   2*(1-2/π))
push!(Δ_static_curve,  1.0)
gs = push!(collect(gs), Inf)


plt = plot(xlabel = L"\Delta_0", ylabel = L"1/g", legend = :topright, xlims = (0.0, 1.05), ylims = (0.0, 1.1)  )

plot!(plt, Δ_boundary_curve, 1.0 ./ gs, label = L"\;V''(0;\Delta_0)=0")
plot!(plt, Δ_decay_curve, 1.0 ./ gs, label = L"\;V(\Delta_0;\Delta_0)=0")
plot!(plt, Δ_static_curve, 1.0 ./ gs, label = L"\;F(\Delta_0;\Delta_0)=\Delta_0")

display(plt)
```
@@

\label{#fig-scs-phase-diagram}
@@center ![Phase diagram of the one-replica DMFT solutions](/assets/images/blog/scs_dmft_phase_diagram.png) @@
@@center *Figure 3: Phase diagram in the $(\Delta_0,1/g)$ plane. The orange curve is defined by $V''(0;\Delta_0)=0$ and separates the single-well and double-well shapes of the potential. The light blue curve is defined by $V(\Delta_0;\Delta_0)=0$ and to the DMFT solution with $C(\tau) \to 0$ for $\tau \to \infty$. The green curve is defined by $F(\Delta_0;\Delta_0)=\Delta_0$ and corresponds to the static solution. Therefore below the light blue curve the autocorrelation is periodic and sign changing whereas between the light blue and the green curve it is periodic with positive sign. Above the green line there are no non-trivial solutions to the DMFT equation except for $\Delta(\tau) = 0$. See also [\cite{CS2018}].* @@


## What's next?

I will deliberately leave two questions for a second post. First, for $g>1$, the one-replica DMFT equation admits several formal solutions, and one has to understand which of them corresponds to a stable attractor of the original neural network dynamics. We will find out that the only physical DMFT solution is the decaying one. Secondly, one can ask whether this solution is genuinely chaotic. This leads naturally to the computation of the Lyapunov exponent.


<!--### Close to the transition

For completeness, let us also record what the decaying solution looks like just above the transition. Put

$$
\sigma=g-1,
\qquad
0<\sigma\ll1.
$$

Then the amplitude of the nontrivial solution is small. Expanding $\tanh(g h)$ and using Wick's theorem gives, to leading order, a quartic potential of the form

$$
V(\Delta;\Delta_0)
\simeq
-\frac{\sigma^2}{6}\Delta^2
+
\frac{1}{6}\Delta^4,
$$

where the self-consistency condition fixes

$$
\Delta_0=\sigma+O(\sigma^2).
$$

The zero-energy equation

$$
\frac{1}{2}\dot\Delta^2+V(\Delta;\Delta_0)=0
$$

then gives

$$
\boxed{
\Delta(\tau)
\simeq
\sigma\,
\operatorname{sech}\left(\frac{\sigma\tau}{\sqrt 3}\right),
\qquad
\sigma=g-1\to0^+.
}
$$

Thus the equal-time variance turns on continuously as

$$
\Delta_0\sim g-1,
$$

while the correlation time diverges as

$$
\tau_c\sim \frac{\sqrt 3}{g-1}.
$$

This is the first signature of the transition: as $g\downarrow1$ from above, the fluctuating solution becomes weaker and slower, until it merges with the zero solution at $g=1$.

What remains is the question we have deliberately postponed. For $g>1$, the one-replica DMFT equation admits several formal solutions: static ones, periodic ones, and the decaying one. Which one is the physical attractor of the original network? And how do we show that the decaying solution is genuinely chaotic rather than merely irregular? These require the two-replica stability analysis and the Lyapunov exponent, which will be the subject of the next post.
-->


<!-- FOOTNOTES -->

[^1]: Remove neuron $i$ and denote by $h_j^{(i)}(t)$ the trajectories of the remaining network. These cavity trajectories are independent of the couplings $J_{ij}$ entering the removed neuron. Consequently, $$ \eta_i^{\mathrm{cav}}(t)= \frac{1}{\sqrt N} \sum_j J_{ij}\phi(g h_j^{(i)}(t)) $$ is, at large $N$, a Gaussian process with covariance $N^{-1}\sum_j\phi(g h_j^{(i)}(t))\phi(g h_j^{(i)}(t'))$. Adding neuron $i$ back perturbs each of the other trajectories only by $O(N^{-1/2})$. Adding neuron $i$ back perturbs neuron $j$ through the connection $J_{ji}/\sqrt N$. To first order in this perturbation, $$ \delta S_j(t) = S_j(t)-S_j^{(i)}(t) \simeq \frac{1}{\sqrt N} \int ds \, R_{j}(t,s) \, J_{ji} \, S_i(s)\,,$$ where $R_j(t,s)$ measures the response of neuron $j$ at time $t$ to a perturbation at time $s$. Each individual perturbation is therefore of order $N^{-1/2}$. Substituting this into the input to neuron $i$, $\eta_i(t)=\eta_i^{\rm cav} (t)+\delta\eta_i(t)$, gives $$ \delta\eta_i(t) \simeq \frac{1}{N} \sum_{j} J_{ij} J_{ji} \int ds \, R_{j}(t,s) \, S_i(s)\,.$$ The resulting correction to the input of neuron $i$ involves products $J_{ij}J_{ji}$. For fully asymmetric couplings, $J_{ij}$ and $J_{ji}$ are independent and centered, so these feedback contributions add incoherently and vanish in the large $N$ limit. Hence the true input and the cavity input coincide to leading order as $N\to\infty$. If reciprocal couplings were correlated, this cancellation would not occur and an additional retarded self-interaction term would survive in the effective dynamics. 

[^2]: Differentiating \eqref{eq::Sol_ODE_Delta} gives $$ \dot\Delta(\tau) = \frac{1}{2} \left[ \int_{\tau}^{+\infty}d\tau'\,e^{-(\tau'-\tau)}C(\tau') - \int_{-\infty}^{\tau}d\tau'\,e^{-(\tau-\tau')}C(\tau') \right].$$ At $\tau=0$, using the evenness of $C$, the two integrals are equal.

[^3]: Note that $$ \partial_\Delta F(\Delta; \Delta_0) = g^2 \mathrm{sign}(\Delta) \int Dz \left[ \int Dx\, \phi'\left( g\sqrt{\Delta_0-|\Delta|}\,x + g\sqrt{|\Delta|}\,z \right) \right]^2, $$ i.e. the derivative with respect to $\Delta$ of $F$ is form the same with the difference that the derivative of the function $\phi$ is used. Therefore the potential can be written as $$ V(\Delta;\Delta_0) = -\frac{\Delta^2}{2} + \int_0^\Delta du\,F(u;\Delta_0).$$ Therefore denoting by $\Phi(x) = \int_{-\infty}^x dy \, \phi(x)$ the primitive of $\phi$ we can write the potential in terms of the integrated outputs $$ V(\Delta;\Delta_0) = -\frac{\Delta^2}{2} + \frac{1}{g^2} \left. \int Dz \left[ \int Dx\, \Phi\left( g\sqrt{\Delta_0-|u|}\,x + g\sqrt{|u|}\,z \right) \right]^2\right|_{u=0}^{u=\Delta}$$


[^4]: For a saturating non-linearity as $\phi(x)=\tanh x$ one finds $$ 0\le \frac{\partial F}{\partial \Delta} \le g^2 \le 1\,. $$ Moreover, since $F(0;\Delta_0)=0$, for $0<\Delta\le\Delta_0$, we have that $$ F(\Delta;\Delta_0)<\Delta\,. $$ Hence $$ \frac{\partial V}{\partial \Delta} = -\Delta+F(\Delta;\Delta_0) <0 \,. $$ So, below the transition, the potential is strictly decreasing as we move from $0$ to $\Delta_0$. A particle released from rest at $\Delta=\Delta_0$ feels a positive force which moves it to values larger than $\Delta_0$. This is forbidden for an autocorrelation, because any admissible solution must satisfy $|\Delta(\tau)|\le\Delta_0.$ Thus no nonzero solution is possible for $g\le1$. The only admissible stationary DMFT solution is $\Delta(\tau)= 0$. 




@@notoc
References
@@

[1] \biblabel{SCS1988}{Sompolinsky, Crisanti and Sommers (1988)}  H. Sompolinsky, A. Crisanti and H. J. Sommers, "Chaos in Random Neural Networks", *Physical Review Letters* **61**, 259--262 (1988).

[2] \biblabel{CS2018}{Crisanti and Sompolinsky (2018)} A. Crisanti and H. Sompolinsky, "Path Integral Approach to Random Neural Networks", arXiv:1809.06042 (2018).

[3] \biblabel{CavagnaCastellani2005}{Castellani and Cavagna (2005)} T. Castellani and A. Cavagna, "Spin-Glass Theory for Pedestrians", *Journal of Statistical Mechanics: Theory and Experiment* **2005**, P05012 (2005).


{{ blogcomments }}

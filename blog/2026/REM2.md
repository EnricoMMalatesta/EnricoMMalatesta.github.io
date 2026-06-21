[**Could be worse**](https://www.youtube.com/watch?v=GPfKEhFh9ww) - Random topics in statistical physics, spin glasses, optimization and average case hardness

# The Random Energy Model -- Part II


@def title = "The Random Energy Model -- Part II"
@def authors = "E. Malatesta"
@def published = "15 March 2026"
@def pt_lang = false
@def rss_pubdate = Date(2026, 3, 15)
@def rss = "REM2"
@def rss_description = """REM part 2"""

{{ published }} | **{{ authors }}**

@@post-top-nav
[← All blog posts](/pages/blog/)
@@

In a previous [post](/blog/2025/REM1) we rigorously solved the Random Energy Model (REM) unveiling the presence of a ''condensation'' or ''freezing'' transition at an inverse temperature 
$$
\beta_c \equiv \frac{e_0}{J} = \sqrt{\frac{2 \ln 2}{J}} \,.
$$
Namely, for $\beta<\beta_c$ the Gibbs measure is spread over exponentially many configurations as the corresponding entropy of the model is positive, while for $\beta > \beta_c$ it condenses onto a sparse set of states corresponding to the lowest available energy $-e_0 = - \sqrt{2J \ln 2}$. This is evident in the expression of the energy as a function of temperature
$$
e(\beta) = \frac{\partial (\beta f)}{\partial \beta} =
\begin{cases}
- \beta J & \; \mathrm{if} \;\beta < \beta_c \\
- \sqrt{2 J \ln 2} & \; \mathrm{if} \; \beta \ge \beta_c \,.
\end{cases}
$$
The goal of this post is to show explicitly how this condensation occurs by computing the participation ratio
$$
Y_N(\beta)=\sum_{i=1}^{2^N} w_i^2
%= \frac{Z_N(2\beta)}{Z_N(\beta)^2}.
$$
where $w_i$ represents the Boltzmann weight corresponding to the energy level $E_i$
$$
w_i \equiv \mu_\beta(E_i) = \frac{e^{-\beta E_i}}{Z_N(\beta)} 
$$
Intuitively if the Boltzmann measure is equidistributed over all the energy levels we expect that $Y_N(\beta) = 2^{-N}$ on average, i.e. it vanishes exponentially fast with $N$. If instead $Y_N(\beta)$ remains finite in the large $N$ limit it means that the measure is focused only on a $O(1)$ number of energy levels. 


@@toc-title
Contents
@@
\toc

## The high temperature phase

When $\beta < \beta_c$ it is easy to compute the participation ratio by using the identity
$$
\label{eq::YN_partition_function}
Y_N(\beta) = \sum_{i=1}^{2^N} \frac{e^{- 2\beta E_i}}{Z_N^2(\beta)} = \frac{Z_N(2\beta)}{Z_N(\beta)^2}
$$
which expresses it in terms of the partition function found in the previous [post](/blog/2025/REM1) which we report here for convenience
$$
Z_N(\beta) = e^{-N\beta f_N(\beta)} = \begin{cases}
e^{N \left(\ln 2 + \frac{\beta^2 J}{2} \right)} &  \; \mathrm{if} \;\beta< \beta_c\\
e^{N \beta \sqrt{2 J \ln 2} } & \; \mathrm{if} \; \beta \ge \beta_c \,.
\end{cases}
$$
We have two subcases, depending if $2\beta< \beta_c$ or $2\beta > \beta_c$, i.e. when the numerator in~\eqref{eq::YN_partition_function} is respectivelt in the high phase or the condensed phase. 

When $2\beta< \beta_c$ we have
$$
Y_N(\beta) = \frac{e^{N \left(\ln 2 + 2\beta^2 J \right)}}{e^{N \left(2\ln 2 + \beta^2 J \right)}} = e^{-N\left( \ln 2 - \beta^2 J \right)} < e^{-N\left( \ln 2 - \frac{\beta_c^2 J}{4} \right)} = e^ {- \frac{N}{2}\ln 2 } \to 0
$$
whereas if $\beta< \beta_c$ but $2\beta > \beta_c$ the numerator in~\eqref{eq::YN_partition_function} is in the low temperature phase
$$
Y_N(\beta) = \frac{e^{2 N \beta \sqrt{2J\ln 2}}}{e^{N \left(2\ln 2 + \beta^2 J \right)}} = e^{-N J \left( \beta_c^2 + \beta^2 - 2\beta \beta_c\right)} = e^{-N J \left( \beta - \beta_c\right)^2} \to 0
$$
Notice how the two expressions match at $\beta = \frac{\beta_c}{2}$. 

## The condensed phase: extreme energies and the distribution of Boltzmann weights

In order to derive the expression of the participation ratio in the low temperature phase, we need to better control the energy fluctuations near the lowest energy $E_i \simeq - N e_0$. 

To this end, we introduce a more general observable, denoted by $P(w)$, in terms of which the participation ratio can be easily expressed. $P(w)$ represents the *disorder-averaged density of Gibbs weights*
$$
P(w) = \mathbb{E}\left[\sum_{i} \delta(w-w_i) \right]
$$
meaning that $P(w)\, dw$ is the *expected number of energy levels having Gibbs weights in $[w, w+dw]$*. The (expected) participation ratio is expressed in terms of $P(w)$ as
$$
\label{eq::exp_part_ratio_p(w)}
\mathbb E[Y_N(\beta)]
=
E\left[\sum_i w_i^2\right] = \int_0^1 dw\, P(w)\,w^2 \,.
$$
In the following sections we will derive the expression of $P(w)$ in the condensed phase by controlling the energy fluctuations near the lowest available energy. 

#### Zooming near the lower edge of the spectrum

To study the relevant configurations at low temperature, we therefore write
$$
E_i=-Ne_0+ s_i,
$$
where $s_i=O(1)$ measures the small energy shift from the lower edge.

Using the Gaussian density of energy levels
$$
\rho_N(E)=\frac{1}{\sqrt{2\pi NJ}}e^{-\frac{E^2}{2NJ}},
$$
one finds
$$
\label{eq::scaling_edge_rho}
\rho_N(-Ne_0+s) = \frac{1}{\sqrt{2\pi NJ}}e^{-\frac{(N e_0 - s)^2}{2NJ}} = \frac{1}{\sqrt{2\pi NJ}}e^{-\frac{N e_0^2}{2J} + \frac{e_0 s }{J}} \left(1 + O\left(\frac{1}{N}\right)\right) \simeq A_N 2^{-N} e^{\beta_c s},
$$
up to a prefactor $A_N= (2\pi N J)^{-1/2}$ that depends on $N$ but not on $s$. This tells us that with respect to the extreme energy $-N e_0$ the energy levels are exponentially distributed in the shift variable $s$; namely among the $2^N$ energy levels, the states with a shift $s=O(1)$ occur with a density proportional to $e^{\beta_c s}\,ds$.


#### The Gumbel law

We are now going to show that in the condensed phase the dominant states live at a distance $O(1)$ from the spectral edge $-N e_0$, and their statistics is governed by extreme-value theory.

Let $E$ be the minimum among the $2^N$ energies. The probability density of $E$ is given by the probability of sampling an energy level $E$ times the probability that all the other $2^N-1$ energy levels have an energy larger than $E$:
$$
\rho_{\min}(E)=2^N \rho_N(E)\left[\int^{\infty}_E dE' \rho_N(E')\right]^{2^N-1} = 2^N \rho_N(E) \, \left[1 - H\left(-\frac{E}{\sqrt{NJ}}\right)\right]^{2^N-1}.
$$
where $H(x) = \frac{1}{2} \, \mathrm{Erfc}\left(\frac{x}{\sqrt{2}}\right)$. Now we impose that the minimal energy scales as $E=-Ne_0+s$; we use the fact that
$$
H(x) \simeq \frac{e^{-x^2/2}}{\sqrt{2\pi} x} \qquad \mathrm{if} \; x \to \infty
$$ 
so we have
$$
 \begin{split}
\rho_{\min}(-N e_0 + s) &\simeq 2^N \rho_N(-N e_0 + s) \, \left[1 - \frac{e^{-\frac{(-N e_0 + s)^2}{2NJ}}}{\sqrt{2\pi N J} \beta_c}\right]^{2^N} \\ 
&= 2^N \rho_N(-N e_0 + s) \, \left[1 - \frac{1}{\beta_c}\rho_N(-N e_0 + s)\right]^{2^N}.
\end{split}
$$
Using the asymptotic form of the density near the edge \eqref{eq::scaling_edge_rho}, we therefore find
$$
\rho_{\min}(s) = A_N e^{\beta_c s} \exp\left[- \frac{A_N}{\beta_c} e^{\beta_c s}\right]
$$
Introducing the rescaled variable
$$
\eta=\beta_c s+\ln\left(\frac{A_N}{\beta_c}\right),
$$
one gets
$$
\rho_{\min}(\eta)= \rho_{\min}(s(\eta)) \left|\frac{ds}{d\eta}\right| = e^{\eta-e^\eta}.
$$
which is a Gumbel distribution. The Gumbel law appears in the REM because below $T_c$ the Gibbs measure is controlled by the lowest energies, and the minimum of a large collection of independent Gaussian variables belongs to the Gumbel universality class of extreme-value theory (\cite{gumbel1958}).

#### From extreme energies to Gibbs weights

Now consider again the scaling $E_i=-Ne_0+u_i$ but directly applied to the Gibbs weight $ w_i=\frac{e^{-\beta E_i}}{Z_N(\beta)}$. As I will show here this will allow us to derive the distribution of Gibbs weights $P(w)$

Factoring out the common contribution $e^{\beta N e_0}$, we can write
$$
w_i=\frac{e^{-\beta s_i}}{e^{-\beta s_i}+Z_{\ne i}},
$$
where
$$
Z_{\ne i}=\sum_{j\neq i} e^{-\beta s_j}.
$$
We can find the distribution of $w_i$ by expressing it in terms of the distribution of the shift variables $s_i$ found previously in \eqref{eq::scaling_edge_rho}; solving for $s_i$ we find
$$
\label{eq::s(w)}
\beta s_i= - \ln\left(\frac{w_i Z_{\ne i}}{1-w_i}\right) = -\ln Z_{\ne i} - \ln w_i + \ln (1-w_i) \,.
$$
Since the $s_i$'s are distributed with intensity proportional to $e^{\beta_c s}$, we can transform variables from $s$ to $w$. Conditioned on the value of $Z_{\ne i}$, this gives
$$
P(w\mid Z_{\ne i}) \propto e^{\beta_c s}\left|\frac{ds}{dw}\right|.
$$
Using \eqref{eq::s(w)} we can write $e^{\beta_c s}=\left(\frac{w Z_{\ne i}}{1-w}\right)^{-\beta_c/\beta}$ and 
$$
\left|\frac{ds}{dw}\right|= T \left( \frac{1}{w} + \frac{1}{1-w} \right)= \frac{T}{w(1-w)},
$$
so that we get
$$
P(w\mid Z_{\ne i}) \propto (Z_{\ne i})^{-m} w^{-1-m}(1-w)^{-1+m}
$$
where
$$
m=\frac{T}{T_c}<1\,.
$$

Averaging over the random variable $Z_{\ne i}$ only affects the overall prefactor. We can compute it by simply noticing that $\sum_i w_i = 1$ so that 
$$
\int_0^1 dw \, w \, P(w) = 1\,.
$$
The normalization can be computed using the [beta function](https://en.wikipedia.org/wiki/Beta_function). The final result is
$$
\label{eq::P(w)}
\boxed{
    \begin{split}
P(w)&=\frac{w^{-m-1}(1-w)^{m-1}}{\Gamma(m)\Gamma(1-m)}\,,\\
m&=\frac{T}{T_c} \,.
\end{split}
}
$$
This is the main result. The expression of $P(w)$ gives a very concrete picture of the low-temperature phase of the REM. First, notice that $P(w)$ diverges near $w\to 0$ as
$$
P(w)\sim w^{-m-1}\quad (w\to 0)
$$
so that 
$$
\int_0^{\epsilon} d w\, P(w) \sim \int d w\, w^{-m-1} = \infty.
$$
i.e. there are infinitely many configurations carrying extremely small Gibbs weights. At the same time, $P(w)$ diverges near $w=1$ as
$$
P(w)\sim (1-w)^{m-1}\quad (w\to 1).
$$
so that the expected number of configurations with weight between $[1-\epsilon, 1]$ is, when $\epsilon$ is small equal to
$$
\int_{1-\epsilon}^{1} P(w) d w \propto \int dw \, (1-w)^{m-1} = \frac{\epsilon^m}{m}
$$
This shows that a single configuration may carry a finite fraction of the whole Gibbs measure. This is the real signature of condensation: below $T_c$, the measure is no longer democratically spread over exponentially many states, but becomes strongly uneven, with a few configurations carrying a macroscopic fraction of the total Boltzmann measure and many others contributing only weakly. In other words, below $T_c$ the measure is dominated by rare states sitting close to the lower edge of the energy spectrum.

These two singular behaviors therefore encode the full geometry of the frozen phase: a small number of dominant low-energy configurations coexist with a very large background of increasingly less important ones. The fact that $p(w)$ is not integrable is not a pathology, but precisely the mathematical expression of this accumulation of arbitrarily small Gibbs weights. This is precisely the signature of condensation: the Gibbs measure is not evenly spread, but concentrated on a sparse random set of configurations.



#### Recovering the participation ratio

Plugging in the explicit form of $P(w)$ in equation \eqref{eq::exp_part_ratio_p(w)} we can recover the expected participatio ratio as
$$
\mathbb E[Y_N(\beta)]
=
\frac{1}{\Gamma(m)\Gamma(1-m)}
\int_0^1 dw\, w^{1-m}(1-w)^{m-1}.
$$
Using again the [Beta integrals](https://en.wikipedia.org/wiki/Beta_function) the integral can be performed yielding to
$$
\mathbb E[Y_N(\beta)] = 1-m = 1-\frac{T}{T_c}.
$$
We have so far found the participatio ratio to be
$$
\lim_{N\to\infty} Y_N(\beta) =
\begin{cases}
0, & T>T_c,\\
1-\frac{T}{T_c}, & T < T_c
\end{cases}
$$
A finite participation ratio in the condensed phase signals that the Gibbs measure has condensed onto a small number of relevant low-energy states. This is why the low-temperature phase of the REM provides perhaps the simplest example of a condensed, or ''*glassy*'', Boltzmann measure.

The REM is special because its energies are independent, but the same structure survives in more complicated mean-field models where the energies are correlated. In particular, in the one-step replica-symmetry-breaking phase of $p$-spin models, one recovers the same distribution $P(w)$ we have derived here[^1]. The REM thus captures, in its simplest possible form, a structure that persists far beyond the independent-energy setting.

[^1]: For the readers expert in replica theory, the parameter $m$ appearing in the expression \eqref{eq::P(w)} of $P(w)$ can be identified with the Parisi's breaking parameter. 


@@notoc
References
@@

[1] \biblabel{derrida1980}{Derrida (1980)} Derrida, Bernard, "Random-energy model: Limit of a family of disordered models", Physical Review Letters 45.2 (1980): 79.

[2] \biblabel{gumbel1958}{Gumbel (1958)} Gumbel, E. "Statistics of extremes (New York: Columbia Univ. press)." (1958).

{{ blogcomments }}

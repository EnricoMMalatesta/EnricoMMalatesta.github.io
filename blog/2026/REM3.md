[**Could be worse**](https://www.youtube.com/watch?v=GPfKEhFh9ww) - Random topics in statistical physics, spin glasses, optimization and average case hardness

# The Random Energy Model via replica method


@def title = "The Random Energy Model -- Part III"
@def authors = "E. Malatesta"
@def published = "15 March 2025"
@def pt_lang = false
@def rss_pubdate = Date(2026, 3, 15)
@def rss = "REM3"
@def rss_description = """REM part 3"""

{{ published }} | **{{ authors }}**


In a previous [post](/blog/2025/REM1), we solved the random energy model (REM) by computing the microcanonical entropy, applying Markov’s inequality, and then performing a Legendre transform to obtain the canonical entropy. In this post, we show how the same results can be derived using the replica method, a heuristic yet powerful technique rooted in statistical physics. 


Applying the replica method to the REM is not merely a matter of rederiving known results by a different (and non-rigorous) route. It also gives access to additional information, such as the order parameters and their physical interpretation, and, more importantly, provides a framework that extends to many disordered models for which no equally direct solution is available. Historically, this application of the replica method also provided an early check on the method’s own validity.


@@toc-title
Contents
@@
\toc

## Main definitions

$$
Z^n = \sum_{i_1, \dots, i_n = 1}^{2^N} e^{- \beta \sum_{a = 1}^n E_{i_a}} = \sum_{i_1, \dots, i_n = 1}^{2^N} e^{- \beta \sum_k E_k \sum_{a = 1}^n \delta_{k i_a}}
$$

$$
\mathbb{E}\left[ Z^n \right] = \sum_{i_1, \dots, i_n = 1}^{2^N} \prod_{k} e^{- \frac{N J \beta^2}{2} \left(\sum_{a = 1}^n \delta_{k i_a}\right)^2} = \sum_{i_1, \dots, i_n = 1}^{2^N} e^{- \frac{N J \beta^2}{2} \sum_{a, b = 1}^n \delta_{i_a i_b} }
$$


Imposing
$$
q_{ab} = \delta_{i_a i_b}
$$
Denote by $\mathcal{M}_n$ the set of symmetric matrices having diagonal 1 and elements either 0 or 1. We can therefore rewrite the sum over the configurations $i_1\,, \dots \,, i_n \in [2^N]$ as the sum over the averaged replicated partition function as the sum over $q\in \mathcal{M}_n$ as follows
$$
\mathbb{E}\left[ Z^n \right] = \sum_{q \in \mathcal{M}_n} \mathcal{N}_N(q) \, e^{- \frac{N J \beta^2}{2} \sum_{a, b = 1}^n q_{ab} }
$$
where
$$
\mathcal{N}_N(q) \equiv \sum_{i_1, \dots, i_n = 1}^{2^N} \prod_{a<b}\delta(q_{ab} - \delta_{i_a i_b})
$$
counts the possible arrangements of indices $i_1\,, \dots \,, i_n \in [2^N]$ such that $q_{ab} = \delta_{i_a i_b}$. 


@@notoc
References
@@

[1] \biblabel{derrida1980}{Derrida (1980)} Derrida, Bernard, "Random-energy model: Limit of a family of disordered models", Physical Review Letters 45.2 (1980): 79.

{{ blogcomments }}

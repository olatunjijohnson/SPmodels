# Bounded-memory point prediction

## Purpose and exactness

A prediction grid may contain many more rows than the fitted dataset. Ordinary
prediction forms a training-by-new cross-covariance matrix. For \(n\) training
rows and \(m\) prediction rows this has \(nm\) entries, even though pointwise
means and marginal standard errors do not require covariance between distinct
prediction rows.

`sp_predict_chunked()` partitions the new rows and calls the same fitted-model
prediction equations on each partition. It changes memory scheduling only. It
does not approximate distances, covariance, latent conditioning or nonlinear
response integration.

## Partition and equivalence

Let the ordered new-row indices be partitioned into disjoint consecutive sets

\[
I_1\cup\cdots\cup I_J=\{1,\ldots,m\},\qquad
I_j\cap I_k=\varnothing\quad(j\ne k),                         \tag{CH-01}
\]

with at most \(c\) rows per set. For any set \(I_j\), the ordinary Gaussian
conditional mean is

\[
\widehat y_{I_j}=X_{I_j}\widehat\beta+o_{I_j}+
C_{I_j}^{\mathsf T}V^{-1}(y-o-X\widehat\beta),               \tag{CH-02}
\]

and its pointwise variances are the diagonal entries of equation Q-10. No
column of \(C_{I_j}\) enters the mean or marginal variance of a row outside
\(I_j\). Concatenating results in original order therefore gives

\[
(\widehat y_{I_1}^{\mathsf T},\ldots,
\widehat y_{I_J}^{\mathsf T})^{\mathsf T}=\widehat y_*        \tag{CH-03}
\]

and the same equality holds for the diagonal predictive variances.

For binomial and Poisson models, equations CH-01--CH-03 apply first to the
approximate latent mean and variance in Chapter 5. Logistic-normal quadrature
or lognormal moments are rowwise transformations of those quantities, so the
response means and marginal standard errors are also unchanged.

## Memory and limitations

The unchunked cross covariance needs \(O(nm)\) numeric storage. A chunk needs
at most \(O(nc)\), plus output vectors of length \(m\). Computation still
evaluates all \(nm\) distances and therefore has the same leading arithmetic
cost. Chunking does not make dense model fitting scalable in \(n\).

The method must not be used as a substitute for joint prediction. Covariance
between rows in different chunks is not constructed. Regional targets,
simultaneous excursions and dependent grid simulations must continue to use
`sp_predict_draws()` or a future explicitly validated scalable joint method.

Known and marginal unseen grouping levels retain ordinary pointwise behavior.
Rows sharing a new level across chunks would have joint covariance, but this
does not alter either row's marginal mean or variance; that cross-row
covariance is simply outside this function's declared output.

## Algorithm and checks

1. Validate the fitted object, non-empty prediction data and positive chunk
   size.
2. Form consecutive index sets by \(\lceil i/c\rceil\).
3. Call `predict()` with identical target, scale and unseen-level policy for
   each set.
4. Concatenate means and standard errors in original row order and restore row
   identifiers.
5. Record the requested chunk size and number of chunks.

`tests/test-scalable-prediction.R` compares chunk sizes that do not divide the
row count against one-shot Gaussian and binomial predictions, including
observation uncertainty, response integration, fixed correlated grouped
slopes and a shared unseen level. Equality is checked numerically, along with
ordering, metadata and invalid-input failures.

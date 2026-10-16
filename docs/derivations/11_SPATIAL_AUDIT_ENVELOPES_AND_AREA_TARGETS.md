# Spatial audit utilities, variogram envelopes and area targets

## Purpose and prerequisites

This chapter derives and specifies four supporting operations: auditable
distance summaries, projected UTM selection, Monte Carlo variogram envelopes,
and exact fractional-area aggregation. These operations do not introduce a new
response model. They make spatial inputs and model checks explicit enough for
research reports and MSc teaching.

The reader should know Euclidean distance, empirical quantiles, coordinate
reference systems, Monte Carlo simulation and elementary polygon area.

## Pairwise distance summary

For two-dimensional projected coordinates
\(\boldsymbol s_i=(s_{i1},s_{i2})\), the unordered-pair distance is

\[
d_{ij}=\sqrt{(s_{i1}-s_{j1})^2+(s_{i2}-s_{j2})^2},
\qquad 1\leq j<i\leq n.                                      \tag{D-01}
\]

There are \(n(n-1)/2\) pairs. A repeated location is defined by exact equality
of both supplied coordinate values, not by an undocumented proximity
tolerance. Its distance is zero. The caller must choose one of three policies:

1. include zero-distance pairs in every summary;
2. exclude them while retaining their count; or
3. stop with an error.

The returned mean and sample standard deviation use the retained distances.
Quantiles use R's type-8 definition. If \(p\in[0,1]\), set

\[
h=(m+1/3)p+1/3,                                               \tag{D-02}
\]

where \(m\) is the retained pair count, and interpolate between the adjacent
ordered distances at position \(h\), with the usual endpoint truncation. The
coordinate unit is never guessed for an ordinary data frame; it is inherited
only from a fitted model or projected `sf` CRS.

Implementation: `sp_distance_summary()` in `R/geometry-utils.R`.

## UTM coordinate-reference selection

### Ordinary zone rule

After transforming all geometry vertices to WGS84 longitude/latitude, the
ordinary six-degree zone for longitude \(\lambda\) is

\[
z(\lambda)=\min\{60,\max[1,
\lfloor(\lambda+180)/6\rfloor+1]\}.                           \tag{D-03}
\]

Thus an exact interior six-degree boundary belongs to the eastern zone and
longitude 180 degrees is capped at zone 60.

### Standard high-latitude exceptions

The implemented rule assigns zone 32 for 56--64 degrees north and 3--12
degrees east. From 72--84 degrees north, it assigns zones 31, 33, 35 and 37 to
longitude bands 0--9, 9--21, 21--33 and 33--42 degrees east respectively.
These rules prevent the ordinary formula from selecting the narrow zones that
are replaced in Norway and Svalbard.

UTM selection is rejected south of 80 degrees south or north of 84 degrees
north. A longitude span above 180 degrees is treated as an antimeridian
crossing and rejected because an arithmetic bounding-box centre is then
misleading.

Every coordinate is assigned a zone. Multiple observed zones cause an error
by default. With the explicit centre policy, the zone is instead evaluated at

\[
\lambda_c=(\lambda_{\min}+\lambda_{\max})/2,
\qquad
\varphi_c=(\varphi_{\min}+\varphi_{\max})/2.                  \tag{D-04}
\]

Likewise, a geometry spanning the equator errors by default. The explicit
centre policy uses the sign of \(\varphi_c\), with zero assigned north. The
WGS84 EPSG code is

\[
\operatorname{EPSG}=
\begin{cases}
32600+z,&\text{northern hemisphere},\\
32700+z,&\text{southern hemisphere}.
\end{cases}                                                    \tag{D-05}
\]

The result records all observed zones, centre, policies and equator status so
the choice remains auditable. It selects a CRS; the user still performs the
visible transformation with `sf::st_transform()`.

Implementation: `sp_utm_crs()` in `R/geometry-utils.R`.

## Residual-variogram Monte Carlo envelope

Let \(\widehat\gamma_k\) be the observed mean semivariance in distance bin
\(k\), defined in derivation 03. SPmodels supports two null-generating rules.

### Permutation null

For the no-location-association null, keep coordinates and distance bins fixed
and randomly permute the selected residual vector:

\[
\boldsymbol e^{(b)}=\pi_b(\boldsymbol e),
\qquad b=1,\ldots,B.                                          \tag{D-06}
\]

This preserves the exact observed residual distribution but breaks its
association with location. It addresses whether residual spatial arrangement
is more structured than expected under exchangeability.

### Fitted-model null

The parametric rule draws

\[
\boldsymbol Y^{(b)}\sim
N(\boldsymbol X\widehat{\boldsymbol\beta}+\boldsymbol o,
  \widehat{\boldsymbol V})                                    \tag{D-07}
\]

with fitted parameters held fixed. Mean residuals subtract the fitted mean.
Conditional residuals additionally subtract the fixed-parameter spatial
conditional adjustment

\[
\widehat{\boldsymbol K}\widehat{\boldsymbol V}^{-1}
(\boldsymbol Y^{(b)}-\boldsymbol X\widehat{\boldsymbol\beta}-
\boldsymbol o).                                                \tag{D-08}
\]

This checks behaviour under the fitted model but does not include parameter
re-estimation uncertainty.

For every simulation, use the same pairs and bins to obtain
\(\widehat\gamma_k^{(b)}\). A pointwise level \(1-\alpha\) envelope is

\[
[q_{\alpha/2,k},q_{1-\alpha/2,k}],                            \tag{D-09}
\]

where the endpoints are type-8 empirical quantiles of the simulated bin
values. Pointwise limits are not a simultaneous confidence band.

For a global curve comparison, let \(\bar\gamma_k\) and \(s_k\) be the
simulated mean and standard deviation in every populated, variable bin. Define

\[
T_{\mathrm{obs}}=\sum_k
\left(\frac{\widehat\gamma_k-\bar\gamma_k}{s_k}\right)^2,
\qquad
T_b=\sum_k
\left(\frac{\widehat\gamma_k^{(b)}-\bar\gamma_k}{s_k}\right)^2.\tag{D-10}
\]

The finite-simulation p-value is

\[
p=\frac{1+\sum_{b=1}^B I(T_b\geq T_{\mathrm{obs}})}{B+1}.     \tag{D-11}
\]

The added one prevents a reported p-value of zero. A seed is local to the call,
and raw simulated curves are returned so every limit and statistic can be
reconstructed.

Implementation: `sp_variogram_envelope()` in `R/assessment.R`.

## Exact fractional-area polygon targets

Let prediction cell \(C_i\) and target region \(A_g\) be projected polygons.
Their exact intersection area is

\[
a_{ig}=|C_i\cap A_g|,\qquad
f_{ig}=a_{ig}/|C_i|.                                          \tag{D-12}
\]

The fraction \(f_{ig}\) is retained for auditing. With optional non-negative
cell density \(d_i\), define raw mass weight

\[
w_{ig}=a_{ig}d_i;                                             \tag{D-13}
\]

without a density, \(d_i=1\). As in derivation 06, weights are normalized
within each region and the target is calculated inside every joint draw:

\[
T_g^{(b)}=
\frac{\sum_i w_{ig}Y_i^{(b)}}{\sum_i w_{ig}}.                 \tag{D-14}
\]

A boundary-crossing cell can therefore contribute to several regions without
duplicating its full weight. Target polygons must not overlap, all geometry
must share one known projected CRS, and intersections with zero area are
discarded. The caller explicitly chooses whether unmatched cells and declared
regions with no positive-area intersection cause an error or are recorded and
dropped. Dropping is never silent.

`sp_target_geometry()` then joins the already calculated, unrounded target
table to point or polygon `sf` geometry by exact unique identifiers or an
explicit row-order contract. Its plot method is only a view of those values.

Implementation: `sp_polygon_area_weights()` and `sp_target_geometry()` in
`R/spatial-targets.R`.

## Code and validation map

| Equations | Implementation | Independent checks |
| --- | --- | --- |
| D-01--D-02 | `sp_distance_summary()` | Complete four-location distance enumeration and all duplicate policies |
| D-03--D-05 | `sp_utm_crs()` | North/south EPSG codes, six-degree edge, longitude 180, Norway, Svalbard, multi-zone, equator, polar and antimeridian cases |
| D-06--D-11 | `sp_variogram_envelope()` | Exact seeded reconstruction of every permutation curve, standardized statistics and p-value; parametric reproducibility |
| D-12--D-14 | `sp_polygon_area_weights()` | Direct split-cell areas/fractions, density weighting, draw-level aggregation, unmatched cell and empty-region policies |
| Join contract | `sp_target_geometry()` | Reordered identifier match, attached numerical values and mismatch failure |

The relevant tests are `tests/test-geometry-utils.R`,
`tests/test-variogram-envelope.R` and `tests/test-spatial-targets.R`.

## Limitations

UTM is a practical local projection choice, not a universal answer for a study
spanning zones, the equator, polar areas or the antimeridian. Centre policies
must be scientifically justified. Parametric variogram envelopes condition on
fitted parameters. Fractional-area density weighting assumes the supplied
density is constant within each prediction cell. The target plot uses base
`sf` display and does not choose palettes or cartographic annotations.

# upsis

`upsis` updates fitted `brms` models after modest additions or removals of data
using Pareto-smoothed importance sampling (PSIS). Updating with PSIS gives an
effecient way to "upsize" the dataset used for model fitting, without a full
re-fit.

The package computes the change in log likelihood implied by the requested data
update, smooths the resulting importance weights with `loo::psis()`, and
stratified-resamples posterior draws back into the fitted `brmsfit` object. The
updated model can then be used with familiar `brms` post-processing tools.

```r
library(brms)
library(upsis)

fit <- brm(y ~ x, data = old_data)
result <- upsis(fit, data_add = new_data)

summary(result$updated_model)
result$pareto_k
```

## When to refit instead

`upsis` is an approximation to refitting the model on the changed data. It is
best suited to small, weakly influential updates that do not change the model
structure.

Prefer a full refit when:

- the Pareto `k` diagnostic is high;
- the new data are highly influential or much larger than the original data;
- the update changes data-dependent priors, spline knots, basis expansions, or
  other model components chosen during the original fit;
- the update introduces unsupported new grouping levels or otherwise requires
  new parameters.

## Development

Run the lightweight test suite with:

```r
devtools::test()
```

The optional `brms` integration test is skipped by default. To run it locally,
set `RUN_BRMS_TESTS=true` before testing.

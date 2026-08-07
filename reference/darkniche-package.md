# darkniche: Climatic Dark Diversity and Directional Niche Coverage

Which species are absent from a region although its climate is inside
their tolerance? `darkniche` answers that question end to end: it
estimates niches and regional climate availabilities on one shared
environmental lattice, scores climatic admissibility with a
*directional* metric family that expresses containment (a question
symmetric overlap statistics such as Schoener's *D* provably cannot
express), screens each pair against three analytical null models,
calibrates the admissibility threshold against documented introduced
ranges, and reports the resulting dark diversity with the two guardrails
the underlying studies found to be indispensable: an area adjustment on
the dark-diversity ratio, and a threshold-sweep diagnostic that reveals
when the direction of an aggregate null-model test is an artefact of the
threshold rather than a property of the region.

## Workflow

1.  [`dn_space()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_space.md),
    [`dn_density()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_density.md)
    – one lattice, one density per species and per region.

2.  [`dn_metrics()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_metrics.md)
    – the directional family (`DNC_q`, `DNC_auc`, `NSI_q`, `CFI`,
    `CFI_w`, `NOV`) alongside `D`, `I` and `BC`.

3.  [`dn_null_breadth()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_breadth.md),
    [`dn_null_shift()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_shift.md),
    [`dn_null_relocate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_relocate.md)
    – null screens;
    [`dn_ses()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_ses.md)
    standardises them.

4.  [`dn_match()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_match.md),
    [`dn_auc()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_auc.md),
    [`dn_delong()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_delong.md),
    [`dn_calibrate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_calibrate.md)
    – validation against documented introduced ranges.

5.  [`dn_dark()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_dark.md),
    [`dn_cddr()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_cddr.md),
    [`dn_report()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_report.md)
    – accounting.

6.  [`dn_aggregate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_aggregate.md),
    [`dn_sweep()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_sweep.md),
    [`dn_direction_census()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_direction_census.md)
    – region-level inference and its diagnostic.

## Guardrails

[`dn_aggregate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_aggregate.md)
suppresses the effect size, p-value and enrichment ratio whenever the
expected count falls below five;
[`dn_cddr()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_cddr.md)
returns the area-adjusted residual beside the raw ratio and documents
that raw CDDR must not rank regions;
[`dn_sweep()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_sweep.md)
prints a reversal warning from the object itself, computed only from
thresholds the test considers reliable.

## Example data

[`dn_data()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_data.md)
gives the published results for all 619 WWF terrestrial ecoregions, so
every example and test runs on real output.

## References

Rahimi, E. & Jung, C. (in review) Directional niche coverage and the
climatic dark diversity of Iranian ecoregions.

Rahimi, E. & Jung, C. (in review) Climatic dark diversity of the world's
terrestrial ecoregions and the threshold dependence of aggregate
null-model tests.

## See also

Useful links:

- <https://github.com/ehsanrahimi666/darkniche>

- <https://ehsanrahimi666.github.io/darkniche/>

- Report bugs at <https://github.com/ehsanrahimi666/darkniche/issues>

## Author

**Maintainer**: Ehsan Rahimi <ehsanrahimi666@gmail.com>

Authors:

- Chuleui Jung

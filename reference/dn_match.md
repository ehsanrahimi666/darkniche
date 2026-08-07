# Realm- and area-matched negatives for introduced-range validation

For each positive (region a species actually reached), samples `ratio`
negatives from unreached regions in the same biogeographic realm with
the nearest log-areas – matching without which negatives are dominated
by climatically implausible regions and every metric rejects them
equally.

## Usage

``` r
dn_match(positives, regions, reached, ratio = 5L, seed = 1L)
```

## Arguments

- positives:

  data.frame with `species`, `region`.

- regions:

  data.frame with `region`, `realm`, `area`.

- reached:

  list mapping species to reached regions.

- ratio:

  negatives per positive (default 5).

- seed:

  RNG seed.

## Value

data.frame `species`, `region`, `label`.

## Examples

``` r
regs <- data.frame(region = paste0("r", 1:8), realm = rep(c("PA", "NT"), 4),
                   area = 10^runif(8, 4, 6))
pos <- data.frame(species = "sp1", region = "r1")
dn_match(pos, regs, reached = list(sp1 = "r1"), ratio = 3)
#>   species region label
#> 1     sp1     r1  TRUE
#> 2     sp1     r5 FALSE
#> 3     sp1     r3 FALSE
#> 4     sp1     r7 FALSE
```

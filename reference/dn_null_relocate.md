# Geographic relocation null (N3)

Null coverage obtained by stamping the species' occupied-cell footprint
at random positions on the geographic grid and recomputing coverage from
the relocated range. Requires the geographic cell adjacency of the
analysis grid; supply `project` to map relocated geographic cells to
environmental scores.

## Usage

``` r
dn_null_relocate(
  cells,
  grid_dim,
  valid,
  project,
  sp,
  Ze,
  q = 0.95,
  statistic = c("auc", "q"),
  B = 199L,
  seed = 1L
)
```

## Arguments

- cells:

  integer matrix (n x 2) of occupied geographic cell indices.

- grid_dim:

  geographic grid dimensions `c(nrow, ncol)`.

- valid:

  logical matrix marking terrestrial cells.

- project:

  function taking an (n x 2) cell-index matrix and returning (n x 2)
  environmental scores (NA rows for cells without climate).

- sp:

  a `dn_space`.

- Ze:

  (k x P) region densities.

- q, statistic, B, seed:

  as in
  [`dn_null_shift()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_shift.md).

## Value

(k x B) matrix of null coverage values (NA where a relocation could not
be placed).

## Examples

``` r
# Requires a geographic grid and a projection function; sketch:
if (FALSE) { # \dontrun{
nul <- dn_null_relocate(cells, grid_dim = dim(valid), valid = valid,
                        project = my_lookup, sp = sp, Ze = Ze, B = 99)
} # }
```

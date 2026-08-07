# Export a table as CSV, Markdown or LaTeX

Export a table as CSV, Markdown or LaTeX

## Usage

``` r
dn_table_export(x, file, digits = 3, caption = NULL)
```

## Arguments

- x:

  data.frame.

- file:

  output path; format from extension (`.csv`, `.md`, `.tex`).

- digits:

  rounding for numeric columns (default 3).

- caption:

  optional caption (md/tex).

## Value

invisibly, the path.

## Examples

``` r
tf <- tempfile(fileext = ".md")
dn_table_export(data.frame(a = 1:2, b = c(0.123456, 2)), tf,
                caption = "demo")
cat(readLines(tf), sep = "\n")
#> **demo**
#> 
#> | a | b |
#> | --- | --- |
#> | 1 | 0.123 |
#> | 2 | 2 |
```

# WIDclocks

Easy calculation of WID clocks (Barrett J., Herzog C. & Kim. Y.-N. et al., 2021).

**Please note: this code is available for research use only. Commercial use of the code or any data related to it is prohibited.**

THE WID_clocks function uses a methylation beta matrix and returns the WID general, epithelial and immune clocks.
It also returns the WID-relative-epithelial-age (WID_rea) and WID-relative-immune-age (WID_ria) (∆WID epithelial or WID immune and WID general clock).

# Installation

```
devtools::install_github("chiaraherzog/WIDclocks")
```

# Use

```
res <- WID_clocks(beta_matrix)
```

## Cell composition

The package requires previous installation of the EpiDISH package to infer the cellular composition of the samples. In case EpiDISH is not automatically installed, please use the following code:

```{r echo=TRUE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EpiDISH")
```




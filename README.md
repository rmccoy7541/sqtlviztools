# sQTLviztools

An R package and associated scripts for visualisation splicing quantitative trait loci (sQTLs).

## Installation

1. clone the sqtlviztools repo

```
git clone https://github.com/jackhump/sqtlviztools.git 
```

2. install the R package and any dependencies

in R:
```{r}
install.packages(c(
"shiny", "DT", "shinycssloaders", "shinyjs", "dplyr", "ggplot2", "leafcutter", "reshape2", "gridExtra", "intervals", "foreach", "grid", "gtable", "ggrepel", "ggbeeswarm", "stringr"
))

```
on the command line: 

```
R CMD INSTALL sqtlviztools
```

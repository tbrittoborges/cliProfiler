# cliProfiler
Author: You Zhou, Kathi Zarnack    

---

# Introduction
Cross-linking immunoprecipitation (CLIP) is a technique that combines UV 
cross-linking and immunoprecipitation to analyse protein-RNA interactions or to 
pinpoint RNA modifications (e.g. m6A). CLIP-based methods, such as iCLIP and 
eCLIP, allow precise mapping of RNA modification sites or RNA-binding protein 
(RBP) binding sites on a genome-wide scale. These techniques help us to unravel 
post-transcriptional regulatory networks. In order to make the visualization of 
CLIP data easier, we develop `r Biocpkg("cliProfiler")` package. The 
`r Biocpkg("cliProfiler")` includes seven functions which allow users easily 
make different profile plots.

The `r Biocpkg("cliProfiler")` package is available at
[https://bioconductor.org](https://bioconductor.org) and can be
installed via `BiocManager::install`:

```{r BiocManager, eval=FALSE}
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("cliProfiler")
```

A package only needs to be installed once. Load the package into an
R session with

```{r initialize, results="hide", warning=FALSE, message=FALSE}
library(cliProfiler)
```

# How to use it
Documentation (vignette and user manual) is available at the **cliProfiler** 
Bioconductor landing page at http://bioconductor.org/packages/cliProfiler.

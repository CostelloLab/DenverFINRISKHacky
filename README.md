# DenverFINRISKHacky
Hackathon for the DREAM FINRISK Challenge cardiovascular event prediction from microbiome 

Placeholder README for the Hackathon taking place in Denver 14th to 18th November 2022.

## Things to consider:

### Installed software on your computer (preferably laptop):

- R (4.2.1 or preferably newer)
- R Studio suggested if you want a bit of assistance
- Highly recommend git, so we can merge and share code easily
- GitHub (where you're currently at) will serve as our main platform

### R packages & functions

#### Installation from CRAN or BioConductor

Install packages from the terminal with `install.packages("pckgName")` (replace `pckgName` with the exact name of the package you want to install.
Note that this only works for Central R Archive Network packages (CRAN). 

For Bioconductor you need to use:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
```

#### Recommended R packages

Packages (and for what purpose they could be used):

- `ggplot2` (for visualizations beyond base visualizations)
- `ComplexHeatmap` (for advanced heatmap visualizations)
- `glmnet` (for regularized regression with L1/L2 norms)
- ...

#### R functions to check out

name of R `function`; description of the purpose (if it needs a non-base package, `package name` indicated here)

- `write.table`; for writing output
- `read.table`; for reading output
- `setwd`; for setting current working directory, notice that R understands this as your root directory
- `|>`; R's native pipe operator (R >= 4.3.0)
- `plot`; an ordinary x-axis vs. y-axis plotting
- `boxplot`; a boxplot
- `ls`; list objects in current R workspace
- `<-`; placement operator (recommended for using when not talking about function parameters)
- `=`; placement operator or function parameters (recommended to use only with function parameters)
- ... 

### Advanced stuff to consider

- Learning how to use Singularity containers: https://docs.sylabs.io/guides/3.5/user-guide/introduction.html (the DREAM Challenge comes with an example container)


### Literature regarding the topic
For example, reading 
- Salosensaari A, Laitinen V, Havulinna AS et al., Taxonomic signatures of cause-specific mortality risk in human gut microbiome. Nat Commun 12, 2671 (2021). https://doi.org/10.1038/s41467-021-22962-y 
would be useful.


### Synapse and the DREAM

Registration on Synapse platform is highly encourated:


Overview to the Challenge phases
[https://www.synapse.org/#!Synapse:syn27130803/wiki/619280](https://www.synapse.org/#!Synapse:syn27130803/wiki/619280)

Data description is essential for understanding what we're trying to model:
[https://www.synapse.org/#!Synapse:syn27130803/wiki/619274](https://www.synapse.org/#!Synapse:syn27130803/wiki/619274)

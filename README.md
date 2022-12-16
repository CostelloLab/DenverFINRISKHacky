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

### Objects in train_test.RData

- train_x / test_x : rank inverse normal distributed counts
- train_raw / test_raw : read counts before above normalization
- train_p / test_p : phenotype (clinical) info for train and test sets respectively
- train_y / test_y : Surv-response vector for train and test sets respectively 

### Short citations / insight

- "Valtimokovettuma tautiin eli ateroskleroosiin ja siihen liittyviin syd�n- ja verisuoniongelmiin on mikrobiomilla tutkittuja yhteyksi�. Asteikolla 1-10 mitattuna, mit� enemm�n suojaavia tekij�it� l�ydet��n mikrobiomista, sit� korkeampi on suojaava vaikutus n�ihin sairauksiin.   Jotkut bakteerit, kuten Erysipelotrichaceae tuottaa trimetyyliamiinia (TMA) koliinista, jota l�ytyy kananmunan keltuaisesta, sek� L-kartiniinia, jota saadaan punaisesta lihasta ja kalasta. TMA hapetetaan maksassa trimetyyliamiini-typpioksidiksi TMAO, johon liittyy kohonnut riski sairastua ateroskleroosiin ja muihin syd�n- ja verisuonitauteihin."
(https://www.riskigeeni.net/mikrobiomi.html)

-> Some bacteria such as Erysipelotrichaeae produce trimethylaminine (TMA) from coline and L-cartinine; TMA is oxidated in the liver to TMAO, which involves a heightened risk to suffer from ateroscloreosis and other cardiovascular diseases

### Summary of literature findings



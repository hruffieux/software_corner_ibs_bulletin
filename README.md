## Software corner IBS bulletin on atlasqtl

### Data

The genotyping data are protected. We are therefore using a synthetic dataset 
emulating the real data. The expression and replicated genotyping data are in 
[data/replicated_data.RData](https://github.com/hruffieux/software_corner_ibs_bulletin/blob/master/data/replicated_data.RData). The ready-to-use data with simulated genetic associations 
are in 
[data/prepared_data.RData](https://github.com/hruffieux/software_corner_ibs_bulletin/blob/master/data/prepared_data.RData).

The file [data/prepared_data.RData](https://github.com/hruffieux/software_corner_ibs_bulletin/blob/master/data/prepared_data.RData) 
is obtained by updating the real expression data to simulate genetic 
associations between the synthetic genotyping data and the transcript levels. 
This last step is obtained by running the R file 
[scripts/prepare_data.R](https://github.com/hruffieux/software_corner_ibs_bulletin/blob/master/scripts/prepare_data.R) 
after installing the R package **echoseq**:

```R
if(!require(remotes)) install.packages("remotes")
remotes::install_github("hruffieux/echoseq")
```

### Algorithm

The package **atlasqtl** used for the analysis may also be installed with:
 
```R
if(!require(remotes)) install.packages("remotes")
remotes::install_github("hruffieux/atlasqtl")
```

### eQTL analysis

The eQTL analysis can be run using the Rmarkdown script: 
[scripts/atlasqtl_software_corner.Rmd](https://github.com/hruffieux/software_corner_ibs_bulletin/blob/master/scripts/atlasqtl_software_corner.Rmd). 
This file also acts as a simple tutorial for the use of **atlasqtl**.

### Issues

To report an issue, please use the 
[issue tracker](https://github.com/hruffieux/software_corner_ibs_bulletin/issues) 
at github.com.
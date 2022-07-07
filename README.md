## Software corner IBS bulletin: atlasqtl – an R package for variable selection in sparse regression with hierarchically-related responses

### Data

The genotyping data are protected. We are therefore using a synthetic dataset emulating the real data. 
The expression and replicated genotyping data are in 
[data/replicated_data.RData](https://github.com/hruffieux/bayesian_variable_selection_book_chapter/blob/master/data/replicated_data.RData). The ready-to-use data with simulated genetic associations are in 
[data/prepared_data.RData](https://github.com/hruffieux/bayesian_variable_selection_book_chapter/blob/master/data/prepared_data.RData).

**Important note:** these are large files which are stored using Git Large File Storage. To clone these
files along with the repository, please install Git LFS, e.g., using Homebrew:

``` bash
brew install git-lfs
```

and, within the cloned repository, initialise it for your account by using:

``` bash
git lfs install
```

and retrieve all large files:

``` bash
git lfs pull
```

(Alternatively, since there are only two such files, they can be downloaded manually from the Github interface.)

The file [data/prepared_data.RData](https://github.com/hruffieux/software_corner_ibs_bulletin/blob/master/data/prepared_data.RData) 
is obtained by updating the real expression data to simulate genetic associations between the synthetic 
genotyping data and the transcript levels. This last step is obtained by running the R file 
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
[scripts/atlasqtl_example.Rmd](https://github.com/hruffieux/software_corner_ibs_bulletin/blob/master/scripts/atlasqtl_example.Rmd). 
This file also provides step-by-step guidance for the use and settings of the **atlasqtl** for our example.

### Issues

To report an issue, please use the 
[issue tracker](https://github.com/hruffieux/software_corner_ibs_bulletin/issues) 
at github.com.
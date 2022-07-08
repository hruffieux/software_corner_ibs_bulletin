rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

data_dir <- file.path(CORE_DIR, "software_corner_ibs_bulletin/data/")
main_dir <- file.path(CORE_DIR, "software_corner_ibs_bulletin/scripts/")

setwd(main_dir)

RNGkind("L'Ecuyer-CMRG") # to ensure reproducibility when using parallel processes
my_seed <- 1
set.seed(my_seed)

## -----------------------------------------------------------------------------
## Please install the following R package: echoseq using:
## devtools::install_github("hruffieux/echoseq")
## -----------------------------------------------------------------------------

require(echoseq)

bool_save <- TRUE

## -------------------------------------------------------------
## --- This RData file was generated using genotyping data -----
## -- It can't be regenerated as we do not provide the real ----
## -- genotyping data but the expression data and synthetic ----
## -- genotyping data obtained by emulating the real can be ----
## -- downloaded from data/replicated_data.RData and loaded as 
## -- done below -
#
load(file.path(data_dir, "replicated_data.RData"))

p <- ncol(X) # genetic variants (single nucleotide polymorphisms, SNPs)
q <- ncol(Y_original) # B cell expression data

stopifnot(nrow(X) == nrow(Y_original)) # same number of samples for SNPs and expression levels
n <- nrow(X)
rownames(Y_original) <- rownames(X)

## simulates association pattern between the SNPs and the expression levels
#
p0 <- 5  # number of active SNPs
ind_p0 <- sort(sample(1:p, p0, replace = FALSE))

q0 <- 250  # number of expression levels under genetic control
ind_q0 <- sort(sample(1:q, q0, replace = FALSE))

# vec_prob_sh: vector of length p0 providing the probabilities with which each
# active SNP will be associated with an additional active expression level
vec_prob_sh <- rbeta(p0, shape1 = 1, shape2 = 2)

max_tot_pve <- 0.5 # maximum total proportion of outcome variance explained by the SNPs.

obj_snps <- convert_snps(X)
obj_expr <- convert_phenos(Y_original)

# add some noise so these propensities differ slightly depending on the condition
vec_prob_sh_type <- vec_prob_sh + rnorm(p0, sd = 0.1)

# make sure that between zero and one
vec_prob_sh_type[vec_prob_sh_type < .Machine$double.eps | vec_prob_sh_type > 1 - .Machine$double.eps] <- vec_prob_sh[vec_prob_sh_type < .Machine$double.eps | vec_prob_sh_type > 1 - .Machine$double.eps]

obj_data <- generate_dependence(obj_snps, obj_expr, ind_q0, ind_p0,
                                vec_prob_sh_type, family = "gaussian",
                                max_tot_pve = max_tot_pve)

stopifnot(all.equal(rownames(obj_data$snps), rownames(obj_data$phenos)))

X <- obj_data$snps
Y <- obj_data$phenos
pat <- obj_data$pat

if (bool_save) {
  save(X, Y, pat, file = file.path(data_dir, "prepared_data.RData"))
}

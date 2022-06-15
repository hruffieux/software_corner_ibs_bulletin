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
## -- downloaded from data/replicated_data.R and loaded as 
## -- done below -
#
load(file.path(data_dir, "replicated_data.RData"))

vec_type <- names(list_expr)

p <- unique(sapply(list_repl_snps, ncol))
stopifnot(length(p) == 1) # same number of of SNPs for all conditions

q <- unique(sapply(list_expr, ncol))
stopifnot(length(q) == 1) # same number of of levels for all conditions

## simulates association pattern between the SNPs and the expression levels
# same active SNPs and active levels for all conditions
# but a given active SNP may be associated with different active levels across the conditions
#
p0 <- 5  # number of active SNPs
ind_p0 <- sort(sample(1:p, p0, replace = FALSE))

q0 <- 250  # number of expression levels under genetic control
ind_q0 <- sort(sample(1:q, q0, replace = FALSE))

# vec_prob_sh: vector of length p0 providing the probabilities with which each
# active SNP will be associated with an additional active expression level
vec_prob_sh <- rbeta(p0, shape1 = 1, shape2 = 2)

max_tot_pve <- 0.5 # maximum total proportion of outcome variance explained by the SNPs.

list_data <- NULL
for (type in vec_type) {

  replicated_snps <- list_repl_snps[[type]]
  original_expr <- list_expr[[type]]

  n <- unique(c(nrow(replicated_snps), nrow(original_expr)))
  stopifnot(length(n) == 1) # same number of samples for SNPs and expression levels
  
  rownames(original_expr) <- rownames(replicated_snps)

  obj_snps <- convert_snps(replicated_snps)
  obj_expr <- convert_phenos(original_expr)

  # add some noise so these propensities differ slightly depending on the condition
  vec_prob_sh_type <- vec_prob_sh + rnorm(p0, sd = 0.1)
  
  # make sure that between zero and one
  vec_prob_sh_type[vec_prob_sh_type < .Machine$double.eps | vec_prob_sh_type > 1 - .Machine$double.eps] <- vec_prob_sh[vec_prob_sh_type < .Machine$double.eps | vec_prob_sh_type > 1 - .Machine$double.eps]
  
  obj_data <- generate_dependence(obj_snps, obj_expr, ind_q0, ind_p0,
                                  vec_prob_sh_type, family = "gaussian",
                                  max_tot_pve = max_tot_pve)

  stopifnot(all.equal(rownames(obj_data$snps), rownames(obj_data$phenos)))

  list_data <- append(list_data,
                      list(list("snps" = obj_data$snps, "expr" = obj_data$phenos)))
}
names(list_data) <- vec_type

# Check that same set of SNPs and expression levels across conditions.
list_snp_names <- lapply(list_data, function(ll_type) names(ll_type$snps))
stopifnot(length(unique(list_snp_names)) == 1)

list_expr_names <- sapply(list_data, function(ll_type) names(ll_type$expr))
stopifnot(length(unique(list_expr_names)) == 1)

if (bool_save) {
  save(list_data, file = file.path(data_dir, "prepared_data.RData"))
}

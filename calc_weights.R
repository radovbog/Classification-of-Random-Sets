# Calculates weight that should be given to ratio w.r.t. curvature

weight <- function(NdistC, NdistR) {
  # Input is Ndistance matrices NdistC and NdistR calculated using 
  # only curvature and only ratio, respectively
  
  # Number of rows in NdistC and NdistR
  nC <- nrow(NdistC)
  nR <- nrow(NdistR)
  
  # Compute the column means
  LC <- rowSums(NdistC) / (nC - 1)
  LR <- rowSums(NdistR) / (nR - 1)
  
  # Compute NC and NR
  NC <- mean(LC)
  NR <- mean(LR)
  
  # Compute variance components
  VarC <- sum((NdistC - NC)^2) / (nC * (nC - 1))
  VarR <- sum((NdistR - NR)^2) / (nR * (nR - 1))
  
  # Compute final weight
  w <- sqrt(VarC / VarR)
  
  return(w)
}


# ####################################################
# ###################### MAIN ########################
# ####################################################

# Capture the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Base directory containing the input data
#base_dir <- args[1]
base_dir <- "~/Desktop/"
#base_dir<-"C:/Users/Spravce/Desktop"

# The radius of the osculating circle
#radius <- args[2]
radius <- 5

# Number of runs
#n_runs <- args[2]
n_runs <- 50
#partial_run <- args[3]
#partial_run <- 1

# Number of realisations we consider
#n_reals <- args[4]
n_reals <- 200

# Realisations start indexing from 0
n_real <- n_reals - 1
#n_real <-2

# Considered processes
n_classes <- 3
names <- c("Boolean", "Cluster", "Repulsive")
classes <- c("B", "C", "R")

# Considered sample sizes
ss_names <- c("10", "20", "All")

weights <- matrix(0, ,nrow=n_runs, ncol=n_classes)
a_from_w <- matrix(0, ,nrow=n_runs, ncol=n_classes)
# Loop through different runs and set the seed for each run
for(n_run in 1:n_runs){
  #for(n_run in partial_run:partial_run){
  
  set.seed(n_run)
  
  for(ss in 1:3){
    sample_size <- ss_names[ss]
    
    # Load data 
    N.dist.matrix_R <- readRDS(paste0(base_dir, "/", radius,"/", radius, "/", n_run, "_run_Ndist_",n_reals, "_real_r_", radius, "_R_", sample_size ,".rds"))
    N.dist.matrix_C <- readRDS(paste0(base_dir, "/", radius,"/", radius, "/", n_run, "_run_Ndist_",n_reals, "_real_r_", radius, "_C_", sample_size ,".rds"))
    
    w <- weight(N.dist.matrix_C, N.dist.matrix_R)
    weights[n_run, ss] <- w
    
    a <- w/(w+1)
    a_from_w[n_run, ss] <-a
    
  }
} 

saveRDS(weights, file =paste0(n_run, "_weights_R_",n_reals, "_real_r_", radius, ".rds"))
saveRDS(a_from_w, file =paste0(n_run, "_a_weights_R_",n_reals, "_real_r_", radius,".rds"))
write.table(weights, file =paste0(n_run, "_weights_R_",n_reals, "_real_r_", radius,".txt"), row.names = FALSE, col.names = FALSE)
write.table(a_from_w, file =paste0(n_run, "_weights_R_",n_reals, "_real_r_", radius,".txt"), row.names = FALSE, col.names = FALSE)


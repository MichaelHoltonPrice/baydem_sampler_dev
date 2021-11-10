library(baydem)

rm(list=ls())

source("baydem_mcmc_candidates.R")

data_dir <- file.path(".","outputs")
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

set.seed(1005)
N <- 200
th_sim <-
  c(
    pi1 = 0.2,
    pi2 = 0.8,
    mu1 = 775,
    mu2 = 1000,
    sig1 = 35,
    sig2 = 45
  )

tau_min <- 600
tau_max <- 1300

calibration_curve <- "intcal20"
calib_df <- load_calib_curve(calibration_curve)
sim_spec <- list(model_spec=list(density_type = "trunc_gauss_mix",
                                 th=c(th_sim,tau_min,tau_max),
                                 error_spec=list(type="unif_fm",
                                                 min=.0021,max=.0028),
                                 is_AD=T),
                 N=N,
                 calib_curve=calibration_curve,
                 seed=1002)

sim <- simulate_rc_data(sim_spec)
rc_meas <- sim$data$rc_meas

density_model <- list(type="trunc_gauss_mix",
                      tau_min=600,
                      tau_max=1300,
                      K=2)
hp <-
  list(
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 3,
    # The gamma distribution rate parameter for sigma, yielding a mode of 300
    alpha_r = (3 - 1) / 300,
    # Spacing for the measurement matrix (years)
    dtau = 5
  )

prepped_inputs <-
  prep_for_trunc_gauss_mix_neg_log_lik(rc_meas,
                                       density_model,
                                       hp,
                                       calibration_curve,
                                       init_seed=20000)

t0 <- Sys.time()
# TODO: implement seed for reproducible chains
all_chain_data <- run_multiple_baydem_chains(prepped_inputs,
                                             num_chains=4,
                                             num_samp=400,
                                             verbose=TRUE)
t1 <- Sys.time()
print(difftime(t1, t0))

# coda appears to want samples to be the first index and variables the second
# index, so create a list of mcmc objects from the transpose of TH.
mcmc_list <- lapply(all_chain_data$chain_data_list,
                    function(chain_data){coda::mcmc(t(chain_data$chain$TH))})
print(coda::gelman.diag(mcmc_list))
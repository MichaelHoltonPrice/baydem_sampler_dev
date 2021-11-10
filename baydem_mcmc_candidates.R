sample_baydem_chain <- function(inputs, mcmc_control) {
  # Check that the prior and likelihood can be calculated at the initial point
  if (density_model$type == "trunc_gauss_mix") {
    log_prior0 <-  calc_trunc_gauss_mix_log_prior(inputs$th0,
                                                  inputs$density_model$tau_min,
                                                  inputs$density_model$tau_max,
                                                  inputs$hp)
    if (!is.finite(log_prior0)) {
      stop("log_prior cannot be evaluated at starting vector")
    }

    ll <- log_lik(inputs$th0, inputs)
    if (!is.finite(ll)) {
      stop("log likelihood cannot be evaluated at starting vector")
    }
  } else {
    stop(paste0("Unsupported density model type = ", density_model$type))
  }


  theta <- inputs$th0
  eta   <- log_prior0 + ll

  num_samp_tot <- mcmc_control$num_warmup + mcmc_control$num_samp
  accept_vect <- rep(NA,num_samp_tot)
  theta_mat <- matrix(NA,length(theta),num_samp_tot)
  for(n in 1:num_samp_tot) {
    print("--")
    print(n)
    delta_theta <- rnorm(length(theta))*mcmc_control$prop_scale
    theta_prop <- theta + delta_theta

    log_prior_prop <-
      calc_trunc_gauss_mix_log_prior(theta_prop,
                                     inputs$density_model$tau_min,
                                     inputs$density_model$tau_max,
                                     inputs$hp)
    if (is.finite(log_prior_prop)) {
      ll_prop <- log_lik(theta_prop, inputs)
      eta_prop <- log_prior_prop + ll_prop
      if (is.finite(eta_prop)) {
        alpha <- min(1,exp(eta_prop-eta))
        accept <- runif(1) < alpha
      } else {
        accept <- FALSE
      }
    } else {
      accept <- FALSE
    }

    accept_vect[n] <- accept
    print(accept)
    if(!accept) {
      theta_prop <- theta
      eta_prop   <- eta
    }

    # Get ready for next sample
    theta <- theta_prop
    theta_mat[,n] <- theta
    eta   <- eta_prop
    print(theta)
    print(eta)
  }
  return(theta_mat)
}

# TODO: support K=1
is_trunc_gauss_mix_valid <- function(th, tau_min, tau_max) {
  K <- (length(th) + 1) / 3

  # The mixture weights without the first one
  pi_reduced <- th[1:(K-1)]
  pi_full <- c(1-sum(pi_reduced), pi_reduced)

  # Probabilities must lie strictly on 0 to 1
  if (any(pi_full <= 0)) {
    return(FALSE)
  }

  if (any(pi_full >= 1)) {
    return(FALSE)
  }

  # mu must be ordered and all all values must lie between tau_min and tau_mix
  mu <- th[K:(2*K-1)]
  if (is.unsorted(mu)) {
    return(FALSE)
  }
  if (mu[1] < tau_min) {
    return(FALSE)
  }
  if (mu[length(mu)] > tau_max) {
    return(FALSE)
  }

  # s must be strictly positive
  s <- th[(2*K):(3*K-1)]
  if (any(s <= 0)) {
    return(FALSE)
  }
  return(TRUE)
}

# TODO: initialize tau_min and tau_max in the prep function
prep_for_trunc_gauss_mix_neg_log_lik <- function(rc_meas,
                                                 density_model,
                                                 hp,
                                                 calibration_curve,
                                                 init_seed=NA) {
  # TODO: Add save input_init_seed and init_seed to outputs
  input_init_seed <- init_seed
  if (is.na(init_seed)) {
    init_seed <- sample.int(1000000,1)
  }

  K <- density_model$K
  tau_min <- density_model$tau_min
  tau_max <- density_model$tau_max
  dtau    <- hp$dtau
  tau <- seq(tau_min, tau_max, by = dtau)
  th0 <- init_trunc_gauss_mix(K,
                              1,
                              tau_min,
                              tau_max,
                              input_seed=init_seed)

  N <- length(rc_meas$phi_m)

  calibrations <-
    Bchron::BchronCalibrate(ages=round(rc_meas$trc_m),
                            ageSds=round(rc_meas$sig_trc_m),
                            calCurves=rep(calibration_curve,N))

  # The minimum and maximum calendar date for each observation
  min_date_AD <-
    as.numeric(unlist(lapply(calibrations,
                             function(calib){min(1950 - calib$ageGrid)})))
  max_date_AD <-
    as.numeric(unlist(lapply(calibrations,
                             function(calib){max(1950 - calib$ageGrid)})))

  # Place min_AD_vect and max_AD_vect on the same sample spacing as tau
  if(dtau != 1) {
    for (n in 1:N) {
          min_date_AD[n] <- min_date_AD[n] - ( min_date_AD[n] %% dtau)
          max_date_AD[n] <- max_date_AD[n] + (-max_date_AD[n] %% dtau)
    }
  }

  # Ensure all dates are on tau_min to tau_max
  min_date_AD[min_date_AD < tau_min] <- tau_min
  max_date_AD[max_date_AD > tau_max] <- tau_max

  # At least for now, store M in a list for each observation
  M_list <- list()
  subset_matrix <- matrix(NA,N,2)
  calib_df <- load_calib_curve(calibration_curve)
  for (n in 1:N) {
    ind_min <- which(tau == min_date_AD[n])
    ind_max <- which(tau == max_date_AD[n])
    tau_n <- tau[ind_min:ind_max]
    # M has been multiplied by dtau
    M_n <- calc_meas_matrix(tau_n,
                            rc_meas$phi_m[n],
                            rc_meas$sig_m[n],
                            calib_df)
    M_list[[n]] <- M_n[1,]
    subset_matrix[n,] <- c(ind_min, ind_max)

  }
  return(list(th0=th0[2:length(th0)],
              tau=tau,
              density_model=density_model,
              hp=hp,
              calibration_curve=calibration_curve,
              calib_df=calib_df,
              M_list=M_list,
              subset_matrix=subset_matrix))
}

calc_trunc_gauss_mix_log_prior <- function(th, tau_min, tau_max, hp) {
  if (!is_trunc_gauss_mix_valid(th, tau_min, tau_max)) {
    return(Inf)
  }

  K <- (length(th) + 1) / 3

  # The mixture weights without the first one
  pi_reduced <- th[1:(K-1)]
  pi_full <- c(1-sum(pi_reduced), pi_reduced)
  pi_prior <- gtools::ddirichlet(pi_full, rep(hp$alpha_d, K))
  pi_log_prior <- log(pi_prior)

  mu_log_prior <- -K * log(tau_max - tau_min)

  s <- th[(2*K):(3*K-1)]
  s_log_prior <- sum(dgamma(s, shape=hp$alpha_s, rate=hp$alpha_r, log=TRUE))

  log_prior <- pi_log_prior + mu_log_prior + s_log_prior
  if (is.na(log_prior)) {
    stop("log_prior should not be NA")
  }

  if (!is.finite(log_prior)) {
    stop("log_prior should be finite")
  }
  return(pi_log_prior + mu_log_prior + s_log_prior)
}

calc_baydem_log_lik <- function(th, inputs, return_vect=FALSE) {
  K <- (length(th) + 1) / 3

  pi_reduced <- th[1:(K-1)]
  th <- c(1-sum(pi_reduced), th)
  N <- length(inputs$M_list)
  ll <- 0
  f <- calc_gauss_mix_pdf(th,
                          inputs$tau,
                          inputs$density_model$tau_min,
                          inputs$density_model$tau_max)
  f <- f / sum(f) / inputs$hp$dtau

  log_lik_vect <- rep(NA,N)
  for (n in 1:N) {
    ind_n <- inputs$subset_matrix[n,1]:inputs$subset_matrix[n,2]
    f_n <- f[ind_n]
    log_lik_vect[n] <- log(sum(inputs$M_list[[n]] * f_n))
  }

  if (return_vect) {
    return(log_lik_vect)
  } else {
    return(sum(log_lik_vect))
  }
}

run_multiple_baydem_chains <- function(prepped_inputs,
                                       num_chains=4,
                                       num_samp=1000,
                                       seed_vect=c(),
                                       verbose=FALSE) {
  input_seed_vect <- seed_vect

  if (length(seed_vect) == 0) {
    seed_vect <- sample.int(1000000, num_chains)
  }

  chain_data_list <- list()

  for (cc in 1:num_chains) {
    chain_data_list[[cc]] <- run_single_baydem_chain(prepped_inputs,
                                                     num_samp,
                                                     seed_vect[cc],
                                                     verbose)
  }

  # Number of samples per chain (loo notation)
  I <- nrow(chain_data_list[[1]]$log_lik_mat)
  # Number of observations
  N <- ncol(chain_data_list[[1]]$log_lik_mat)

  log_lik_array <- array(NA,c(I,num_chains,N))

  for (cc in 1:num_chains) {
    log_lik_array[,cc,] <- chain_data_list[[cc]]$log_lik_mat
  }
  loo_val <- loo::loo(log_lik_array)

  return(list(chain_data_list=chain_data_list,
              log_lik_array=log_lik_array,
              loo=loo_val))
}

run_single_baydem_chain <- function(prepped_inputs,
                                    num_samp=1000,
                                    seed=NA,
                                    verbose=FALSE) {
  input_seed <- seed
  if (is.na(seed)) {
    seed <- sample.int(1000000,1)
  }

  thinning <- 40
  tot_samp <- thinning * num_samp
  warmup <- baydem_warmup_tempering(prepped_inputs,
                                    verbose=verbose)
  if (verbose) {
    print("Sampling chain")
  }
  chain <- baydem_sample_mcmc(warmup$th,
                              tot_samp,
                              1,
                              prop_spec=list(prop_spec_type="normal",
                                             prop_scale=warmup$prop_scale),
                              prepped_inputs,
                              return_TH=TRUE,
                              verbose=FALSE)

  TH <- chain$TH[,seq(thinning, tot_samp, by=thinning)]

  # Calculate the log-likelihood for each sample and each observation. This is
  # used to call loo::loo for model selection.
  if (verbose) {
    print("Calculating log-likelihood matrix for chain")
  }
  N <- length(prepped_inputs$M_list)
  log_lik_mat <- matrix(NA, num_samp, N)
  for (s in 1:num_samp) {
    log_lik_mat[s,] <- calc_baydem_log_lik(TH[,s],
                                           prepped_inputs,
                                           return_vect=TRUE)
  }
  return(list(warmup=warmup, chain=chain, log_lik_mat=log_lik_mat))
}

# TODO: update documentation for baydem
#' @title
#' Sample a likelihood function at a given temperature for the posterior
#' distribution of a baydem density model.
#'
#' @description
#' Given a posterior density and temperature, do MCMC sampling of
#' the parameter vector. Let eta be the value of the log posterior.
#' The function that is sampled is exp(eta/temp), where temp is the
#' temperature; that is, the likelihood^(1/temp) is sampled. The proposal
#' is either drawn from the prior or made using an independent normal draw for
#' each variable in the parameter vector, th.
#'
#' @param th0 The initial value of the parameter vector
#' @param num_samp The number of samples to make
#' @param temp The temperature
#' @param prop_spec A specification of the proposal distribution
#' @param log_post0 The value of the log posterior for the starting parameter
#'   vector (default: NA, not used).
#' @param prop_spec A specification of the proposal distribution
#' @param return_TH Whether to return the sampled parameter vectors as a matrix
#'   (default: FALSE, not returned).
#' @param verbose Whether to print out progress information
#'
#' @return A list
#'
#' @export
baydem_sample_mcmc <- function(th0,
                               num_samp,
                               temp,
                               prop_spec,
                               prepped_inputs,
                               log_post0=NA,
                               return_TH=FALSE,
                               verbose=FALSE) {
  if (is.na(log_post0)) {
    log_post0 <- calc_baydem_log_post(th0, prepped_inputs)
  }

  if(!is.finite(log_post0)) {
    stop('The log posterior is not finite at th0')
  }

  th       <- th0
  log_post <- log_post0

  if(return_TH) {
    TH <- matrix(NA, length(th), num_samp)
  }

  if (verbose) {
    print("Sampling chain")
  }
  for(n in 1:num_samp) {
    if (prop_spec$prop_spec_type == "normal") {
      th_prop <- th + rnorm(length(th))*prop_spec$prop_scale}
    else if (prop_spec$prop_spec_type == "prior") {
      th_prop <- sample_baydem_post(prepped_inputs)
    } else {
      stop(paste0("Unrecognized prop spec type = ", prop_spec$prop_spec_type))

    }
    log_post_prop <- calc_baydem_log_post(th_prop, prepped_inputs)

    if(!is.finite(log_post_prop)) {
      accept <- F
    } else {
      alpha <- min(1,exp((log_post_prop-log_post)/temp))
      accept <- runif(1) < alpha
    }

    if(!accept) {
      th_prop <- th
      log_post_prop   <- log_post
    }

    # Get ready for next sample
    th        <- th_prop
    log_post  <- log_post_prop
    if (return_TH) {
      TH[,n] <- th
    }
    if (verbose) {
      print(log_post_prop)
    }
  }

  if (!return_TH) {
    return(list(th=th, log_post=log_post))
  } else {
    return(list(th=th, log_post=log_post, TH=TH))
  }
}

calc_baydem_log_post <- function(th, prepped_inputs) {

  if (prepped_inputs$density_model$type != "trunc_gauss_mix") {
    stop(paste0("Unrecognized density model type = ",
                prepped_inputs$density_model$type))
  }

  # Check that the parameter vector is valid. If not, return Inf
  tau_min <- prepped_inputs$density_model$tau_min
  tau_max <- prepped_inputs$density_model$tau_max
   if (!is_trunc_gauss_mix_valid(th, tau_min , tau_max)) {
    return(Inf)
  }

  # Calculate the log prior
  K <- (length(th) + 1) / 3

  pi_reduced <- th[1:(K-1)]
  pi_full <- c(1-sum(pi_reduced), pi_reduced)
  pi_prior <- gtools::ddirichlet(pi_full, rep(hp$alpha_d, K))
  pi_log_prior <- log(pi_prior)

  mu_log_prior <- -K * log(tau_max - tau_min)

  s <- th[(2*K):(3*K-1)]
  s_log_prior <- sum(dgamma(s, shape=hp$alpha_s, rate=hp$alpha_r, log=TRUE))
  log_prior <- pi_log_prior + mu_log_prior + s_log_prior

  if (is.na(log_prior)) {
    return(Inf)
  }

  if (!is.finite(log_prior)) {
    return(Inf)
  }

  # Calculat the log-likelihood
  log_lik <- calc_baydem_log_lik(th, prepped_inputs)
  log_post <- log_prior + log_lik
  return(log_post)
}

sample_baydem_post <- function(prepped_inputs) {
  # Number of mixtures
  K <- prepped_inputs$density_model$K

  # Sample mixture weights
  pi_full <- gtools::rdirichlet(1, rep(prepped_inputs$hp$alpha_d, K))
  pi_reduced <- pi_full[2:length(pi_full)]

  # Sample mu (then order)
  mu <- runif(K,
              prepped_inputs$density_model$tau_min,
              prepped_inputs$density_model$tau_max)
  mu <- sort(mu)

  # Sample s
  s <- rgamma(K,
              shape=prepped_inputs$hp$alpha_s,
              rate=prepped_inputs$hp$alpha_r)
  th <- c(pi_reduced, mu, s)
  return(th)
}

temp_vect_builder <- function(ratio_in_exp, num_temp) {
  temp_vect <- rep(NA, num_temp)
  temp_vect[1] <- 1
  for (n in 2:num_temp) {
    temp_vect[n] <- 1/ (1/temp_vect[n-1] - ratio_in_exp)
  }
  return(temp_vect)
}

baydem_warmup_tempering <- function(prepped_inputs,
                                    verbose=FALSE) {

  samp_per_cyc <- 5
  num_cyc <- 10000

  # Create the temperature vector, which is ten temperatures from 1 to
  # 100 evenly ordered in log space (constant ratios of adjacent temperatures).
  temp_vect <- exp(seq(log(1), log(100), len=20))
  #temp_vect <- temp_vect_builder(1/20, 20)

  # Create prop_spec_list, which is the list of proposal distribution
  # specifications. It is the prior density for the hottest temperature and
  # scaled normal draws for all other temperatures.

  # Number of mixtures
  K <- prepped_inputs$density_model$K

  # lo and hi scales for the mixture weights
  pi_scale_lo <- 0.01 / K / 2
  pi_scale_hi <- 0.10 / K / 2

  # lo and hi scales for the means
  tau_range <- prepped_inputs$density_model$tau_max -
    prepped_inputs$density_model$tau_min
  mu_scale_lo <- tau_range / 1000 / 2
  mu_scale_hi <- tau_range / 10 / 2

  # Lo and hi scales for the standard deviations
  s_mode <-  (prepped_inputs$hp$alpha_s - 1) / prepped_inputs$hp$alpha_r
  s_scale_lo <- s_mode / 1000 / 2
  s_scale_hi <- s_mode / 10 / 2

  # Add scales for all but the last temperature
  num_temp <- length(temp_vect)
  num_to_add <- num_temp - 1
  prop_spec_list <- list()

  for (temp_num in 1:num_to_add) {
    ratio <- (temp_num - 1) / (num_to_add - 1)
    pi_scale <- pi_scale_lo + (pi_scale_hi - pi_scale_lo) * ratio
    mu_scale <- mu_scale_lo + (mu_scale_hi - mu_scale_lo) * ratio
     s_scale <-  s_scale_lo + ( s_scale_hi -  s_scale_lo) * ratio
    prop_scale <- c(rep(pi_scale, K-1),
                    rep(mu_scale, K),
                    rep( s_scale, K))
    prop_spec_list[[temp_num]] <- list(prop_spec_type = "normal",
                                       prop_scale=prop_scale)
  }
  prop_spec_list[[num_temp]] <- list(prop_spec_type = "prior")

  # Iterate over number of cycles
  chain_list <- list()

  # Store the parameter values for the the coldest temperature. The standard
  # deviations of TH1 are used to set the proposal scale for the ultimate
  # sampling for this chain.
  TH1 <- matrix(NA, 3*K-1, num_cyc*samp_per_cyc)
  for(cc in 1:num_cyc) {
    if(verbose) {
      print("**********")
      print("Warmup tempering")
      print(paste0("Cycle ", cc, " of ", num_cyc))
    }
    # Extend chains for this cycle
    for(temp_num in 1:num_temp) {
      if (cc == 1) {
        th0 <- sample_baydem_post(prepped_inputs)
        log_post <- NA
      } else {
        th0      <- chain_list[[temp_num]]$th
        log_post <- chain_list[[temp_num]]$log_post
      }
      return_TH <- (temp_num == 1) || (temp_num == 2)
      chain_list[[temp_num]] <- baydem_sample_mcmc(th0,
                                                   samp_per_cyc,
                                                   temp_vect[temp_num],
                                                   prop_spec_list[[temp_num]],
                                                   prepped_inputs,
                                                   log_post=log_post,
                                                   return_TH=return_TH)
    }

    # Swaps for the coldest temperature must be handled below
    ind_lo <- (cc-1)*samp_per_cyc + 1
    ind_hi <- ind_lo + samp_per_cyc - 1
    TH1[,ind_lo:ind_hi] <- chain_list[[1]]$TH

    # Randomly choose an adjacent pair of temperatures to attempt a swap
    temp_num <- sample(1:(num_temp-1),1)
    log_post_cold <- chain_list[[temp_num    ]]$log_post
    log_post_hot  <- chain_list[[temp_num + 1]]$log_post
    temp_cold <- temp_vect[temp_num]
    temp_hot  <- temp_vect[temp_num + 1]
    swap_prob <- exp(-(log_post_cold-log_post_hot)*(1/temp_cold-1/temp_hot))

    log_post_vect <- unlist(lapply(chain_list,
                                   function(chain) {chain$log_post}))
    if (verbose) {
      print(log_post_vect)
    }

    if(swap_prob >= 1) {
      accept_swap <- T
    } else {
      accept_swap <- runif(1) <= swap_prob
    }

    if (verbose) {
      if (accept_swap) {
        print(paste0("Swap ", temp_num))
      } else {
        print("No swap")
      }
    }

    if(accept_swap) {
      cold_chain <- chain_list[[temp_num]]
      chain_list[[temp_num]]     <- chain_list[[temp_num + 1]]
      chain_list[[temp_num + 1]] <- cold_chain
      if (temp_num == 1) {
        TH1[,ind_hi] <- chain_list[[2]]$TH[,samp_per_cyc]
      }
    }

    if (verbose) {
      print("Current parameter vector")
      print(chain_list[[1]]$th)
    }
  }

  # Extract the final parameter vector for a temperature of 1, the coldest
  # temperature
  th <- chain_list[[1]]$th

  # Use the last 2000 samples of TH1 to set the proposal distribution scale
  # for subsequent sampling. The rescaling by 2.4 / sqrt(num_param) comes from
  # Harrio, Saksman, and Tamminen 2001 -- An Adaptive Metropolis Algorithm
  prop_scale <- 2.4 * apply(TH1[,8001:10000],1,sd) / sqrt(length(th))
  return(list(th=th, prop_scale=prop_scale))
}
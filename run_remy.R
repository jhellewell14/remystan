library(cmdstanr)
library(bayesplot)
library(ggoutbreak)
available_mcmc(pattern = "_nuts_")


## Set a load of ground truth parameters for simulation from Rob's code
sim_dt <- data.table::data.table(variable = c("R0", "inc_mean", "gt_mean", "delay_mean", "i0", "gt_sd", "inc_sd", "delay_sd"),
                                 true = c(2.2, 3.5, 3.5, 6, 6, 1, 1, 2))

## Runs Rob's simulation code to give linelist data
source("05_simulation.R")

## Convert linelist into onset case counts
cases <- data.table::as.data.table(sim_cases)
min_date <- min(c(cases$OnsetDate, cases$Date_of_Interview))
cases$obs_date <- max(cases$Date_of_Interview) + 1
onsets <- cases[, .(onset = as.integer(OnsetDate - min_date),
                    case_notification = as.integer(Date_of_Interview - min_date))]
onsets <- onsets[, .(onsets = .N), onset]
maxons <- max(onsets$onset)
onsets <- onsets[data.table::CJ(onset = 0:maxons), on = "onset"][, .(onset, count = data.table::nafill(onsets, fill = 0))]
plot(onsets)


## Set up stan data object
gtgam <- ggoutbreak:::.reparam_gamma(mean = sim_dt[variable == "gt_mean", true], sd = sim_dt[variable == "gt_sd", true])
lead_days <- 0
data <- list(ngen = 9,
             max_time = nrow(onsets) + lead_days,
             gt_shape = gtgam$shape,
             gt_rate = gtgam$rate,
             # i0 = sim_dt[variable == "i0", true],
             ip_shape = gtgam$shape,
             del_mean = log(sim_dt[variable == "delay_mean", true]),
             del_sd = log(sim_dt[variable == "delay_sd", true]),
             sigma = 0.25,
             R_init = sim_dt[variable == "R0", true],
             obs = c(rep(0, lead_days), onsets$count),
             lead_days = lead_days,
             init = c(0, 0.5),
             init_mat = matrix(c(1, -0.99, -0.99, 1), nrow = 2)
             )

## Simulation model for playing with
mod <- cmdstan_model("remy_multi.stan")
res <- mod$sample(data, fixed_param = TRUE, iter_sampling = 1, chains = 1)
as.vector(res$draws("R"))
plot(res$draws("reports"))
plot(as.vector(res$draws("reports")))


## Model for fitting
fmod <- cmdstan_model("remy_fit.stan")
init_fn <- function(){list("sigma" = 0.25, R = rlnorm(data$ngen - 1, 0, 0.1))}
fres <- fmod$sample(data, parallel_chains = 4, init = init_fn, iter_warmup = 1000, iter_sampling = 1000)
lp_cp <- log_posterior(fres)
np_cp <- nuts_params(fres)

mcmc_intervals(fres$draws("Rg")) +
  geom_point(data = data.frame(x = paste0("Rg[",1:data$ngen, "]"), value = sim_dt[variable == "R0", true]),
             aes(x = value, y = x), col = "firebrick", size = 2, inherit.aes = FALSE) +
  labs(title = "Reproduction number by generation")

mcmc_hist(fres$draws(c("sigma", "initial")))
fres$summary()
mcmc_parcoord(fres$draws("Rg"), np = np_cp) +
  geom_hline(yintercept = sim_dt[variable == "R0", true], col = "firebrick", lty = 2, linewidth = 1) +
  labs(title = "Time-varying reproduction number")
mcmc_parcoord(fres$draws("onsets")) +
  geom_point(data = data.frame(x = paste0("onsets[",1:data$max_time, "]"), y = data$obs), aes(x = x, y = y),
             inherit.aes = FALSE, col = "firebrick") +
  labs(title = "Unobserved onsets")

## We miss data at the start because the data starts at t = 0 first onset and we use t = 0 first exposure
## I'd like to formulate a modified time distribution for generation 0 that assigns first 1 to 2 cases on day 0
## and the rest in the coming days
## Like an incubation period distribution conditional on the first onset having happened, then where do the rest fall?
mcmc_parcoord(fres$draws("reports")) +
  geom_point(data = data.frame(x = paste0("reports[",1:data$max_time, "]"), y = data$obs), aes(x = x, y = y),
             inherit.aes = FALSE, col = "firebrick") +
  labs(title = "Onsets by date of report")

mcmc_parcoord(fres$draws("dbm")) +
  # geom_point(data = data.frame(x = paste0("reports[",1:data$max_time, "]"), y = data$obs), aes(x = x, y = y),
  #            inherit.aes = FALSE, col = "firebrick") +
  labs(title = "Onsets by date of report")


mcmc_trace(fres$draws(), pars = "sigma", np = np_cp)


plot(res$draws("dbm"))
plot(as.vector(res$draws("dbm")))
plot(as.vector(res$draws("onsets")))
plot(as.vector(res$draws("reports")))

plot(as.vector(res$draws("del_dist")))

gtdist <- matrix(as.vector(res$draws("gt_dist", format = "matrix")[1,]), ncol = data$ngen, nrow = data$max_time)
dim(gtdist)
# Contribution to time point 10 from each generation
plot((gtdist[10,] / sum(gtdist[10,])))
# Generation time distribution for 1st generation
plot(gtdist[, 1], type = 'l')
for(i in 2:ncol(gtdist)){
  lines(gtdist[, i])
}
lines(gtdist2[, 1], col = "red")
lines(gtdist[, 2], col = "red")
lines(gtdist[, 3], col = "blue")
lines(gtdist[, 4], col = "green")
lines(gtdist[, 5], col = "pink")

plot(pgamma(1:26, shape = 12.25, rate = 3.5, lower.tail = TRUE, log.p = FALSE) -
       pgamma(0:25, shape = 12.25, rate = 3.5, lower.tail = TRUE, log.p = FALSE), ylim = c(0, 1))

points((1 - (1 - pgamma(1:26, 12.25, 3.5))^5) - (1 - (1 - pgamma(0:25, 12.25, 3.5))^5), col = 'red')

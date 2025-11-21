options(
  repos = c(
    "ai4ci" = 'https://ai4ci.r-universe.dev/',
    "terminological" = 'https://terminological.r-universe.dev/',
    CRAN = 'https://cloud.r-project.org'
  )
)

library(tidyverse)

# Setup siumulation ----

sim_params = list(
  # A short generation time
  mean_gt = sim_dt[variable == "gt_mean", true],
  sd_gt = sim_dt[variable == "gt_sd", true],
  R0 = sim_dt[variable == "R0", true],
  I0 = sim_dt[variable == "i0", true],
  # Add a longish and very variable delay to symptoms
  p_symptomatic = 1,
  mean_incub = sim_dt[variable == "inc_mean", true],
  sd_incub = sim_dt[variable == "inc_sd", true],
  # and a slightly shorter delay to observation:
  # Only symptomatic cases are observed
  p_detected_given_symptoms = 1,
  mean_obs = sim_dt[variable == "delay_mean", true],
  sd_obs = sim_dt[variable == "delay_sd", true],
  # Observation cutoff:
  T = 30
)

# Run simulation ----

sim_ip = ggoutbreak::make_gamma_ip(
  median_of_mean = sim_params$mean_gt,
  median_of_sd = sim_params$sd_gt
)

sim_params$r0 = ggoutbreak::inv_wallinga_lipsitch(sim_params$R0, sim_ip)

truth = ggoutbreak::sim_branching_process(
  fn_Rt = ~ sim_params$R0,
  fn_ip = ~sim_ip,
  fn_imports = \(t) ifelse(t == 1, sim_params$I0, 0),
  max_time = 40
)

delayed = truth %>%
  ggoutbreak::sim_delay(
    p_fn = ~ sim_params$p_symptomatic,
    delay_fn = ~ rlnorm(
      .x,
      log(sim_params$mean_incub),
      sim_params$sd_incub
    ),
    input = "time",
    output = "symptom"
  ) %>%
  ggoutbreak::sim_delay(
    p_fn = \(t, symptom) {
      ifelse(symptom, sim_params$p_detected_given_symptoms, 0)
    },
    delay_fn = ~ ggoutbreak::rgamma2(
      .x,
      sim_params$mean_obs,
      sim_params$sd_obs
    ),
    input = "symptom_time",
    output = "observation"
  )

observed = delayed %>% dplyr::filter(observation_time < sim_params$T)
traced_contacts = observed %>%
  dplyr::semi_join(observed, by = c("infector" = "id"))

# Visualise simulation ----

if (FALSE) {
  ggoutbreak::plot_cases(
    observed %>%
      rename(infection_time = time, time = symptom_time) %>%
      mutate(time = floor(time))
  ) +
    geom_step(
      data = truth %>% dplyr::filter(time < sim_params$T),
      aes(x = as.Date(time)),
      stat = "bin",
      binwidth = 1,
      inherit.aes = FALSE
    ) +
    ggoutbreak::scale_y_log1p()

  serial_pairs = observed %>%
    inner_join(
      traced_contacts,
      by = c("id" = "infector"),
      suffix = c(".1", ".2")
    ) %>%
    mutate(
      serial_interval.2 = floor(symptom_time.2) - floor(symptom_time.1)
    ) %>%
    glimpse()

  ggplot(serial_pairs) +
    geom_histogram(aes(x = serial_interval.2), binwidth = 1)

  # This is shorter than the generation time by a lot.
  mean(serial_pairs$serial_interval.2)
  # The earlier generations have generation interval as expected
  truth %>%
    filter(generation < 9) %>%
    pull(generation_interval) %>%
    mean(na.rm = TRUE)
  # but this is lower in the whole set due to truncation of longer SI in
  mean(truth$generation_interval, na.rm = TRUE)
}

# Reformat to match pegasus ----

# Match cases and contacts
sim_cases = observed %>%
  transmute(
    person_id = id,
    tracer_id = ifelse(infector %in% traced_contacts$infector, infector, -1),
    OnsetDate = as.Date(floor(symptom_time)),
    Date_of_Interview = as.Date(floor(observation_time))
  )

# Sim contacts is only positive contacts because the BPM does not have non-infected
# all contacts have a tracer_id
sim_contacts = traced_contacts %>%
  transmute(
    person_id = id,
    tracer_id = infector,
    OnsetDate = as.Date(floor(symptom_time)),
    Date_of_Interview = as.Date(floor(observation_time))
  )

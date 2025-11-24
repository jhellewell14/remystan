
data {
  int<lower=0> ngen; // number of generations
  int<lower=1> max_time; // length of time series
  real<lower=0> gt_shape;
  real<lower=0> gt_rate;
  // int<lower=0> i0;
  real<lower=0> ip_shape;
  // real<lower=0> ip_rate; // going to assume ip_rate = gt_rate
  real<lower=0> del_mean;
  real<lower=0> del_sd;
  real<lower=0> R_init;
  real<lower=0> sigma;
}

generated quantities {
  // Generate discrete delay distribution
  vector[max_time] del_dist;
  for(i in 1:max_time){
    del_dist[max_time - i + 1] = lognormal_cdf(i| del_mean, del_sd);
  }

  // Generate discrete generation distribution in real time for each generation
  matrix[max_time, ngen] gt_dist;
  for(i in 1:ngen){
    for(j in 1:max_time){
      // For gt_dist col 1, generation 0, I want to use a different distribution to account for the fact that our
      // time series starts at first observed case
      gt_dist[j, i] = gamma_cdf(j + (ip_shape / gt_rate)| gt_shape * i, gt_rate) -
      gamma_cdf(j + (ip_shape / gt_rate) - 1| gt_shape * i, gt_rate);
    }
  }
  // Will be fitting distribution shape and rate parameters, they will need decent priors
  // I think take gt shape and ip shape from multivariate normal to ensure correlation

  // Discrete Geometric random walk for generation sizes
  vector[ngen] R;
  R[1] = R_init;
  for(i in 2:ngen){
    R[i] = normal_rng(R[i - 1], sigma);
  }

  vector[ngen] dbm;
  dbm[1] = R[1];
  for(i in 2:ngen){
    dbm[i] = dbm[i - 1] * R[i - 1];
  }

  // Calculates daily infections by summing contributions from each generation
  vector[max_time] onsets = gt_dist * dbm;
  array[max_time] int reports = poisson_rng(onsets .* del_dist);

  // Next step is
  // observed_infections[t] ~ Poisson(infections[t]);
  // This is a Poisson approximation to what should be Poisson Binomial
  // Could also try Binomial approximation

  // If stan did integers I'd do this

  // vector[max_time + 1] infections = rep_vector(0, max_time + 1);
  // for(k in 1:ngen){
  //   infections += to_vector(multinomial_rng(gt_dist[, k], gen_size[k]));
  // }

}

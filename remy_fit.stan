
data {
  int<lower=0> ngen; // number of generations
  int<lower=1> max_time; // length of time series
  real<lower=0> gt_shape;
  real<lower=0> gt_rate;
  real<lower=0> ip_shape;
  // real<lower=0> ip_rate; // going to assume ip_rate = gt_rate
  real<lower=0> del_mean;
  real<lower=0> del_sd;
  array[max_time] int obs;
  int<lower=0> lead_days;
}

parameters {
  real<lower=0> sigma;
  vector<lower=0>[ngen] Rg;
  // real<lower=0> gt_shape;
  // real<lower=0> gt_rate;
  // real<lower=0> ip_shape;
}

transformed parameters {

  vector[max_time] del_dist;
  for(i in 1:max_time){
    del_dist[max_time - i + 1] = lognormal_cdf(i| del_mean, del_sd);
  }

  // Generate discrete generation distribution in real time for each generation
  matrix[max_time, ngen] gt_dist;
  for(i in 1:ngen){
    for(j in 1:max_time){
      gt_dist[j, i] = gamma_cdf(j | gt_shape * i + ip_shape, gt_rate) - gamma_cdf(j - 1| gt_shape * i + ip_shape, gt_rate);
    }
  }
  // Will be fitting distribution shape and rate parameters, they will need decent priors
  // I think take gt shape and ip shape from multivariate normal to ensure correlation

  // Rg (R for each generation) links previous generation size to next generation size
  vector[ngen] dbm;
  dbm[1] = Rg[1];
  for(i in 2:ngen){
    dbm[i] = dbm[i - 1] * Rg[i - 1];
  }

  // Calculates daily infections by summing contributions from each generation
  // Sum over all generations of Prob(gen(i) case on day t * size of gen(i))
  vector[max_time] onsets = gt_dist * dbm;
  vector[max_time] reports = onsets .* del_dist;
}

model {
  //
  obs ~ poisson(reports);
  // rw_noise ~ lognormal(0, sigma);
  sigma ~ gamma(1, 5);
  Rg[1] ~ normal(1, sigma)T[0,];
  for(i in 2:ngen){
    Rg[i] ~ normal(Rg[i - 1], sigma)T[0, ];
  }
  // gt_shape ~ std_normal();
  // gt_rate ~ std_normal();
  // ip_shape ~ std_normal();
}

generated quantities {
  vector[max_time] Rt;
  for(t in 1:max_time){
    Rt[t] = (gt_dist[t, 2:ngen] / sum(gt_dist[t, 2:ngen])) * Rg[2:ngen];
  }
}


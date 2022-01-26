functions {
  real n_beta_given_m_cdf(int n, real m, real alpha_0, real alpha_1, 
                          real sigma, int b, int a, real beta) {
  return Phi(       // Phi is CDF of normal
   (
    ((b - a) * m - beta) - 
     (alpha_0 * (b - a - n) + alpha_1/2 * ((b - n) * (b + n + 1) - a * (a + 1)))
   ) / (sigma * sqrt(b + a - n))
   );
  }
}

data {
  int<lower=0> N;           // number of experiments n = 1, ..., N
  int<lower=0> S;           // number of sites
  int<lower=0> Y;           // number of years
  int site[N];              // site of experiment n
  int year[N];              // year of experiment n
  vector[N] M;              // mean temperature in experiment n
  int n_beta[N];            // leaf out day in year n
  int a;                    // beginning of average temperature window 
  int b;                    // end of average temperature of window
  real alpha_0;             // baseline temperature on day 0
}

parameters {
  real alpha_1[Y];          // daily increase in temperature over the 'spring'
  real<lower = 0> beta[Y];  // thermal sum
  real<lower = 0> sigma[S]; // noise around the temperature measure
  real gamma_0_alpha;       // trend line over time: intercept for temperature
  real gamma_1_alpha;       // trend line over time: daily increase in temperature
  real gamma_0_beta;        // trend line over time: intercept for thermal sum
  real gamma_1_beta;        // trend line over time: change in beta (slope) over time
  real tau_alpha;           // sigma on the trend line for alpha
  real tau_beta;            // sigma on the trend line for beta
}

model {
  for(n in 1:N) {
    target += 
//log likelihood of n_beta | M    
    log(n_beta_given_m_cdf(n_beta[n], M[n], // WTF? Why 1 unit apart?
                           alpha_0, alpha_1[year[n]], 
                           sigma[site[n]], b, a, beta[year[n]]) -
                  n_beta_given_m_cdf(n_beta[n] - 1, M[n], 
                           alpha_0, alpha_1[year[n]], 
                           sigma[site[n]], b, a, beta[year[n]])) +      
//log likelihood of M
               normal_lpdf(M[n] | 
                           alpha_0 + alpha_1[year[n]]/2 * (b + a + 1), 
                           sigma[site[n]] / sqrt(b-a));
  }
//allow alpha_1 and beta to change linearly over time
for(y in 1:Y)   {
  alpha_1[y] ~ normal(gamma_0_alpha + gamma_1_alpha * y, tau_alpha);
  beta[y]    ~ normal(gamma_0_beta  + gamma_1_beta  * y, tau_beta);
  }
}

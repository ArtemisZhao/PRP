// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <math.h>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;

void MH_step(double &bbar, double &k, double hbo, double so,  vector<double> &k_vec, double & prop_k, gsl_rng *gr){
  
  double prop_sd = sqrt(prop_k*prop_k*hbo*hbo+so*so);
  
  // propose a new set of values
  double bbar_new = gsl_ran_gaussian(gr, prop_sd)+hbo;
  double k_new = k_vec[rand()%k_vec.size()];
  
  
  //ratio of proposed vs. current in proposal distribution
  double pr = gsl_ran_gaussian_pdf(bbar_new - hbo, prop_sd)/gsl_ran_gaussian_pdf(bbar - hbo, prop_sd);
  
  // ratio of proposed vs. current in stationary distribution
  double lr = gsl_ran_gaussian_pdf(hbo - bbar_new , sqrt(k_new*k_new*bbar_new*bbar_new + so*so))/gsl_ran_gaussian_pdf(hbo - bbar, sqrt(k*k*bbar*bbar + so*so));
  
  // Acceptance ratio
  double mhr = lr/pr;
  
  
  if(gsl_rng_uniform(gr) <= mhr){
    bbar = bbar_new;
    k = k_new;
  }
}


double run_MCMC(double hbo, double so, double hbr, double sr, vector<double> &k_vec, int nreps, double burnin_rate) {
  
  vector<double> bbar_sample;
  vector<double> k_sample;
  
  // init positions
  double bbar = hbo;
  double k = k_vec[rand()%k_vec.size()];
  
  
  // init gsl_rng
  const gsl_rng_type * T;
  gsl_rng * gr;
  
  gsl_rng_env_setup();
  T = gsl_rng_default;
  gr = gsl_rng_alloc (T);
  
  
  double prop_k = *max_element(k_vec.begin(), k_vec.end());
  
  
  for (int i = 0; i < nreps; i++){
    MH_step(bbar, k, hbo, so, k_vec, prop_k, gr);
    bbar_sample.push_back(bbar);
    k_sample.push_back(k);
  }
  
  /*
   for( int i=0; i< nreps;i++){
   printf("%f %f\n",bbar_sample[i], k_sample[i]);
   }
   exit(0);
   */
  
  int start = int(nreps*burnin_rate);
  int count = 0;
  double pv1_sum = 0;
  for (int i = start; i< nreps; i++){
    double pv1 = gsl_cdf_gaussian_P(hbr/hbo - bbar_sample[i]/hbo, 
                                   sqrt((pow(k_sample[i],2)*pow(bbar_sample[i],2)+sr*sr)/pow(hbo,2)));
  
     pv1_sum += pv1;
     count++;
  }
  
  gsl_rng_free(gr);
  double pval = pv1_sum;
  return pval/count;
}

// [[Rcpp::export]]
double PRP_MCMC(vector<double> beta, vector<double> se, vector<double> k_vec){
  double hbo = beta[0];
  double hbr = beta[1];
  double so = se[0];
  double sr = se[1];
  double pval = run_MCMC(hbo, so, hbr, sr, k_vec, 100000, 0.01);
  return pval;
}
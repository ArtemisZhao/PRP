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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#define MAX_ITERATIONS 100
using namespace std;

struct prior_PRP_rst {
  
  double pvalue;
  double left_bound;
  double right_bound;
};

struct sample_info_params {
  double sd;
  vector<double> bbar_vec;
  vector<double> k_vec;
};

double left_tail_prob(double x, void *params){
  
  struct sample_info_params *p;
  p = (struct sample_info_params *)params;
  double rst = 0;
  for (int i=0; i < p->bbar_vec.size(); i++){
    rst += gsl_cdf_gaussian_P(x - p->bbar_vec[i], sqrt(pow(p->k_vec[i],2)*pow(p->bbar_vec[i],2)+pow(p->sd,2)));
  }
  rst =  rst/p->bbar_vec.size();
  return rst - 0.025;
}

double right_tail_prob(double x, void *params){
  struct sample_info_params *p;
  p = (struct sample_info_params *)params;
  double rst = 0;
  for (int i=0; i < p->bbar_vec.size(); i++){
    rst += gsl_cdf_gaussian_Q(x - p->bbar_vec[i], sqrt(pow(p->k_vec[i],2)*pow(p->bbar_vec[i],2)+pow(p->sd,2)));
  }
  rst = rst/p->bbar_vec.size();
  return rst - 0.025;
}

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

struct prior_PRP_rst run_MCMC(double hbo, double so, double hbr, double sr, vector<double> &k_vec, int nreps, double burnin_rate) {
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
  
  int start = int(nreps*burnin_rate);
  int count = 0;
  double pv1_sum = 0;
  double pv2_sum = 0;
  double max_bbar = -9e9;
  double min_bbar = 9e9;
  vector<double> bbv;
  vector<double> kv;
  for (int i = start; i< nreps; i++){
    
    double pv1 = gsl_cdf_gaussian_P(hbr - bbar_sample[i], sqrt(pow(k_sample[i],2)*pow(bbar_sample[i],2)+sr*sr));
    pv1_sum += pv1;
    pv2_sum += 1-pv1;
    count++;
    
    bbv.push_back(bbar_sample[i]);
    kv.push_back(k_sample[i]);
    
    if(bbar_sample[i]>max_bbar){
      max_bbar = bbar_sample[i];
    }
    
    if(bbar_sample[i]< min_bbar){
      min_bbar = bbar_sample[i];
    }
  }
  
  gsl_rng_free(gr);
  
  double pval = pv1_sum;
  if(pval> pv2_sum)
    pval= pv2_sum;
  
  pval = 2*pval/count;
  
  // compute predictive intervals:
    
  struct sample_info_params params;
  params.sd = sr;
  params.bbar_vec = bbv;
  params.k_vec = kv;
  
  
  int status;
  const gsl_root_fsolver_type *solver_type;
  gsl_root_fsolver *solver;
  gsl_function F1, F2;
  double x_lo, x_hi;
  
  x_lo = min_bbar - 3*sqrt(prop_k*prop_k*min_bbar*min_bbar+sr*sr);
  x_hi = max_bbar;
  
  
  F1.function = &left_tail_prob;
  F1.params = &params;

  solver_type = gsl_root_fsolver_bisection;
  solver = gsl_root_fsolver_alloc(solver_type);
  gsl_root_fsolver_set(solver, &F1, x_lo, x_hi);

  status = GSL_CONTINUE;
  double sol;
  for (int i = 1; i <= MAX_ITERATIONS && status == GSL_CONTINUE; ++i) {
    status = gsl_root_fsolver_iterate(solver);
    /* get the solver's current best solution and bounds */
        sol = gsl_root_fsolver_root(solver);
        x_lo = gsl_root_fsolver_x_lower(solver);
        x_hi = gsl_root_fsolver_x_upper(solver);
        status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
    }
    double left_bound = NAN;
    if(status == GSL_SUCCESS){
        left_bound = sol;
    }

    x_lo = min_bbar;
    x_hi = max_bbar + 3*sqrt(prop_k*prop_k*min_bbar*min_bbar+sr*sr);

    F2.function = &right_tail_prob;
    F2.params = &params;
    gsl_root_fsolver_set(solver, &F2, x_lo, x_hi);
    
    status = GSL_CONTINUE;
    for (int i = 1; i <= MAX_ITERATIONS && status == GSL_CONTINUE; ++i) {
        status = gsl_root_fsolver_iterate(solver);
        /* get the solver's current best solution and bounds */
      sol = gsl_root_fsolver_root(solver);
      x_lo = gsl_root_fsolver_x_lower(solver);
      x_hi = gsl_root_fsolver_x_upper(solver);
      status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
  }
  
  double right_bound = NAN;
  if(status == GSL_SUCCESS){
    right_bound = sol;
  }
  
  gsl_root_fsolver_free(solver);
  
  struct prior_PRP_rst rst;
  rst.pvalue = pval;
  rst.left_bound = left_bound;
  rst.right_bound = right_bound;
  return rst;
}

// int main(int argc, char **argv){
//   vector<double>  k_vec{ 0 , 0.155368927298941 , 0.166005089694258 , 0.173337899205751 , 0.178978966221498 , 0.18388938695547 , 0.188067112764603 , 0.191791542513613 , 0.195271332191419 , 0.198575148690133 , 0.201534613519245 , 0.204481460061963 , 0.207059615021334 , 0.209648317091446 , 0.212212828958339 , 0.21443824547403 , 0.216731336214093 , 0.219075253455047 , 0.221090629773248 , 0.223007024495773 , 0.225129102734274 , 0.227105121873094 , 0.229073533229895 , 0.230927660776478 , 0.232684191110707 , 0.234490009720275 , 0.236500593128432 , 0.23807026226874 , 0.239798064432984 , 0.24137077435528 , 0.242982232444001 , 0.244687134655899 , 0.246358284045983 , 0.247862152521707 , 0.249551446935383 , 0.251016347829005 , 0.252623405238645 , 0.25419738927027 , 0.25548147869798 , 0.25715799531597 , 0.258514034625798 , 0.260067172562903 , 0.26149385602969 , 0.262818252318785 , 0.264347547831445 , 0.265740978439278 , 0.267262201879966 , 0.268573378227942, 0.269859998040611, 0.271219052445701 , 0.272612228503996 };
// 
//   // vector<double> k_vec {0};
//   // vector<double> k_vec{0.272612228503996};
//   
//   ifstream infile(argv[1]);
//   string line;
//  
//   while(getline(infile,line)){
//     string nulls;
//     string id;
//     double hbo;
//     double so;
//     double hbr;
//     double sr;
//   
//     istringstream ins(line);
//     if(ins>>nulls>>id>>hbo>>so>>hbr>>sr){
//       struct prior_PRP_rst rst = run_MCMC(hbo, so, hbr, sr, k_vec, 100000, 0.01);
//       printf("%s  %f %f %f %f     %f  (%f , %f)\n", id.c_str(), hbo, so, hbr, sr, rst.pvalue,  rst.left_bound, rst.right_bound); 
//     }
//   }
// 
// }

// [[Rcpp::export]]
Rcpp::List PRP_MCMC_with_PI(vector<double> beta, vector<double> se, vector<double> k_vec){
  double hbo = beta[0];
  double hbr = beta[1];
  double so = se[0];
  double sr = se[1];
  struct prior_PRP_rst rst = run_MCMC(hbo, so, hbr, sr, k_vec, 100000, 0.01);
  //vector<double> pval = rst.pvalue;
  //vector<double> PI = {rst.left_bound,rst.right_bound};
  return Rcpp::List::create(Rcpp::Named("pval") = rst.pvalue,
                             Rcpp::Named("left_bound") = rst.left_bound,
                             Rcpp::Named("right_bound") = rst.right_bound);
}
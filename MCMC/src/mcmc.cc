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
    double pv2_sum = 0;
    for (int i = start; i< nreps; i++){
        double pv1 = gsl_cdf_gaussian_P(hbr - bbar_sample[i], sqrt(pow(k_sample[i],2)*pow(bbar_sample[i],2)+sr*sr));
        pv1_sum += pv1;
        pv2_sum += 1-pv1;
        count++;
    }

    gsl_rng_free(gr);
    double pval = pv1_sum;
    if(pval> pv2_sum)
        pval= pv2_sum;

    return 2*pval/count;
}



int main(int argc, char **argv){

    vector<double>  k_vec{ 0 , 0.155368927298941 , 0.166005089694258 , 0.173337899205751 , 0.178978966221498 , 0.18388938695547 , 0.188067112764603 , 0.191791542513613 , 0.195271332191419 , 0.198575148690133 , 0.201534613519245 , 0.204481460061963 , 0.207059615021334 , 0.209648317091446 , 0.212212828958339 , 0.21443824547403 , 0.216731336214093 , 0.219075253455047 , 0.221090629773248 , 0.223007024495773 , 0.225129102734274 , 0.227105121873094 , 0.229073533229895 , 0.230927660776478 , 0.232684191110707 , 0.234490009720275 , 0.236500593128432 , 0.23807026226874 , 0.239798064432984 , 0.24137077435528 , 0.242982232444001 , 0.244687134655899 , 0.246358284045983 , 0.247862152521707 , 0.249551446935383 , 0.251016347829005 , 0.252623405238645 , 0.25419738927027 , 0.25548147869798 , 0.25715799531597 , 0.258514034625798 , 0.260067172562903 , 0.26149385602969 , 0.262818252318785 , 0.264347547831445 , 0.265740978439278 , 0.267262201879966 , 0.268573378227942, 0.269859998040611, 0.271219052445701 , 0.272612228503996 };

    // vector<double> k_vec {0};

    // vector<double> k_vec{0.272612228503996};

    ifstream infile(argv[1]);
    string line;


    while(getline(infile,line)){
        string nulls;
        string id;
        double hbo;
        double so;
        double hbr;
        double sr;

        istringstream ins(line);
        if(ins>>nulls>>id>>hbo>>so>>hbr>>sr){
            double pval = run_MCMC(hbo, so, hbr, sr, k_vec, 100000, 0.01);
            printf("%s  %f %f %f %f     %f\n", id.c_str(), hbo, so, hbr, sr, pval); 
        }
    }




}




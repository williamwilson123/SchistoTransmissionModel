// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <Rcpp.h>
#include <Rmath.h>
#include <RcppNumerical.h>
#include <cmath>

using namespace Numer;
using namespace std;
using namespace Rcpp;

// CLASS for integrating the mating probability
class MatingProbIntegrand: public Func
{
private:
  double alpha;
  double k;
public:
  MatingProbIntegrand(double alpha_, double k_) : alpha(alpha_), k(k_) {}
  
  double operator()(const double& x) const
  {
    return (1-cos(x))/ (pow( (1+alpha*cos(x)), (1+k) ));
  }
};


//------------------------------------------------------------------------------
// DEFINE CONSTANTS & STORAGE ELEMENTS
//------------------------------------------------------------------------------

// define global hard-wired parameters
const int na = 51, amax = 51;
// const int na = 50, amax = 50;
const int nw = 10;
const double da = na / static_cast<double>(amax);
const double aged = 1.0/da;
const double shape = 1.11, scale = 20.69;
const double beta1_male = 78.00111563, beta1_female = 65.4666457, beta2_male = 0.08858346, beta2_female = 0.0673592;
const double sr = 0.427;
const double prop_female = sr, prop_male = 1 - sr;
const double muW = 0.2;
const double muS = 12;
const double rho = 0.5;
const double lower = 0.0, upper=2.0*M_PI; 
double err_est;
int err_code;
double res;  // for mating probability integration
const double tolerance = 1e-5; //tolerance for floating point comparisons
int coverage_index = 0;       //to locate time-specific coverage data during MDA

// define parameters that are assigned values at run time
int nsteps;
double R0, R0_weight, NhNs, kW, decay_immunity, protection_immunity, epg_constant, ige_constant;

// treatment parameters
int toggle_tx, toggle_covdat, n_tx, min_tx_age, max_tx_age;
double start_tx, freq_tx, coverage, efficacy;
double tx_coverage;  //for storing time-specific coverage in time for-loop
double cov_weight;
NumericVector age_specific_cov(na); //age-specific coverage

// Indices for demographic groups (WHO Guideline on Control & Elimination of Schistosomiasis, 2022)
IntegerVector sac_indices   = seq(5*(na/amax), 15*(na/amax));     //5-15 year olds
IntegerVector adult_indices = seq(16*(na/amax), na-1);  //16+ year olds
  
//define temporary vectors for mating probability
NumericVector mating_probability_male(na), mating_probability_female(na);

// tolerance measure to check derivatives by end of simulation are near 0
NumericVector row_tol_female(na), row_tol_male(na);
NumericVector col_tol_female(nw+1), col_tol_male(nw+1);
NumericVector tol(4);
double tol_out = 0;

// define global storage vectors and scales
NumericVector age(na);                                                                      // age mid points
NumericVector demographic_weights(na);                                                     // demographic weights
NumericVector normalised_female_exposure(na), normalised_male_exposure(na);                // normalized female/male exposure
NumericVector worm_burden_age_female(na), worm_burden_age_male(na);                       // age-dependent female/male worm burdens
NumericVector worm_burden_age_female_equib(na), worm_burden_age_female_post(na), worm_burden_age_male_equib(na), worm_burden_age_male_post(na);
NumericVector cumulative_worms_age_male(na), cumulative_worms_age_female(na);
NumericVector kpost_female(na), kpost_male(na);
NumericVector kW_age_female(na), kW_age_male(na);
double kW_female, kW_male;
double worm_burden_male, worm_burden_female;
double worm_burden;
double cumulative_worms_male, cumulative_worms_female;
double cumulative_worms;
NumericVector psi_female_age(na), psi_male_age(na);                                         // age-dependent contamination 
double net_contamination;                                                                   // net contaimation
NumericVector epg_age_male(na), epg_age_female(na), epg_age_male_SAC(na), epg_age_female_SAC(na);
NumericVector prevalence_age_female(na), prevalence_age_male(na), prevalence_age_female_SAC(na), prevalence_age_male_SAC(na);
NumericVector ige_age_male(na), ige_age_female(na);
double epg_male, epg_female, epg, epg_male_SAC, epg_female_SAC, epg_SAC;
double prevalence_male, prevalence_female, prevalence, prevalence_male_SAC, prevalence_female_SAC, prevalence_SAC;
double ige_male, ige_female, ige;
double dead_worms;
double wormk;

// age-dependent force of infection for males and females
NumericVector infect_snails_female(na), infect_snails_male(na);

// define state variable matrices (where nw+1 counts dead worms & nw+2 counts cumulative experience of live worms)
NumericMatrix male_states(na, nw+2), updated_male_states(na,nw+2), female_states(na, nw+2), updated_female_states(na,nw+2);
NumericMatrix tmp_male_states(na, nw+2), tmp_female_states(na,nw+2);

// define matrix for counting worms killed following treatment
NumericMatrix male_killed_worms(na, nw), female_killed_worms(na, nw);

// define state variable for cumulative number of dead worms in each class
NumericVector dead_worms_male(na), dead_worms_female(na);

//define state variable acquired immunity
NumericVector male_acquired_immunity(na), female_acquired_immunity(na);

// define matrices of current derivatives
NumericMatrix dWdt_male(na, nw+2), dWdt_female(na,nw+2); //+1 for dead W, +2 for cumulative W

// define matrices for RK4 integration
NumericMatrix k1_male(na,nw+2), k1_female(na,nw+2), k2_male(na,nw+2), 
k2_female(na,nw+2), k3_male(na,nw+2), k3_female(na,nw+2), k4_male(na,nw+2), k4_female(na,nw+2);

// create vectors to store no. worms killed via treatment
NumericVector summed_male_killed_worms(na), summed_female_killed_worms(na); 

// create vector to store treatment times & MDA coverage values
NumericVector treatment_times(n_tx), coverage_data(n_tx); 


//------------------------------------------------------------------------------
// DEFINE FUNCTIONS 
//------------------------------------------------------------------------------

// define useful R function for later use --------------------------------------
Function pnbinom_func("pnbinom");  

// Function to calculate age midpoints -----------------------------------------
void AgeMidpoints() {
  for (int i = 0; i < na; i++) {
    age[i] = i*aged + 0.5*aged;
  }
}

// Function to calculate demographic weights -----------------------------------
void DemographicWeights() {
  for(int i = 0; i < na; i++) {
    double ageplus = age[i] + 0.5*aged;
    double ageminus = age[i] - 0.5*aged;
    demographic_weights[i] = (R::pweibull(ageplus, shape, scale,1,0) -
      R::pweibull(ageminus, shape, scale,1,0)) / R::pweibull(amax,  shape,  scale,1,0);
  }
}

// Function to calculate normalized exposure -----------------------------------
void NormalisedExposure() {
  
  NumericVector female_exposure(na), male_exposure(na);
  
  for (int i = 0; i < na; i++) {
    male_exposure[i]   = beta1_male*age[i]*exp(-beta2_male*age[i]);
    female_exposure[i] = beta1_female*age[i]*exp(-beta2_female*age[i]);
  }

  double acc_male_exposure = 0, acc_female_exposure = 0;

  for (int i = 0; i < na; i++) {
    acc_female_exposure += demographic_weights[i]*female_exposure[i];
    acc_male_exposure   += demographic_weights[i]*male_exposure[i];
  }

  double norm_constant = 1/(prop_female*acc_female_exposure+prop_male*acc_male_exposure);

  for (int i = 0; i < na; i++) {
    normalised_female_exposure[i] = female_exposure[i]*norm_constant;
    normalised_male_exposure[i]   = male_exposure[i]*norm_constant;
  }
  
}

// Function to initialize storage vectors --------------------------------------
void InitialiseEverything() {
  for (int i = 0; i < na; i++) {
    mating_probability_female[i] = 1;
    mating_probability_male[i] = 1;
    worm_burden_age_male[i] = 0;
    worm_burden_age_female[i] = 0;
    male_acquired_immunity[i] = 0;
    female_acquired_immunity[i] = 0;
    dead_worms_female[i] = 0;
    dead_worms_male[i] = 0;
    epg_age_female[i] = 0;
    epg_age_male[i] = 0;
    epg_age_female_SAC[i] = 0;
    epg_age_male_SAC[i] = 0;
    prevalence_age_female[i] = 0;
    prevalence_age_male[i] = 0;
    prevalence_age_female_SAC[i] = 0;
    prevalence_age_male_SAC[i] = 0;
    ige_age_female[i] = 0;
    ige_age_male[i] = 0;
    cumulative_worms_age_female[i] = 0;
    cumulative_worms_age_male[i] = 0;
    age_specific_cov[i] = 0;
  }
  worm_burden_female = 0;
  worm_burden_male = 0;
  worm_burden = 0;
  epg_female = 0;
  epg_male = 0;
  epg = 0;
  epg_female_SAC = 0;
  epg_male_SAC = 0;
  epg_SAC = 0;
  ige_female = 0;
  ige_male = 0;
  ige = 0;
  dead_worms = 0;
  cumulative_worms_male = 0;
  cumulative_worms_female = 0;
  cumulative_worms = 0;
}

// Function to generate initial state variables based on input parameters ------
void SmartInitialiseStates() {

  for(int i = 0; i < na; i++) {
    for (int j = 0; j < nw; j++) {
      male_states(i,j) = ( (R0/rho)*muS*normalised_male_exposure[i]) /(muW*nw*NhNs*pow(R0,R0_weight));
      female_states(i,j) = ( (R0/rho)*muS*normalised_female_exposure[i])/(muW*nw*NhNs*pow(R0,R0_weight));
    }
    // cumulative dead worms
    male_states(i,nw)   = 0;
    female_states(i,nw) = 0;
    // cumulative live worms
    male_states(i,nw+1)   = 0;
    female_states(i,nw+1) = 0;
  }
}

// Function to calculate age & sex-specific worm burdens -----------------------
void WormBurden() {

  // Incident worm burden
  // by age
  for(int i = 0; i < na; i++) {
    double accmale_age = 0;
    double accfemale_age = 0;
    for (int j = 0; j < nw; j++) {
      accmale_age += male_states(i,j);
      accfemale_age += female_states(i,j);
    }
    worm_burden_age_male[i] = accmale_age;
    worm_burden_age_female[i] = accfemale_age;
  }

  // by sex
  double acc_worm_burden_male = 0;
  double acc_worm_burden_female = 0;
  for(int i = 0; i < na; i++) {
    acc_worm_burden_male += worm_burden_age_male[i]*demographic_weights[i];
    acc_worm_burden_female += worm_burden_age_female[i]*demographic_weights[i];
  }

  worm_burden_male = acc_worm_burden_male;
  worm_burden_female = acc_worm_burden_female;

  //total population
  worm_burden = prop_female*worm_burden_female + prop_male*worm_burden_male;

}


// Function to calculate age & sex-specific cumulative worm burdens ------------
void CumulativeWormBurden() {
  
  // by age
  double accmale_age   = 0;
  double accfemale_age = 0;
  
  for(int i = 0; i < na; i++) {
    accmale_age   += male_states(i,nw+1);
    accfemale_age += female_states(i,nw+1);
    cumulative_worms_age_male[i]   = accmale_age;
    cumulative_worms_age_female[i] = accfemale_age;
  }
  
  // by sex
  double acc_worm_burden_male   = 0;
  double acc_worm_burden_female = 0;
  
  for(int i = 0; i < na; i++) {
    acc_worm_burden_male   += cumulative_worms_age_male[i]*demographic_weights[i];
    acc_worm_burden_female += cumulative_worms_age_female[i]*demographic_weights[i];
  }
  cumulative_worms_male   = acc_worm_burden_male;
  cumulative_worms_female = acc_worm_burden_female;
  
  // total population
  cumulative_worms = prop_female * cumulative_worms_female + prop_male * cumulative_worms_male;
  
}


// Function to calculate infection intensity (epg) -----------------------------
void EPG(){
  // by age
  for (int i=0; i<na; i++) { // total population
    epg_age_female[i] = worm_burden_age_female[i]*mating_probability_female[i]*epg_constant;
    epg_age_male[i] = worm_burden_age_male[i]*mating_probability_male[i]*epg_constant;
  }
  for (int i=5; i<15; i++) { // school-aged children
    epg_age_female_SAC[i] = epg_age_female[i];
    epg_age_male_SAC[i] = epg_age_male[i];
  }
  
  // by sex
  double acc_epg_male = 0;
  double acc_epg_male_SAC = 0; 
  double acc_epg_female = 0; 
  double acc_epg_female_SAC = 0;
  
  for (int i = 0; i < na; i++) {
    acc_epg_male += epg_age_male[i]*demographic_weights[i];
    acc_epg_female += epg_age_female[i]*demographic_weights[i];
    acc_epg_male_SAC += epg_age_female_SAC[i]*demographic_weights[i];
    acc_epg_female_SAC += epg_age_male_SAC[i]*demographic_weights[i];
  }
  epg_male = acc_epg_male;
  epg_female = acc_epg_female;
  epg_male_SAC = acc_epg_male_SAC;
  epg_female_SAC = acc_epg_female_SAC;
  
  //weighted mean for given timepoint
  epg = prop_female*epg_female + prop_male*epg_male; //total population
  epg_SAC = prop_female*epg_female_SAC + prop_male*epg_male_SAC; //SAC
}


// Function to calculate prevalence from epg -----------------------------------
void Prevalence() {
  // by age
  // calculate probability of 0 eggs
    prevalence_age_female = pnbinom_func(0, Named("size")=kW_age_female, _["mu"]=epg_age_female);
    prevalence_age_male   = pnbinom_func(0, Named("size")=kW_age_male,   _["mu"]=epg_age_male);
  // calculate 1 minus that probability to get prevalence
  prevalence_age_female = 1 - prevalence_age_female;
  prevalence_age_male   = 1 - prevalence_age_male;
  // calculate prevalence in school-aged children
  for (int i=5; i<15; i++) { 
    prevalence_age_female_SAC[i] = prevalence_age_female[i];
    prevalence_age_male_SAC[i] = prevalence_age_male[i];
  }
  
  //by sex
  double acc_prevalence_male = 0; 
  double acc_prevalence_male_SAC = 0; 
  double acc_prevalence_female = 0;
  double acc_prevalence_female_SAC = 0;
  
  for(int i = 0; i < na; i++) {
    acc_prevalence_male += prevalence_age_male[i]*demographic_weights[i];
    acc_prevalence_female += prevalence_age_female[i]*demographic_weights[i];
    acc_prevalence_male_SAC += prevalence_age_female_SAC[i]*demographic_weights[i];
    acc_prevalence_female_SAC += prevalence_age_male_SAC[i]*demographic_weights[i];
  }
  prevalence_male = acc_prevalence_male;
  prevalence_female = acc_prevalence_female;
  prevalence_male_SAC = acc_prevalence_male_SAC;
  prevalence_female_SAC = acc_prevalence_female_SAC;
  
  //weighted mean for given timepoint
  prevalence = prop_female*prevalence_female + prop_male*prevalence_male; //total population
  prevalence_SAC = prop_female*prevalence_female_SAC + prop_male*prevalence_male_SAC; //SAC
}


// Function to calculate IgE from acquired immunity variable -------------------
void IgE(){
  // by age
  for (int i = 0; i < na; i++) {
    ige_age_female[i] = female_acquired_immunity[i]*ige_constant;
    ige_age_male[i] = male_acquired_immunity[i]*ige_constant;
  }
  //by sex
  double acc_ige_male = 0; 
  double acc_ige_female = 0;
  for(int i = 0; i < na; i++) {
    acc_ige_male += ige_age_male[i]*demographic_weights[i];
    acc_ige_female += ige_age_female[i]*demographic_weights[i];
  }
  ige_male = acc_ige_male;
  ige_female = acc_ige_female;
  //total population
  ige = prop_female*ige_female + prop_male*ige_male;
}


// Function to calculate cumulative dead worms ---------------------------------
void DeadWorms() {
  double acc_dead_worms_male   = 0;
  double acc_dead_worms_female = 0;
  for(int i = 0; i < na; i++) {
    acc_dead_worms_male   += male_states(i,nw)   * demographic_weights[i];
    acc_dead_worms_female += female_states(i,nw) * demographic_weights[i];
    dead_worms_male[i]     = male_states(i,nw);
    dead_worms_female[i]   = female_states(i,nw);
  }
  //total population
  dead_worms = prop_male*acc_dead_worms_male + prop_female*acc_dead_worms_female;
}


// Function to calculate mating probabilities ----------------------------------
void MonogamousMatingFunc() {
  for(int i = 0; i < na; i++) {
    if (worm_burden_age_female[i] > 0) {
      double Wtmp_female = worm_burden_age_female[i];
      double ktmp_female = kW_age_female[i];
      double alpha = Wtmp_female/(ktmp_female+Wtmp_female);
      MatingProbIntegrand f(alpha, ktmp_female);
      res = integrate(f, lower, upper, err_est, err_code);
      mating_probability_female[i] = 1 - ( pow((1-alpha),(1+ktmp_female)) / (2*M_PI) )*res;
    } else {
      mating_probability_female[i] = 1;
    }
    if (worm_burden_age_male[i] > 0) {
      double Wtmp_male = worm_burden_age_male[i];
      double ktmp_male = kW_age_male[i];
      double alpha = Wtmp_male/(ktmp_male+Wtmp_male);
      MatingProbIntegrand f(alpha, ktmp_male);
      res = integrate(f, lower, upper, err_est, err_code);
      mating_probability_male[i] = 1 - ( pow((1-alpha),(1+ktmp_male)) / (2*M_PI) )*res;
    } else {
      mating_probability_male[i] = 1;
    }
  }
}


// Function to calculate derivatives -------------------------------------------
List Derivs() {

  double acc_net_contamination = 0;

  for (int i = 0; i < na; i++) {
    psi_female_age[i] = normalised_female_exposure[i]*demographic_weights[i]*worm_burden_age_female[i]*mating_probability_female[i];
    psi_male_age[i] = normalised_male_exposure[i]*demographic_weights[i]*worm_burden_age_male[i]*mating_probability_male[i];
    acc_net_contamination += prop_female*psi_female_age[i] + prop_male*psi_male_age[i];
  }

  net_contamination = acc_net_contamination;

  for (int i = 0; i < na; i++) {
    infect_snails_female[i] = (R0*muS*normalised_female_exposure[i] * net_contamination * exp(-protection_immunity * female_acquired_immunity[i])) / (rho * pow(R0,R0_weight) * NhNs * net_contamination + muS/muW);
    infect_snails_male[i]   = (R0*muS*normalised_male_exposure[i]   * net_contamination * exp(-protection_immunity * male_acquired_immunity[i]  )) / (rho * pow(R0,R0_weight) * NhNs * net_contamination + muS/muW);
  }

  double acc_infect_snails_female;
  double acc_infect_snails_male;
  
  for(int i = 0; i < na; i++){
    acc_infect_snails_female += infect_snails_female[i];
    acc_infect_snails_male   += infect_snails_male[i];
  }

  // differential equations
  for (int j = 0; j < (nw+2); j++) { // nw live worm groups, 1 dead worm group, 1 cumulative worm group
    for (int i = 0; i < na; i++) {   //host age groups

      if (j == 0) {  //first worm group
        if (i == 0) {
          // first age group
          dWdt_male(i,j)   = infect_snails_male[i]   - (nw*muW+da)*tmp_male_states(i,j);
          dWdt_female(i,j) = infect_snails_female[i] - (nw*muW+da)*tmp_female_states(i,j);

        } else {
          // subsequent age groups
          dWdt_male(i,j)   = infect_snails_male[i]   - (nw*muW+da)*tmp_male_states(i,j)   + da*tmp_male_states(i-1,j);
          dWdt_female(i,j) = infect_snails_female[i] - (nw*muW+da)*tmp_female_states(i,j) + da*tmp_female_states(i-1,j);
        }

      } else if (j > 0 & j < nw) {  //subsequent worm groups
        if (i == 0) {
          // first age group
          dWdt_male(i,j)   = nw*muW*tmp_male_states(i,j-1)   - (nw*muW+da)*tmp_male_states(i,j);
          dWdt_female(i,j) = nw*muW*tmp_female_states(i,j-1) - (nw*muW+da)*tmp_female_states(i,j);
        } else {
          // subsequent age groups
          dWdt_male(i,j)   = nw*muW*tmp_male_states(i,j-1)   - (nw*muW+da)*tmp_male_states(i,j)   + da*tmp_male_states(i-1,j);
          dWdt_female(i,j) = nw*muW*tmp_female_states(i,j-1) - (nw*muW+da)*tmp_female_states(i,j) + da*tmp_female_states(i-1,j);
        }

      } else if (j == nw) {  //dead worm group
        if (i == 0) {
          // first age group
          dWdt_male(i,j)   = nw*muW*tmp_male_states(i,j-1) - da*tmp_male_states(i,j);
          dWdt_female(i,j) = nw*muW*tmp_female_states(i,j-1) - da*tmp_female_states(i,j);
        } else {
          // subsequent age groups
          dWdt_male(i,j)   = nw*muW*tmp_male_states(i,j-1) + da*tmp_male_states(i-1,j) - da*tmp_male_states(i,j);
          dWdt_female(i,j) = nw*muW*tmp_female_states(i,j-1) + da*tmp_female_states(i-1,j) - da*tmp_female_states(i,j);
        }
      } else if (j == nw+1) {  //cumulative worm group
        dWdt_male(i,j)   = infect_snails_male[i]   - da*tmp_male_states(i,j);
        dWdt_female(i,j) = infect_snails_female[i] - da*tmp_female_states(i,j);
      }

    }
  }

  // return derivatives
  List derivatives;
  derivatives["dWdt_male"]   = dWdt_male;
  derivatives["dWdt_female"] = dWdt_female;
  return derivatives;

}


// Function to implement RK4 integration ---------------------------------------
List RK4(NumericMatrix male_states, NumericMatrix female_states, double stepsize){
  
  // set current states in a temporary matrix (use Rcpp::clone() to ensure changes in tmp states don't affect states)
  tmp_male_states   = Rcpp::clone(male_states);
  tmp_female_states = Rcpp::clone(female_states);

  // calculate derivatives for temporary states
  List derivs_result1 = Derivs();

  // store derivatives
  k1_male   = Rcpp::clone(as<NumericMatrix>(derivs_result1["dWdt_male"]));
  k1_female = Rcpp::clone(as<NumericMatrix>(derivs_result1["dWdt_female"]));

  for (int i = 0; i < na; i++) {
    // for (int j = 0; j < (nw+1); j++) {
    for (int j = 0; j < (nw+2); j++) {
      // move temporary matrices on by stepsize/2
      tmp_male_states(i,j) = male_states(i,j) + (static_cast<double>(stepsize)/2)*k1_male(i,j);
      tmp_female_states(i,j) = female_states(i,j) + (static_cast<double>(stepsize)/2)*k1_female(i,j);
    }
  }
  // calculate derivatives for temporary states and store
  List derivs_result2 = Derivs();

  // store derivatives
  k2_male   = Rcpp::clone(as<NumericMatrix>(derivs_result2["dWdt_male"]));
  k2_female = Rcpp::clone(as<NumericMatrix>(derivs_result2["dWdt_female"]));


  for (int i = 0; i < na; i++) {
    for (int j = 0; j < (nw+2); j++) {
      // move temporary matrices on by stepsize/2
      tmp_male_states(i,j) = male_states(i,j) + (static_cast<double>(stepsize)/2)*k2_male(i,j);
      tmp_female_states(i,j) = female_states(i,j) + (static_cast<double>(stepsize)/2)*k2_female(i,j);
    }
  }
  // calculate derivatives for temporary states and store
  List derivs_result3 = Derivs();

  // store derivatives
  k3_male   = Rcpp::clone(as<NumericMatrix>(derivs_result3["dWdt_male"]));
  k3_female = Rcpp::clone(as<NumericMatrix>(derivs_result3["dWdt_female"]));

  for (int i = 0; i < na; i++) {
    for (int j = 0; j < (nw+2); j++) {
      // move temporary matrices on by stepsize
      tmp_male_states(i,j) = male_states(i,j) + stepsize*k3_male(i,j);
      tmp_female_states(i,j) = female_states(i,j) + stepsize*k3_female(i,j);
    }
  }
  // calculate derivatives for temporary states and store
  List derivs_result4 = Derivs();
  k4_male   = Rcpp::clone(as<NumericMatrix>(derivs_result4["dWdt_male"]));
  k4_female = Rcpp::clone(as<NumericMatrix>(derivs_result4["dWdt_female"]));

  //update the state matrices
  for (int i = 0; i < na; i++) {
    for (int j = 0; j < (nw+2); j++) {
      updated_male_states(i,j) = male_states(i,j) + (static_cast<double>(stepsize)/6)*( k1_male(i,j) + 2*k2_male(i,j) + 2*k3_male(i,j) + k4_male(i,j) );
      updated_female_states(i,j) = female_states(i,j) + (static_cast<double>(stepsize)/6)*( k1_female(i,j) + 2*k2_female(i,j) + 2*k3_female(i,j) + k4_female(i,j) );
    }
  }

  // Create a List to store both updated states and return it
  List result;
  result["updated_male_states"]   = updated_male_states;
  result["updated_female_states"] = updated_female_states;
  return result;
}


// Function to update the state matrices ---------------------------------------
void UpdateStates() {
  male_states   = Rcpp::clone(updated_male_states);
  female_states = Rcpp::clone(updated_female_states);
}


// Function to update all model outputs using state matrices -------------------
void UpdateEverything() {
  // calculate current worm burden
  WormBurden();
  // calculate current mating probabilities
  MonogamousMatingFunc();
  // calculate cumulative dead worms
  DeadWorms();
  // calculate current epgs
  EPG();
  // calculate current prevalence
  Prevalence();
  // calculate current IgE
  IgE();
}


// Function to calculate group-specific treatment coverage ---------------------
NumericVector WeightCoverage(double cov_weight, double coverage, IntegerVector sac_indices, IntegerVector adult_indices, const int na) {

  // Demographic group coverage
  double presac_coverage = 0;
  double sac_coverage    = coverage;
  double adult_coverage  = sac_coverage * cov_weight;
  
  // Create age-specific coverage vector
  NumericVector final_coverage(na);
  
  // pre-SAC
  for (int i = 0; i < min(sac_indices); i++) {
    final_coverage[i] = 0.0;
  }
  // SAC
  for (int i = 0; i < sac_indices.size(); i++) {
    final_coverage[sac_indices[i]] = sac_coverage;
  }
  // adults
  for (int i = 0; i < adult_indices.size(); i++) {
    final_coverage[adult_indices[i]] = adult_coverage;
  }
  
  return final_coverage;
}


// Function to simulate community treatment  -----------------------------------
void TreatmentEvent() {
  
  // Treat between minimum & maximum treatment ages
  for (int i = min_tx_age; i < max_tx_age; i++) {

    //initialise storage for worms killed by treatment
    double tmp_summed_male_killed_worms   = 0;
    double tmp_summed_female_killed_worms = 0;

    for (int j = 0; j < nw; j++) {
      // calculate worms killed by treatment
      male_killed_worms(i,j)   = male_states(i,j)  * age_specific_cov[i] * efficacy;
      female_killed_worms(i,j) = female_states(i,j)* age_specific_cov[i] * efficacy;
      // update worm burden following treatment
      updated_male_states(i,j)   = male_states(i,j)   - male_killed_worms(i,j);
      updated_female_states(i,j) = female_states(i,j) - female_killed_worms(i,j);
      // sum the worms killed by treatment
      tmp_summed_male_killed_worms   += male_killed_worms(i,j);
      tmp_summed_female_killed_worms += female_killed_worms(i,j);
    }
    // then store killed worms for each host age group
    summed_male_killed_worms[i]   = tmp_summed_male_killed_worms;
    summed_female_killed_worms[i] = tmp_summed_female_killed_worms;
    
    // and add treatment-killed worms to dead worm states
    updated_male_states(i,nw)   = male_states(i,nw)   + summed_male_killed_worms[i];
    updated_female_states(i,nw) = female_states(i,nw) + summed_female_killed_worms[i];
  }
  
}


// Function to update k DURING treatment event ---------------------------------
void DynamickTreatment() {
  
  for (int i = 0; i < na; i++) {
    // get worm burden after treatment
    worm_burden_age_female_post[i] = worm_burden_age_female[i];
    worm_burden_age_male_post[i]   = worm_burden_age_male[i];
    // calculate age-specific k post-treatment
    kpost_female[i] = (kW*worm_burden_age_female_post[i])/((1+kW)*worm_burden_age_female_equib[i]-kW*worm_burden_age_female_post[i]);
    kpost_male[i]   = (kW*worm_burden_age_male_post[i])  /((1+kW)*worm_burden_age_male_equib[i]  -kW*worm_burden_age_male_post[i]);
    // update age-specific k
    kW_age_female[i] = pow(worm_burden_age_female[i], 2) * ( pow( (worm_burden_age_female_equib[i] - worm_burden_age_female_post[i]), 2)) /(((pow(worm_burden_age_female_equib[i], 2))/kW) * pow((worm_burden_age_female[i] - worm_burden_age_female_post[i]), 2) + (pow(worm_burden_age_female_post[i], 2) / kpost_female[i]) * (pow((worm_burden_age_female[i] - worm_burden_age_female_equib[i]), 2)));
    kW_age_male[i]   = pow(worm_burden_age_male[i], 2)   * ( pow( (worm_burden_age_male_equib[i]   - worm_burden_age_male_post[i]), 2))   /(((pow(worm_burden_age_male_equib[i], 2))/kW)   * pow((worm_burden_age_male[i]   - worm_burden_age_male_post[i]), 2)   + (pow(worm_burden_age_male_post[i], 2)   / kpost_male[i])   * (pow((worm_burden_age_male[i]   - worm_burden_age_male_equib[i]), 2)));
  }
  
  // calculate weighted mean k for given timepoint
  double acc_kW_female = 0, acc_kW_male = 0;
  
  for (int i = 0; i < na; i++) {
    acc_kW_female += kW_age_female[i]*demographic_weights[i];
    acc_kW_male   += kW_age_male[i]  *demographic_weights[i];
  }
  
  kW_female = acc_kW_female;
  kW_male   = acc_kW_male;
  wormk     = prop_female*kW_female + prop_male*kW_male;
  
}


// Function to update k AFTER treatment event ----------------------------------
void DynamickPostTreatment() {
  
  for (int i = 0; i < na; i++) {
    kW_age_female[i] = ( kW * worm_burden_age_female[i] ) / ( (1+kW) * worm_burden_age_female_equib[i] - kW * worm_burden_age_female[i]);
    kW_age_male[i]   = ( kW * worm_burden_age_male[i]   ) / ( (1+kW) * worm_burden_age_male_equib[i]   - kW * worm_burden_age_male[i]  );
  }
    
  // calculate weighted mean k at given timepoint
  double acc_kW_female = 0, acc_kW_male = 0;
  
  for (int i = 0; i < na; i++) {
    acc_kW_female += kW_age_female[i]*demographic_weights[i];
    acc_kW_male   += kW_age_male[i]  *demographic_weights[i];
  }
  
  kW_female = acc_kW_female;
  kW_male   = acc_kW_male;
  wormk     = prop_female*kW_female + prop_male*kW_male;
  
}


//------------------------------------------------------------------------------
// DEFINE TRANSMISSION MODEL FUNCTION 
//------------------------------------------------------------------------------

// [[Rcpp::export]]
List RunModel(NumericVector theta, int runtime, double stepsize, int alltimes, 
              NumericVector tx_pars, NumericVector tx_times, NumericVector coverage_data) {
  
  // Define parameters to be varied
  R0                  = std::exp(theta["R0"]);
  R0_weight           = std::exp(theta["R0_weight"]);
  NhNs                = std::exp(theta["NhNs"]);
  kW                  = std::exp(theta["kW"]);
  decay_immunity      = std::exp(theta["decay_immunity"]);
  protection_immunity = std::exp(theta["protection_immunity"]);
  epg_constant        = std::exp(theta["epg_constant"]);
  ige_constant        = std::exp(theta["ige_constant"]);
  
  // Define treatment-related parameters
  toggle_tx     = tx_pars["toggle_tx"]; //toggle for treatment on (1) or off (0)
  toggle_covdat = tx_pars["toggle_covdat"]; //toggle for MDA times & coverage informed by data (1) or by user (0)
  start_tx      = tx_pars["start_tx"]; //model time to start treatments
  n_tx          = tx_pars["n_tx"]; //no. treatment events (when toggle_covdat == 0)
  freq_tx       = tx_pars["freq_tx"]; //no. years btwn multiple treatment events (when toggle_covdat == 0)
  coverage      = tx_pars["coverage"]; //treatment coverage (when toggle_covdat == 0)
  efficacy      = tx_pars["efficacy"]; //treatment efficacy
  min_tx_age    = tx_pars["min_tx_age"]; //minimum treatment age 
  max_tx_age    = tx_pars["max_tx_age"]; //maximum treatment age
  cov_weight    = tx_pars["cov_weight"]; //weighting of coverage (0, coverage only in SAC & 1, coverage only in adults)
  
  // if treatment model toggled on
  if (toggle_tx == 1) {
    
    // get treatment times
    treatment_times = Rcpp::clone(tx_times);
    
    // error handling for treatment times
    for (int m = 0; m < treatment_times.size(); m++) {
      if (std::remainder(treatment_times[m], stepsize) > tolerance) {
        Rcpp::stop("Treatment event time is not a multiple of stepsize");
      }
    }
    
  }
  
  // calculate number of time steps
  nsteps = static_cast<int>(runtime/stepsize);  //force nsteps to be integer
  
  // error handling for model run time
  if (std::remainder(runtime, stepsize) > tolerance) { 
    Rcpp::stop("Model stepsize does not divide evenly into model run time");
  }
  
  // storage vector for output (defined here once no. time steps assigned a value)
  NumericVector worm_burden_out(nsteps);
  NumericVector cumulative_worms_out(nsteps);
  NumericMatrix cumulative_worms_age_female_out(nsteps, na), cumulative_worms_age_male_out(nsteps, na);
  NumericMatrix worm_burden_age_female_out(nsteps, na), worm_burden_age_male_out(nsteps, na);
  NumericMatrix male_acquired_immunity_out(nsteps, na), female_acquired_immunity_out(nsteps, na);
  NumericVector time_out(nsteps);
  NumericMatrix dead_worms_male_out(nsteps, na), dead_worms_female_out(nsteps, na);
  NumericMatrix mating_prob_male_out(nsteps, na), mating_prob_female_out(nsteps, na);
  NumericMatrix epg_age_male_out(nsteps, na), epg_age_female_out(nsteps, na), epg_age_male_SAC_out(nsteps, na), epg_age_female_SAC_out(nsteps, na);
  NumericMatrix prevalence_age_male_out(nsteps, na), prevalence_age_female_out(nsteps, na);
  NumericMatrix ige_age_male_out(nsteps, na), ige_age_female_out(nsteps, na);
  NumericVector epg_out(nsteps), epg_SAC_out(nsteps), prevalence_out(nsteps), prevalence_SAC_out(nsteps), ige_out(nsteps), dead_worms_out(nsteps);
  NumericVector kW_out(nsteps), kW_female_out(nsteps), kW_male_out(nsteps);
  NumericMatrix kW_female_age_out(nsteps, na), kW_male_age_out(nsteps, na);

  // initialize k parameters
  kW_female = kW; //all females
  kW_male   = kW; //all males
  wormk     = kW; //overall population
  for (int i = 0; i < na; i++) {
    // age & sex specific
    kW_age_female[i] = kW;
    kW_age_male[i]   = kW;
  }
  
  // run initialization functions
  // age midpoints
  AgeMidpoints();
  // demographic weights
  DemographicWeights();
  // normalized exposure function
  NormalisedExposure();
  // worm burdens, mating probabilities, acquired immunity, epgs, IgEs and dead worms
  InitialiseEverything();
  // initialize state matrices
  SmartInitialiseStates();
  // update outputs based on initial states
  UpdateEverything();
  
  if (toggle_tx == 1) {
    // get age-specific coverage
    age_specific_cov = WeightCoverage(cov_weight, coverage, sac_indices, adult_indices, na);
  }
  
  // store output for time 0
  time_out[0] = 0;
  worm_burden_out[0] = worm_burden;
  cumulative_worms_out[0] = cumulative_worms;
  cumulative_worms_age_female_out(0,_) = cumulative_worms_age_female;
  cumulative_worms_age_male_out(0,_) = cumulative_worms_age_male;
  epg_out[0] = epg;
  epg_SAC_out[0] = epg_SAC;
  prevalence_out[0] = prevalence;
  prevalence_SAC_out[0] = prevalence_SAC;
  ige_out[0] = ige;
  dead_worms_out[0] = dead_worms;
  worm_burden_age_female_out(0,_) = worm_burden_age_female;
  worm_burden_age_male_out(0,_) = worm_burden_age_male;
  male_acquired_immunity_out(0,_) = male_acquired_immunity;
  female_acquired_immunity_out(0,_) = female_acquired_immunity;
  dead_worms_male_out(0,_) = dead_worms_male;
  dead_worms_female_out(0,_) = dead_worms_female;
  mating_prob_female_out(0,_) = mating_probability_female;
  mating_prob_male_out(0,_) = mating_probability_male;
  epg_age_male_out(0,_) = epg_age_male;
  epg_age_female_out(0,_) = epg_age_female;
  epg_age_male_SAC_out(0,_) = epg_age_male_SAC;
  epg_age_female_SAC_out(0,_) = epg_age_female_SAC;
  prevalence_age_male_out(0,_) = prevalence_age_male;
  prevalence_age_female_out(0,_) = prevalence_age_female;
  ige_age_male_out(0,_) = ige_age_male;
  ige_age_female_out(0,_) = ige_age_female;
  kW_female_age_out(0,_)  = kW_age_female;
  kW_male_age_out(0,_)    = kW_age_male;
  kW_out[0] = kW;
  
  // loop through time steps
  for (int h = 0; h < nsteps; h++) {
    
    // update time for output (and use for immunity)
    time_out[h] = h*stepsize;
      
    // perform RK4 integration & store updated state variables
    List updated_states   = RK4(male_states, female_states, stepsize);
    updated_male_states   = as<NumericMatrix>(updated_states["updated_male_states"]);
    updated_female_states = as<NumericMatrix>(updated_states["updated_female_states"]);

    // update state variables
    UpdateStates();
    
    // update worm burden, epg, etc
    UpdateEverything();
    
    // Treatment toggled on
    if (toggle_tx == 1) {
      
      // before first treatment event (-1 is arbitrary)
      if (time_out[h] == treatment_times[0]-1) {
        
        // ensure coverage index initialised to 0
        coverage_index = 0;
        
        // For dynamic k: get worm burden @ equilibrium
        for (int i = 0; i < na; i++) {
          worm_burden_age_female_equib[i] = worm_burden_age_female[i];
          worm_burden_age_male_equib[i]   = worm_burden_age_male[i];
        }
        
        // during treatment period
      } else if (time_out[h] >= treatment_times[0]) {
        
        // check if current time matches a treatment time (using Rcpp sugar syntax)
        if ( is_true(any( Rcpp::abs(time_out[h] - treatment_times) < tolerance )) ) {

          if (toggle_covdat == 0) {
            
            // treat the population
            TreatmentEvent();

          } else {
            // get time-specific coverage
            tx_coverage = coverage_data[coverage_index];
            // treat the population
            TreatmentEvent();
            // update the coverage index for next treatment event 
            coverage_index += 1;
          }
          
          // update states after treatment
          UpdateStates();
          
          // update worm burden, epg, etc after changing the states
          UpdateEverything();
          
          // Calculate dynamic k following treatment event
          DynamickTreatment();
          
        } else if (is_true(all(Rcpp::abs(time_out[h] - treatment_times) >= tolerance))) {

          // Calculate dynamic k when not a treatment event
          DynamickPostTreatment();
          
        } //end of if() where time doesn't match a treatment time
        
      }  //end of if() where times >= treatment times
      
    } //end of if() toggle treatment statement
    
    
    // calculate updated worm burdens & write to output matrices (also for updating acquired immunity)
    // incident worm burden
    worm_burden_out[h] = worm_burden;
    worm_burden_age_female_out(h,_) = worm_burden_age_female;
    worm_burden_age_male_out(h,_) = worm_burden_age_male;
    
    // Calculate & store cumulative worm burden
    CumulativeWormBurden();
    cumulative_worms_out[h] = cumulative_worms;
    cumulative_worms_age_female_out(h,_) = cumulative_worms_age_female;
    cumulative_worms_age_male_out(h,_)   = cumulative_worms_age_male;
    
    // write cumulative dead worms to output matrices (also for updating acquired immunity)
    dead_worms_out[h]          = dead_worms;  //dead worms over time
    dead_worms_male_out(h,_)   = dead_worms_male;
    dead_worms_female_out(h,_) = dead_worms_female;
    
    // calculate updated mating probabilities & write to output matrices (required before epg update)
    mating_prob_female_out(h,_) = mating_probability_female;
    mating_prob_male_out(h,_)   = mating_probability_male;
    
    // calculate updated epgs & write to output matrices
    epg_out[h] = epg;
    epg_SAC_out[h] = epg_SAC;
    epg_age_female_out(h,_) = epg_age_female;
    epg_age_male_out(h,_) = epg_age_male;
    
    // save NBD shape parameter to output
    kW_out[h] = wormk;
    kW_female_age_out(h,_)= kW_age_female;
    kW_male_age_out(h,_)  = kW_age_male;
    kW_female_out[h]= kW_female;
    kW_male_out[h]  = kW_male;
    
    // calculate updated prevalence & write to output matrices
    prevalence_out[h] = prevalence;
    prevalence_SAC_out[h] = prevalence_SAC;
    prevalence_age_female_out(h,_) = prevalence_age_female;
    prevalence_age_male_out(h,_) = prevalence_age_male;
   
    // calculated updated acquired immunity (by direct integration over past dead worms)
    // can change immunological stimulus to current worms by dead_worms_X_out -> worm_burden_age_X
    // required before IgE update
    for (int i = 0; i < na; i++) {
      // set accumulators
      double acc_male_acquired_immunity   = 0;
      double acc_female_acquired_immunity = 0;
      
      for (int z = 0; z < i; z++) {
        // z is an indicator for past age looping up to the present
        // past steps is the time steps in the past to loop over
        // if less than zero loop back over whole current time series
        int past_steps = static_cast<int>(((time_out[h] - (age[i]-age[z]))/stepsize)+(1/stepsize));
        int hh;

        if (past_steps < 0) {
          hh = h;
        } else {
          hh = past_steps;
        }
        // max past steps is h; so either loop back to h or past_steps
        acc_male_acquired_immunity   += exp(-decay_immunity*(age[i]-age[z]))*dead_worms_male_out(hh,z);
        acc_female_acquired_immunity += exp(-decay_immunity*(age[i]-age[z]))*dead_worms_female_out(hh,z);
      }
      male_acquired_immunity[i]   = acc_male_acquired_immunity;
      female_acquired_immunity[i] = acc_female_acquired_immunity;
    }
    
    // write new acquired immunity to output matrices
    male_acquired_immunity_out(h,_)   = male_acquired_immunity;
    female_acquired_immunity_out(h,_) = female_acquired_immunity;
    
    // calculate updated IgE & write to output matrices
    IgE();
    ige_out[h] = ige;  //IgE over time
    ige_age_female_out(h,_) = ige_age_female;
    ige_age_male_out(h,_)   = ige_age_male;
    
  } // end of time loop

  // get value of largest derivative (used to check equilibrium reached ok)
  for (int i = 0; i < na; i++) {
    row_tol_female[i] = max(abs(dWdt_female(i,_)));
    row_tol_male[i]   = max(abs(dWdt_male(i,_)));
  }
  for (int j = 0; j < (nw+1); j++) {
    col_tol_female[j] = max(abs(dWdt_female(_,j)));
    col_tol_male[j]   = max(abs(dWdt_male(_,j)));
  }
  
  tol[0] = max(row_tol_female);
  tol[1] = max(row_tol_male);
  tol[2] = max(col_tol_female);
  tol[3] = max(col_tol_male);
  
  tol_out = max(tol);
  
  // return all outputs over time
  if (alltimes == 1) { 
    //use nested sub-lists as work-around for Rcpp's 20-item list limit
    return Rcpp::List::create(
      Rcpp::Named("age_out") = Rcpp::List::create(
          Rcpp::Named("age") = age,
          Rcpp::Named("epg_age_male") = epg_age_male_out,
          Rcpp::Named("epg_age_female") = epg_age_female_out,
          Rcpp::Named("worm_burden_age_male") = worm_burden_age_male_out,
          Rcpp::Named("worm_burden_age_female") = worm_burden_age_female_out,
          Rcpp::Named("cumulative_worms_age_female") = cumulative_worms_age_female_out,
          Rcpp::Named("cumulative_worms_age_male") = cumulative_worms_age_male_out,
          Rcpp::Named("prevalence_age_male") = prevalence_age_male_out,
          Rcpp::Named("prevalence_age_female") = prevalence_age_female_out,
          Rcpp::Named("IgE_age_male") = ige_age_male_out,
          Rcpp::Named("IgE_age_female") = ige_age_female_out,
          Rcpp::Named("kW_age_female") = kW_female_age_out,
          Rcpp::Named("kW_age_male") = kW_male_age_out
        ), 
        Rcpp::Named("time_out") = Rcpp::List::create(
          Rcpp::Named("time") = time_out,
          Rcpp::Named("worm_burden") = worm_burden_out,
          Rcpp::Named("cumulative_worm_burden") = cumulative_worms_out,
          Rcpp::Named("epg") = epg_out,
          Rcpp::Named("prevalence") = prevalence_out,
          Rcpp::Named("IgE") = ige_out,
          Rcpp::Named("cumulative_dead_worms_female") = dead_worms_female_out,
          Rcpp::Named("cumulative_dead_worms_male") = dead_worms_male_out,
          Rcpp::Named("dead_worms") = dead_worms_out,
          Rcpp::Named("kW_female") = kW_female_out,
          Rcpp::Named("kW_male") = kW_male_out,
          Rcpp::Named("kW") = kW_out
        )
      );
    
    // return outputs only at last timepoint
  } else if (alltimes == 0) {  
    int l = nsteps -1;
    return Rcpp::List::create(
      Rcpp::Named("time") = time_out(l),
      Rcpp::Named("worm_burden") = worm_burden_out(l),
      Rcpp::Named("epg") = epg_out(l),
      Rcpp::Named("worm_burden_age_male") = worm_burden_age_male_out(l,_),
      Rcpp::Named("worm_burden_age_female") = worm_burden_age_female_out(l,_),
      Rcpp::Named("cumulative_worms_age_female") = cumulative_worms_age_female_out(l,_),
      Rcpp::Named("cumulative_worms_age_male") = cumulative_worms_age_male_out(l,_),
      Rcpp::Named("epg_age_male") = epg_age_male_out(l,_),
      Rcpp::Named("epg_age_female") = epg_age_female_out(l,_),
      Rcpp::Named("acquired_immunity_male") = male_acquired_immunity_out(l,_),
      Rcpp::Named("acquired_immunity_female") = female_acquired_immunity_out(l,_),
      Rcpp::Named("IgE_age_male") = ige_age_male_out(l,_),
      Rcpp::Named("IgE_age_female") = ige_age_female_out(l,_),
      Rcpp::Named("cumulative_dead_worms_female") = dead_worms_female_out(l,_),
      Rcpp::Named("cumulative_dead_worms_male") = dead_worms_male_out(l,_),
      Rcpp::Named("age") = age,
      Rcpp::Named("mating_prob_male") = mating_prob_male_out(l,_),
      Rcpp::Named("mating_prob_female") = mating_prob_female_out(l,_),
      Rcpp::Named( "tol") = tol_out);
  }
  
}

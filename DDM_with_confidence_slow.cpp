// in r: install rcpp library and RCPPZiggurat library, then source this script and the function becomes available (source('DDM_KSfitting'))
#include <Rcpp.h>
// [[Rcpp::depends(RcppZiggurat)]]
#include <Ziggurat.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;
//' DDM trials generation with fixed post-decision accumulation time
//' 
//' @description This function generates trials using the Drift Diffusion Model. Additional
//' evidence is accumulated for a fixed time after the decision boundary has been reached.
//' Note: here, decision boundaries are a and -a and starting point is equal to a*z.
//' 
//' @param v Drift rate
//' @param a Decision boundary, evidence at t = 0 is equal to a*z. Decision boundaries are a and -a
//' @param ter Non-decision time
//' @param z Starting point bias in evidence. For unbiased DDM, use z = 0
//' @param ntrials Number of DDM trials to generate
//' @param s Scaling parameter
//' @param dt Time step size
//' @param t2time Post-decision accumulation time
//' @param postdriftmod V-ratio. If equal to 1, post-decision drift rate is equal to decision drift rate
//' @export
// [[Rcpp::export]]
NumericMatrix DDM_with_confidence_slow(double v, double a, double ter, double z, int ntrials, double s, double dt, double t2time, double postdriftmod) {

  // initialize output
  NumericMatrix DATA(ntrials,6);

  // loop over trials
  for (int i = 0; i < ntrials; i++) {

    // initalize variables
    int resp = -1;
    int acc = 0;
    double evidence = a*z;
    double t = 0;

    // Decisional processing
    while (t > -1){

      t = t + dt;
      evidence = evidence + v * dt + s * sqrt(dt) * zigg.norm();

      if (evidence >= a){
        resp = 1;
        evidence=a;
        if (v > 0){
          acc = 1;
        }
        break;
      } else if (evidence <= -1*a) {
        resp = -1;
        evidence=-a;
        if (v < 0){
          acc = 1;
        }
        break;
      }
    }

    DATA(i,0) = (t + ter);
    DATA(i,1) = resp;
    DATA(i,2) = acc;

    //Post-decisional processing
    double v_post = v * postdriftmod;
    for (int j = 0; j < t2time/dt; j++){
      t = t + dt;
      evidence = evidence + v_post * dt + s * sqrt(dt) * zigg.norm();
    }
    DATA(i,3) = evidence;
    DATA(i,4) = (t + ter);
    if (resp == 1){
      DATA(i,5) = evidence-a;
    }
    if(resp == -1){
      DATA(i,5) = -evidence+a;
    }
  }

  return DATA; //RT, resp,accuracy, evidence2, rt2, confidence
}

#include <Rcpp.h>
using namespace Rcpp;

//' Simulate a Multi-Armed Bandit Environment
//'
//' This function simulates a multi-armed bandit environment with Gaussian rewards and non-stationary means.
//'
//' @param N Integer. The number of arms in the bandit.
//' @param T Integer. The number of time steps for the simulation.
//' @param sigma Numeric. The standard deviation of the rewards.
//' @param delta Numeric. The standard deviation of the mean updates.
//' @param mu Numeric vector. The initial means for each arm.
//' @return A numeric matrix containing the simulated rewards for each time step and arm.
//' @export
//' @return A numeric matrix of rewards for each arm across all rounds.
//' @export

 // [[Rcpp::export]]
 NumericMatrix simulate_bandit(int N, int T, double sigma, double delta, NumericVector mu) {

   NumericMatrix rewards(T, N);

   NumericMatrix means(T, N);


   for (int i = 0; i < N; i++) {
     means(0, i) = mu[i];
   }


   for (int t = 0; t < T; t++) {
     for (int i = 0; i < N; i++) {

       rewards(t, i) = R::rnorm(means(t, i), sigma);


       if (t < T - 1) {
         means(t + 1, i) = means(t, i) + R::rnorm(0, delta);
       }
     }
   }


   return rewards;
 }




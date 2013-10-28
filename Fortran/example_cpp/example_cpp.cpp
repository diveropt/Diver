#include <cstdlib>
#include <cmath>
#include <limits>
#include "devo.hpp"

const int         nPar                = 2;                            // Dimensionality of the parameter space
const double      lowerbounds[nPar]   = {-6.0, -6.0};                 // Lower boundaries of parameter space
const double      upperbounds[nPar]   = { 6.0,  6.0};                 // Upper boundaries of parameter space
const char        path[]              = "example_cpp/output/example"; // Path to save samples, resume files, etc 
const int         nDerived            = 0;                            // Number of derived quantities to output
const int         nDiscrete           = 0;                            // Number of parameters that are to be treated as discrete
const int         discrete[]          = {NULL};                       // Indices of discrete parameters, Fortran style, i.e. starting at 1!!
const bool        partitionDiscrete   = false;                        // Split the population evenly amongst discrete parameters and evolve separately
const int         maxciv              = 100;                          // Maximum number of civilisations
const int         maxgen              = 100;                          // Maximum number of generations per civilisation
const int         NP                  = 200;                          // Population size (individuals per generation)
const int         nF                  = 1;                            // Size of the array indicating scale factors
const double      F[nF]               = {0.6};                        // Scale factor(s).  Note that this must be entered as an array.
const double      Cr                  = 0.9;                          // Crossover factor
const double      lambda              = 0.8;                          // Mixing factor between best and rand/current
const bool        current             = false;                        // Use current vector for mutation
const bool        expon               = false;                        // Use exponential crossover
const int         bndry               = 3;                            // Boundary constraint: 1=brick wall, 2=random re-initialization, 3=reflection
const bool        jDE                 = true;                         // Use self-adaptive choices for rand/1/bin parameters as per Brest et al 2006
const bool        lambdajDE           = true;                         // Use self-adaptive rand-to-best/1/bin parameters; based on Brest et al 2006
const bool        removeDuplicates    = true;                         // Weed out duplicate vectors within a single generation
const bool        doBayesian          = true;                         // Calculate approximate log evidence and posterior weightings
const double      maxNodePop          = 1.9;                          // Population at which node is partitioned in binary space partitioning for posterior
const double      Ztolerance          = 1.e-2;                        // Input tolerance in log-evidence
const int         savecount           = 1;                            // Save progress every savecount generations
const bool        resume              = false;                        // Restart from a previous run

const double Pi = 3.14159265359;

//Function to be minimized.  Corresponds to -ln(Likelihood).

//Plain Gaussian centred at the origin, good for any number of dimensions.
double gauss(double params[], const int param_dim, int &fcall, bool &quit, const bool validvector)
{
  double result = 0.0;
  for (int i = 0; i<param_dim; i++) result += params[i]*params[i];
  result += 0.5*nPar*log(Pi);
  if (not validvector) result = std::numeric_limits<double>::max();
  fcall += 1;
  quit = false;
  return result;
}

//Gaussian shells, good for any number of dimensions (just remember to expand the subarrays in c).
double gauss_shell(double params[], const int param_dim, int &fcall, bool &quit, const bool validvector)
{
  double result, temp, dist, loclike;
  int i,j;
  double* greater, lesser; 
  const int nRings = 2;                                      // Number of rings
  const double w[nRings] = {0.1,0.1};                        // Gaussian widths of the shells
  const double r[nRings] = {2.0,2.0};                        // Widths of the rings
  const double c[nRings][nPar] = { {-3.5,0.0}, {3.5,0.0} };  // Positions of ring centres

  result = -std::numeric_limits<double>::max()*1e-5;
  for (i = 0; i < nRings; i++)
  {
    temp = 0.0;
    for (j = 0; j < nPar; j++) temp += (params[j]-c[i][j])*(params[j]-c[i][j]);
    dist = pow(pow(temp,0.5)-r[i], 2);
    loclike = -dist / (2.0*w[i]*w[i]) - 0.5 * log(2.0*Pi*w[i]*w[i]);
    if (result > loclike) result = result + log(1.0 + exp(loclike-result)); 
    else result = loclike + log(1.0 + exp(result-loclike));
  }
  if (not validvector) result = std::numeric_limits<double>::max();
  fcall += 1;
  quit = false;
  return -result;
}

//Flat prior function
double flatprior(const double real_params[], const int real_param_dim)
{
  int result = 1.0;
  for (int i = 0; i < real_param_dim; i++) result *= upperbounds[i] - lowerbounds[i];
  return 1.0/result;
}


//Log prior function
double logprior(const double real_params[], const int real_param_dim)
{
  int result = 1.0;
  for (int i = 0; i < real_param_dim; i++) result /= real_params[i] * log(upperbounds[i]/lowerbounds[i]);
  return result;
}


int main(int argc, char** argv)
{
  runde(gauss, nPar, lowerbounds, upperbounds, path, nDerived, nDiscrete, discrete, partitionDiscrete, 
        maxciv, maxgen, NP, nF, F, Cr, lambda, current, expon, bndry, jDE, lambdajDE, removeDuplicates,
        doBayesian, flatprior, maxNodePop, Ztolerance, savecount, resume); 
  //Note that prior, maxNodePop and Ztolerance are just ignored if doBayesian = false
}
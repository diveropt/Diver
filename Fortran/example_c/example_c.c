#include <float.h>
#include <stdbool.h>
#include <stddef.h>
#include "devo.h"

const int         nPar                = 5;                            // Dimensionality of the parameter space
const double      lowerbounds[]       = {-5.,-50.,-5.,-50.,-2.};      // Lower boundaries of parameter space
const double      upperbounds[]       = { 5., 50., 5., 50., 2.};      // Upper boundaries of parameter space
const char        path[]              = "example_c/output/example";   // Path to save samples, resume files, etc 
const int         nDerived            = 0;                            // Number of derived quantities to output
const int         nDiscrete           = 0;                            // Number of parameters that are to be treated as discrete
const int         discrete[]          = {};                           // Indices of discrete parameters, Fortran style, i.e. starting at 1!!
const bool        partitionDiscrete   = false;                        // Split the population evenly amongst discrete parameters and evolve separately
const int         maxciv              = 1;                            // Maximum number of civilisations
const int         maxgen              = 100;                          // Maximum number of generations per civilisation
const int         NP                  = 1000;                         // Population size (individuals per generation)
const int         nF                  = 1;                            // Size of the array indicating scale factors
const double      F[]                 = {0.6};                        // Scale factor(s).  Note that this must be entered as an array.
const double      Cr                  = 0.9;                          // Crossover factor
const double      lambda              = 0.8;                          // Mixing factor between best and rand/current
const bool        current             = false;                        // Use current vector for mutation
const bool        expon               = false;                        // Use exponential crossover
const int         bndry               = 3;                            // Boundary constraint: 1=brick wall, 2=random re-initialization, 3=reflection
const bool        jDE                 = true;                         // Use self-adaptive choices for rand/1/bin parameters as per Brest et al 2006
const bool        lambdajDE           = true;                         // Use self-adaptive rand-to-best/1/bin parameters; based on Brest et al 2006
const bool        removeDuplicates    = true;                         // Weed out duplicate vectors within a single generation
const bool        doBayesian          = false;                        // Calculate approximate log evidence and posterior weightings
const double      maxNodePop          = 1.9;                          // Population at which node is partitioned in binary space partitioning for posterior
const double      Ztolerance          = 1.e-3;                        // Input tolerance in log-evidence
const int         savecount           = 100;                          // Save progress every savecount generations
const bool        resume              = false;                        // Restart from a previous run


//Function to be minimized.  Corresponds to -ln(Likelihood).
//Plain Gaussian centred at the origin. Valid for any number of dimensions.
double gauss(double params[], const int param_dim, int *fcall, bool *quit, const bool validvector)
{
  double result = 0.0;
  for (int i = 0; i<param_dim; i++) result += params[i]*params[i];
  if (!validvector) result = DBL_MAX;
  *fcall += 1;
  *quit = false;
  return result;
}


int main(int argc, char** argv)
{
  runde(gauss, nPar, lowerbounds, upperbounds, path, nDerived, nDiscrete, discrete, partitionDiscrete, 
        maxciv, maxgen, NP, nF, F, Cr, lambda, current, expon, bndry, jDE, lambdajDE, removeDuplicates,
        doBayesian, NULL, maxNodePop, Ztolerance, savecount, resume); 
  //Note that prior, maxNodePop and Ztolerance are just ignored if doBayesian = false
}
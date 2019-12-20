#include <float.h>
#include <stdbool.h>
#include <stddef.h>
#include "diver.h"

const int         nPar                 = 5;                            // Dimensionality of the parameter space
const double      lowerbounds[]        = {-5.,-50.,-5.,-50.,-2.};      // Lower boundaries of parameter space
const double      upperbounds[]        = { 5., 50., 5., 50., 2.};      // Upper boundaries of parameter space
const char        path[]               = "example_c/output/example";   // Path to save samples, resume files, etc
const int         nDerived             = 0;                            // Number of derived quantities to output
const int         nDiscrete            = 0;                            // Number of parameters that are to be treated as discrete
const int         discrete[]           = {};                           // Indices of discrete parameters, Fortran style, i.e. starting at 1!!
const bool        partitionDiscrete    = false;                        // Split the population evenly amongst discrete parameters and evolve separately
const int         maxciv               = 1;                            // Maximum number of civilisations
const int         maxgen               = 100;                          // Maximum number of generations per civilisation
const int         NP                   = 1000;                         // Population size (individuals per generation)
const int         nF                   = 1;                            // Size of the array indicating scale factors
const double      F[]                  = {0.6};                        // Scale factor(s).  Note that this must be entered as an array.
const double      Cr                   = 0.9;                          // Crossover factor
const double      lambda               = 0.8;                          // Mixing factor between best and rand/current
const bool        current              = false;                        // Use current vector for mutation
const bool        expon                = false;                        // Use exponential crossover
const int         bndry                = 3;                            // Boundary constraint: 1=brick wall, 2=random re-initialization, 3=reflection
const bool        jDE                  = true;                         // Use self-adaptive choices for rand/1/bin parameters as per Brest et al 2006
const bool        lambdajDE            = true;                         // Use self-adaptive rand-to-best/1/bin parameters; based on Brest et al 2006
const double      convthresh           = 1.e-6;                        // Threshold for gen-level convergence: smoothed fractional improvement in the mean population value
const int         convsteps            = 10;                           // Number of steps to smooth over when checking convergence
const bool        removeDuplicates     = true;                         // Weed out duplicate vectors within a single generation
const bool        doBayesian           = false;                        // Calculate approximate log evidence and posterior weightings
const double      maxNodePop           = 1.9;                          // Population at which node is partitioned in binary space partitioning for posterior
const double      Ztolerance           = 1.e-3;                        // Input tolerance in log-evidence
const int         savecount            = 100;                          // Save progress every savecount generations
const bool        resume               = false;                        // Restart from a previous run
const bool        outputSamples        = false;                        // Write output .raw and .sam (if nDerived != 0) files
const int         init_pop_strategy    = 0;                            // Initialisation strategy: 0=one shot, 1=n-shot, 2=n-shot with error if no valid vectors found.
const bool        discard_unfit_points = false;                        // Recalculate any trial vector whose fitness is above max_acceptable_value
const int         max_init_attempts    = 10000;                        // Maximum number of times to try to find a valid vector for each slot in the initial population.
const double      max_acceptable_val   = 1e6;                          // Maximum fitness to accept for the initial generation if init_population_strategy > 0, or any generation if discard_unfit_points = true.
const int         seed                 = 1234567;                      // base seed for random number generation; non-positive or absent means seed from the system clock
const int         verbose              = 1;                            // Output verbosity: 0=only error messages, 1=basic info, 2=civ-level info, 3+=population info


//Function to be minimized.  Corresponds to -ln(Likelihood).
//Plain Gaussian centred at the origin. Valid for any number of dimensions.  Minimum value is the number of dimensions.
double gauss(double params[], const int param_dim, int *fcall, bool *quit, const bool validvector, void** context)
{
  double result = 0.0;
  for (int i = 0; i<param_dim; i++) result += params[i]*params[i] + 1.0;
  if (!validvector) result = DBL_MAX;
  *fcall += 1;
  *quit = false;
  return result;
}


int main(int argc, char** argv)
{
  void* context = &gauss; //Not actually used in this example.
  cdiver(gauss, nPar, lowerbounds, upperbounds, path, nDerived, nDiscrete, discrete, partitionDiscrete,
         maxciv, maxgen, NP, nF, F, Cr, lambda, current, expon, bndry, jDE, lambdajDE, convthresh,
         convsteps, removeDuplicates, doBayesian, NULL, maxNodePop, Ztolerance, savecount, resume,
         outputSamples, init_pop_strategy, discard_unfit_points, max_init_attempts, max_acceptable_val, seed, context, verbose);
         //Note that prior, maxNodePop and Ztolerance are just ignored if doBayesian = false
}

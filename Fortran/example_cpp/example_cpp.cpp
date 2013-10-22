#include <string>
#include <limits>
#include "devo.h"

const int         nPar                = 5;           
const double      lowerbounds[nPar]   = {-5.,-50.,-5.,-50.,-2.};
const double      upperbounds[nPar]   = { 5., 50., 5., 50., 2.};
const char        path[]              = "example_cpp/output/example";
const int         nDerived            = 0;
const int         nDiscrete           = 0;
const int         discrete[]          = {NULL};                    // Indices of discrete parameters -- Fortran style, starting at 1!!
const bool        partitionDiscrete   = false;
const int         maxciv              = 1;
const int         maxgen              = 100;
const int         NP                  = 1000;
const int         nF                  = 1;
const double      F[nF]               = {0.6};
const double      Cr                  = 0.9;
const double      lambda              = 0.8;
const bool        current             = false;
const bool        expon               = false;
const int         bndry               = 3;
const bool        jDE                 = true;
const bool        lambdajDE           = true;
const bool        removeDuplicates    = true;
const bool        doBayesian          = false;
const double      maxNodePop          = 1.9;
const double      Ztolerance          = 1.e-3;
const int         savecount           = 100;
const bool        resume              = false;


//Functions to be minimized.  Assumed to be -ln(Likelihood)

//Plain Gaussian centred at the origin. Valid for any number of dimensions.
double gauss(double params[], const int pSize, int &fcall, bool &quit, const bool validvector)
{
  double result = 0.0;
  for (int i = 0; i<pSize; i++) result += params[i]*params[i];
  if (not validvector) result = std::numeric_limits<double>::max();
  fcall += 1;
  quit = false;
  return result;
}


//Example prior distributions

//Flat prior distribution for all parameters
double flatprior(const double X[], const int pSize) 
{
  double prior_weight = 1.0;
  for (int i = 0; i<nPar; i++) prior_weight *= (upperbounds[i]-lowerbounds[i]);
  prior_weight = 1.0 / prior_weight;
  return prior_weight;
}


int main(int argc, char** argv)
{
  runde(gauss, flatprior, nPar, lowerbounds, upperbounds, path, nDerived, nDiscrete, discrete, 
        partitionDiscrete, maxciv, maxgen, NP, nF, F, Cr, lambda, current, expon, bndry, jDE, 
        lambdajDE, removeDuplicates, doBayesian, maxNodePop, Ztolerance, savecount, resume);
}

// Diver pybind11 bindings (Python -> C++ -> Fortran)
#include "python_binding.hpp"

namespace diver
{

  // pointers to std::functions holding Python objective and prior functions
  func_type* func;
  prior_type* prior;

  // Local redirection functions for objective and prior, to effect various type conversions:
  // - provides plain C-style function pointers for Python Callable types
  // - void* --> Python object, allowing arbitrary Python objects to be passed and later used
  //   as 'context', without faffing about with ctypes module and casting
  // - parameter pack class in place of C-style arrays for model parameters, to allow references
  //   to be passed to them instead of triggering a (potentially expensive) copy operation. We don't
  //   want to opaquely bind std::vector<double>, because we want automatic conversion to list[float]
  //   elsewhere.
  // - unpack return pack instead of modifying types passed by reference that are immutable in Python
  double func_local(double pars[], const int nPar, int& fcall, bool& quit, const bool validvector, void*& context)
  {
    // Statically created parameter pack, so all that needs to be done on each call is to update the pointer to the parameter array.
    static params parameters(pars, nPar);
    parameters.update_pointer(pars);
    std::tuple<double, int, bool> result = (*func)(parameters, fcall, quit, validvector, *reinterpret_cast<py::object*>(context));
    fcall = std::get<1>(result);
    quit = std::get<2>(result);
    return std::get<0>(result);
  }

  double prior_local(const double pars[], const int nPar, void*& context)
  {
    // Statically created parameter pack, so all that needs to be done on each call is to update the pointer to the parameter array.
    // Here we just throw away the constness of the pointer, and rely on the constness of the first argument of the prior std::function type instead.
    static params parameters(const_cast<double*>(pars), nPar);
    parameters.update_pointer(const_cast<double*>(pars));
    return (*prior)(parameters,*reinterpret_cast<py::object*>(context));
  }

  // Local redirection function for the diver main program, using std::vector in place of C-style arrays and explicit size integers,
  // local redirection functions for objective and prior functions, and a return pack instead of output arrays passed by reference.
  std::tuple<double, std::vector<double>, std::vector<double>> diver_cpp(
    func_type func_in,
    std::vector<double> lowerbounds,
    std::vector<double> upperbounds,
    const char path[],
    int nDerived,
    std::vector<int> discrete,
    bool partitionDiscrete,
    int maxciv,
    int maxgen,
    int NP,
    std::vector<double> F,
    double Cr,
    double lambda,
    bool current,
    bool expon,
    int bndry,
    bool jDE,
    bool lambdajDE,
    double convthresh,
    int convsteps,
    bool removeDuplicates,
    bool doBayesian,
    prior_type prior_in,
    double maxNodePop,
    double Ztolerance,
    int savecount,
    bool resume,
    bool disableIO,
    bool outputRaw,
    bool outputSam,
    int init_population_strategy,
    bool discard_unfit_points,
    int max_initialisation_attempts,
    double max_acceptable_value,
    int seed,
    py::object& context_in,
    int verbose )
  {
    func = &func_in;
    prior = &prior_in;
    std::vector<double> bestFitParams(lowerbounds.size());
    std::vector<double> bestFitDerived(nDerived);
    void* context = &context_in;
    double min = cdiver(func_local,
                        lowerbounds.size(),
                        &lowerbounds[0],
                        &upperbounds[0],
                        path,
                        nDerived,
                        &bestFitParams[0],
                        &bestFitDerived[0],
                        discrete.size(),
                        &discrete[0],
                        partitionDiscrete,
                        maxciv,
                        maxgen,
                        NP,
                        F.size(),
                        &F[0],
                        Cr,
                        lambda,
                        current,
                        expon,
                        bndry,
                        jDE,
                        lambdajDE,
                        convthresh,
                        convsteps,
                        removeDuplicates,
                        doBayesian,
                        prior_local,
                        maxNodePop,
                        Ztolerance,
                        savecount,
                        resume,
                        disableIO,
                        outputRaw,
                        outputSam,
                        init_population_strategy,
                        discard_unfit_points,
                        max_initialisation_attempts,
                        max_acceptable_value,
                        seed,
                        context,
                        verbose
                        );
    return { min, bestFitParams, bestFitDerived };
  }

}

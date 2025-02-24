// Diver pybind11 bindings (Python -> C++ -> Fortran)
#include "python_binding.hpp"
    #include <iostream>

namespace diver
{
  // Pointers to std::functions holding Python objective and prior functions
  func_type* func;
  prior_type* prior;

  // Local redirection functions for objective and prior, to effect various type conversions:
  // - provides plain C-style function pointers for Python Callable types
  // - void* --> Python object, allowing arbitrary Python objects to be passed and later used
  //   as 'context', without faffing about with ctypes module and casting
  // - numpy array in place of C-style arrays, providing a view of the underlying array rather than a copy
  // - unpack return pack instead of modifying types passed by reference that are immutable in Python
  double func_local(double pars[], const int nPar, int& fcall, bool& quit, const bool validvector, void*& context)
  {
    py::array_t<double> parameters(
            {nPar},               // shape
            {sizeof(double)},     // stride
            pars,                 // data pointer
            py::capsule([](){})); // base; without passing a base, the array_t constructor copies the underlying data
    std::tuple<double, int, bool> result = (*func)(parameters, fcall, quit, validvector, *reinterpret_cast<py::object*>(context));
    fcall = std::get<1>(result);
    quit = std::get<2>(result);
    return std::get<0>(result);
  }

  double prior_local(const double pars[], const int nPar, void*& context)
  {
    py::array_t<double> parameters(
            {nPar},               // shape
            {sizeof(double)},     // stride
            pars,                 // data pointer
            py::capsule([](){})); // base; without passing a base, the array_t constructor copies the underlying data
    return (*prior)(parameters,*reinterpret_cast<py::object*>(context));
  }

  // Local redirection function for the diver main program, using numpy arrays in place of C-style arrays and explicit size
  // integers, as well as local redirection functions for objective and prior functions.
  std::tuple<double, py::array_t<double>, py::array_t<double>> diver_cpp(
    func_type func_in,
    py::array_t<double>& lowerbounds,
    py::array_t<double>& upperbounds,
    const char path[],
    int nDerived,
    py::array_t<int>& discrete,
    bool partitionDiscrete,
    int maxciv,
    int maxgen,
    int NP,
    py::array_t<double>& F,
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
    // Get data pointer and data size for the lowerbounds array
    py::buffer_info info = lowerbounds.request();
    double* lowerbounds_ptr = static_cast<double*>(info.ptr);
    int nPar = info.shape[0];
    // Get data pointers and sizes for other input arrays
    double* upperbounds_ptr = static_cast<double*>(upperbounds.request().ptr);
    info = discrete.request();
    int* discrete_ptr = static_cast<int*>(info.ptr);
    int nDiscrete = info.shape[0];
    info = F.request();
    double* F_ptr = static_cast<double*>(info.ptr);
    int nF = info.shape[0];
    // Create and get data pointers to output arrays.
    py::array_t<double> bestFitParams({nPar}, {sizeof(double)});
    double* bestFitParams_ptr = static_cast<double*>(bestFitParams.request().ptr);
    py::array_t<double> bestFitDerived({nDerived}, {sizeof(double)});
    double* bestFitDerived_ptr = static_cast<double*>(bestFitDerived.request().ptr);

    func = &func_in;
    prior = &prior_in;
    void* context = &context_in;

    double min = cdiver(func_local,
                        nPar,
                        lowerbounds_ptr,
                        upperbounds_ptr,
                        path,
                        nDerived,
                        bestFitParams_ptr,
                        bestFitDerived_ptr,
                        nDiscrete,
                        discrete_ptr,
                        partitionDiscrete,
                        maxciv,
                        maxgen,
                        NP,
                        nF,
                        F_ptr,
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
    return { min, std::move(bestFitParams), std::move(bestFitDerived) };
  }

}

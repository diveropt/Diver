#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <functional>
#include "diver.hpp"
#pragma once

namespace py = pybind11;

namespace diver
{

  /// Parameter pack class, for exposing explicitly to Python in order to allow parameters to be accessesd without copies, and
  /// without needing to bind std::vector<double> opaquely.
  class params
  {
    public:
      /// Constructor
      params(double arr[], long int arr_size) : wrapped_array(&arr[0]), wrapped_array_size(arr_size) {}
      /// Element accessor
      double operator[] (long int i) const { return wrapped_array[i]; }
      /// Element modifier
      double& operator[](long int i) { return wrapped_array[i]; }
      /// Length
      long int size() const { return wrapped_array_size; }
      /// Update wrapped array
      void update_pointer(double arr[]) { wrapped_array = &arr[0]; }

    private:
      /// Pointer to the first element of the wrapped C-style parameter array
      double* wrapped_array;
      /// Size of wrapped array
      int wrapped_array_size;
  };

  /// Std::function wrapper types for objective and prior functions, to allow automatic conversion from Python Callable type
  typedef std::function<std::tuple<double, int, bool>(params&, int, bool, bool, py::object&)> func_type;
  typedef std::function<double(const params&, py::object&)> prior_type;

  /// C++ prototype of main Diver function, specifically engineered for exposing to Python via pybind11.
  /// You could call this from C++ too though, if you prefer its signature to the C-style signature of cdiver (in diver.hpp).
  std::tuple<double, std::vector<double>, std::vector<double>> diver_cpp(func_type, std::vector<double>, std::vector<double>,
   const char[], int, std::vector<int>, bool, int, int, int, std::vector<double>, double, double, bool,
   bool, int, bool, bool, double, int, bool, bool, prior_type, double, double, int, bool, bool, bool, bool, int, bool, int,
   double, int, py::object&, int);

}

PYBIND11_MODULE(diver_cpp, m) {
  m.doc() = "Diver differential evolution C++ interface";
  py::class_<diver::params>(m, "params")
   .def("__getitem__", [](const diver::params& p, long int i) { return p[i]; })
   .def("__setitem__", [](diver::params& p, long int i, double val) { p[i] = val; })
   .def("__len__", [](const diver::params& p) { return p.size(); })
   .def("__iter__", [](diver::params& p) { return py::make_iterator(&p[0], &p[p.size()]); },
    py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */
  m.def("run", &diver::diver_cpp, "Run differential evolution via C++ interface.",
        py::arg("func"),
        py::arg("lowerbounds"),
        py::arg("upperbounds"),
        py::arg("path"),
        py::arg("nDerived"),
        py::arg("discrete"),
        py::arg("partitionDiscrete"),
        py::arg("maxciv"),
        py::arg("maxgen"),
        py::arg("NP"),
        py::arg("F"),
        py::arg("Cr"),
        py::arg("lmbda"),
        py::arg("current"),
        py::arg("expon"),
        py::arg("bndry"),
        py::arg("jDE"),
        py::arg("lambdajDE"),
        py::arg("convthresh"),
        py::arg("convsteps"),
        py::arg("removeDuplicates"),
        py::arg("doBayesian"),
        py::arg("prior"),
        py::arg("maxNodePop"),
        py::arg("Ztolerance"),
        py::arg("savecount"),
        py::arg("resume"),
        py::arg("disableIO"),
        py::arg("outputRaw"),
        py::arg("outputSam"),
        py::arg("init_population_strategy"),
        py::arg("discard_unfit_points"),
        py::arg("max_initialisation_attempts"),
        py::arg("max_acceptable_value"),
        py::arg("seed"),
        py::arg("context"),
        py::arg("verbose")
  );
}


#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <functional>
#include "diver.hpp"
#pragma once

namespace py = pybind11;

namespace diver
{

  /// Std::function wrapper type for objective function, to allow automatic conversion from Python Callable type
  typedef std::function<std::tuple<double, int, bool>(py::array_t<double>&, int, bool, bool, py::object&)> func_type;

  /// C++ prototype of main Diver function, specifically engineered for exposing to Python via pybind11.
  /// You could call this from C++ too though, if you prefer its signature to the C-style signature of cdiver (in diver.hpp).
  std::tuple<double, py::array_t<double>, py::array_t<double>> diver_cpp(func_type, py::array_t<double>&, py::array_t<double>&,
   const char[], int, py::array_t<int>&, bool, int, int, py::array_t<double>&, double, double, bool, bool, int, bool, bool,
   double, int, bool, int, bool, bool, bool, bool, int, py::array_t<double>&, bool, int, double, int, py::object&, int);

}

PYBIND11_MODULE(diver_cpp, m) {
  m.doc() = "Diver differential evolution C++ interface";
  m.def("run", &diver::diver_cpp, "Run differential evolution via C++ interface.",
        py::arg("func"),
        py::arg("lowerbounds"),
        py::arg("upperbounds"),
        py::arg("path"),
        py::arg("nDerived"),
        py::arg("discrete"),
        py::arg("partitionDiscrete"),
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
        py::arg("savecount"),
        py::arg("resume"),
        py::arg("disableIO"),
        py::arg("outputRaw"),
        py::arg("outputSam"),
        py::arg("init_population_strategy"),
        py::arg("initial_guesses"),
        py::arg("discard_unfit_points"),
        py::arg("max_initialisation_attempts"),
        py::arg("max_acceptable_value"),
        py::arg("seed"),
        py::arg("context"),
        py::arg("verbose")
  );
}


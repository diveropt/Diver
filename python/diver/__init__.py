"""
Diver differential evolution
"""

import numpy as np
from numpy.typing import NDArray
from collections import abc
from collections.abc import Callable
from typing import TypedDict
from diver_cpp import run as run_cpp

class options(TypedDict):
  """Diver options typed dictionary"""
  lowerbounds: NDArray[np.float64]
  upperbounds: NDArray[np.float64]
  path: str
  nDerived: int
  discrete: NDArray[np.int32]
  partitionDiscrete: bool
  maxciv: int
  maxgen: int
  NP: int
  F: NDArray[np.float64]
  Cr: float
  lmbda: float
  current: bool
  expon: bool
  bndry: int
  jDE: bool
  lambdajDE: bool
  convthresh: float
  convsteps: int
  removeDuplicates: bool
  doBayesian: bool
  prior: Callable[[NDArray[np.float64], object], float]
  maxNodePop: float
  Ztolerance: float
  savecount: int
  resume: bool
  disableIO: bool
  outputRaw: bool
  outputSam: bool
  init_population_strategy: int
  discard_unfit_points: bool
  max_initialisation_attempts: int
  max_acceptable_value: float
  seed: int
  context: object
  verbose: int

def dummy_prior(p: NDArray[np.float64], context: object) -> float:
  """Flat dummy prior."""
  return 1.0

def defaults(lowerbounds, upperbounds):
  """Return a Diver options typed dictionary populated with passed values for required options and defaults for other options."""
  d: options
  d = {'lowerbounds': lowerbounds,
       'upperbounds': upperbounds,
       'path': 'output',
       'nDerived': 0,
       'discrete': np.array([], dtype=np.float64),
       'partitionDiscrete': False,
       'maxciv': 1,
       'maxgen': 300,
       'NP': max(10*len(upperbounds), 5),
       'F': np.array([0.7]),
       'Cr': 0.9,
       'lmbda': 0.0,
       'current': False,
       'expon': False,
       'bndry': 1,
       'jDE': True,
       'lambdajDE': True,
       'convthresh': 1e-3,
       'convsteps': 10,
       'removeDuplicates': True,
       'doBayesian': False,
       'prior': dummy_prior,
       'maxNodePop': 1.9,
       'Ztolerance': 0.01,
       'savecount': 1,
       'resume': False,
       'disableIO': False,
       'outputRaw': True,
       'outputSam': True,
       'init_population_strategy': 0,
       'discard_unfit_points': False,
       'max_initialisation_attempts': 10000,
       'max_acceptable_value': 1e6,
       'seed': -1,
       'context': None,
       'verbose': 1
       }
  return d

def run(func: Callable[[NDArray[np.float64], int, bool, bool, object], tuple[float, int, bool]], opt: options) -> tuple[float, NDArray[np.float64], NDArray[np.float64]]:
  """
  Run Diver differential evolution.

  Args:
      func [Callable]: A callable object that returns the value of the function to be optimised.
          Args:
              params [numpy.NDArray[np.float64]]: parameter values at which to evaluate the function.
              fcall [int]: number of calls to func so far.
              finish [bool]: True if diver should halt at the end of the current generation.
              validvector [bool]: True if the parameters passed are considered valid (inside upperbounds and lowerbounds options passed when invoking diver.run).
              context [object]: Python object passed as an option when invoking diver.run.
          Returns [tuple[float, int, bool]]:
              function value [float]
              updated fcall [int]
              updated validvector [bool]
      opt [diver.options]: A typed dictionary containing the names and values all options accepted by diver.

  Returns [tuple[float, NDArray[np.float64], NDArray[np.float64]]]:
      minimum function value found [float]
      parameters at minimum [NDArray[np.float64]]
      derived parameters at minimum [NDArray[np.float64]]
  """
  return run_cpp(func=func, **opt)

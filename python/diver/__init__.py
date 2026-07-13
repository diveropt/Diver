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
  """
  Diver options typed dictionary.

  Fields:
    lowerbounds [NDArray[np.float64]]: lower boundaries of parameter space.
    upperbounds [NDArray[np.float64]]: upper boundaries of parameter space.
    path [str]: path to save samples, resume files, etc.
    nDerived [int]: number of derived quantities to output.
    discrete [NDArray[np.int32]]: a 1D array listing all discrete dimensions of parameter space.
    partitionDiscrete [bool]: split the population evenly amongst discrete parameters and evolve separately.
    maxgen [int]: maximum number of generations.
    NP [int]: population size (individuals per generation).
    F [NDArray[np.float64]]: scale factor(s).
    Cr [float]: crossover factor.
    lmbda [float]: mixing factor between best and rand/current.
    current [bool]: use current vector for mutation.
    expon [bool]: use exponential crossover.
    bndry [int]: boundary constraint: 1 -> brick wall, 2 -> random re-initialization, 3 -> reflection.
    jDE [bool]: use self-adaptive choices for rand/1/bin parameters as described in Brest et al 2006.
    lambdajDE [bool]: use self-adaptive choices for rand-to-best/1/bin parameters; based on Brest et al 2006.
    convthresh [float]: threshold for generation-level convergence.
    convsteps [int]: number of steps to smooth over when checking convergence.
    removeDuplicates [bool]: weed out duplicate vectors within a single generation.
    savecount [int]: save progress every savecount generations.
    resume [bool]: restart from a previous run.
    disableIO [bool]: disable all I/O.
    outputRaw [bool]: output raw parameter samples to a .raw file.
    outputSam [bool]: output rounded and derived parameter samples to a .sam file.
    init_population_strategy [int]: initialisation strategy: 0=one shot, 1=n-shot, 2=n-shot with error if no valid vectors found.
    initial_guesses [NDArray[np.float64]]: Initial guesses to include in the starting population. 2D array. In Python the first index is parameter, second is individual.
    discard_unfit_points [bool]: recalculate any trial vector whose fitness is above max_acceptable_value. Likely incompatible with any objective function that makes MPI calls of its own.
    max_initialisation_attempts [int]: maximum number of times to try to find a valid vector for each slot in the initial population.
    max_acceptable_value [float]: maximum fitness to accept for the initial generation if init_population_strategy > 0. Also applies to later generations if discard_unfit_points = .true.
    seed [int]: base seed for random number generation; non-positive or absent means seed from the system clock.
    context [object]: context object, used for passing info from the caller to the objective. Use this for passing a callback object that can be used for I/O, harvesting samples in situ, printing or whatever else you like.
    verbose [int]: output verbosity: 0=only error messages, 1=basic info, 2+=population info.
  """
  lowerbounds: NDArray[np.float64]
  upperbounds: NDArray[np.float64]
  path: str
  nDerived: int
  discrete: NDArray[np.int32]
  partitionDiscrete: bool
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
  savecount: int
  resume: bool
  disableIO: bool
  outputRaw: bool
  outputSam: bool
  init_population_strategy: int
  initial_guesses: NDArray[np.float64]
  discard_unfit_points: bool
  max_initialisation_attempts: int
  max_acceptable_value: float
  seed: int
  context: object
  verbose: int

def defaults(lowerbounds, upperbounds):
  """Return a Diver options typed dictionary populated with passed values for required options and defaults for other options."""
  d: options
  d = {'lowerbounds': lowerbounds,
       'upperbounds': upperbounds,
       'path': 'output',
       'nDerived': 0,
       'discrete': np.array([], dtype=np.int32),
       'partitionDiscrete': False,
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
       'savecount': 1,
       'resume': False,
       'disableIO': False,
       'outputRaw': True,
       'outputSam': True,
       'init_population_strategy': 0,
       'initial_guesses': np.empty((len(upperbounds), 0), dtype=np.float64),
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
      parameter values at minimum [NDArray[np.float64]]
      values of derived quantities at minimum [NDArray[np.float64]]
  """
  return run_cpp(func=func, **opt)

"""
Diver differential evolution
"""

from typing import TypedDict, Callable
from diver_cpp import params
from diver_cpp import run as run_cpp

class options(TypedDict):
  """Diver options typed dictionary"""
  lowerbounds: list[float]
  upperbounds: list[float]
  path: str
  nDerived: int
  discrete: list[int]
  partitionDiscrete: bool
  maxciv: int
  maxgen: int
  NP: int
  F: list[float]
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
  prior: Callable[[params, object], float]
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

def dummy_prior(p: params, context: object) -> float:
  """Flat dummy prior."""
  return 1.0

def defaults(lowerbounds, upperbounds):
  """Return a Diver options typed dictionary populated with passed values for required options and defaults for other options."""
  d: options
  d = {'lowerbounds': lowerbounds,
       'upperbounds': upperbounds,
       'path': 'output',
       'nDerived': 0,
       'discrete': [],
       'partitionDiscrete': False,
       'maxciv': 1,
       'maxgen': 300,
       'NP': max(10*len(upperbounds), 5),
       'F': [0.7],
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

def run(func: Callable[[params, int, bool, bool, object], tuple[float, int, bool]], opt: options) -> tuple[float, list[float], list[float]]:
  """
  Run Diver differential evolution.
  """
  return run_cpp(func=func, **opt)

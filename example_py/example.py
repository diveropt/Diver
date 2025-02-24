# Simple python example for Diver.  'Just works' either serially as
#   python3 example.py
# or using e.g. 4 MPI processes as
#   mpirun -n 4 python3 example.py
import diver
import numpy as np
from numpy.typing import NDArray

def gauss(p: NDArray[np.float64], fcall: int, finish: bool, validvector: bool, context: object) -> tuple[float, int, bool]:
  """ Plain Gaussian centred at the origin. Valid for any number of dimensions.  Minimum value is the number of dimensions."""
  finish = False
  objective = 1e300 if not validvector else np.sum(p**2 + 1)
  return objective, fcall+1, finish

def main():
  N=10
  opts = diver.defaults(upperbounds=[2]*N, lowerbounds=[-2]*N)
  s = diver.run(gauss, opts)
  print(f"Min in {N:>2d}D: ", s[0])
  print("Found at: ", s[1])

if __name__=="__main__":
    main()


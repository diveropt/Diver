# Simple python example for Diver.  'Just works' either serially as
#   python3 example.py
# or using e.g. 4 MPI processes as
#   mpirun -n 4 python3 example.py
import diver

def gauss(p: diver.params, fcall: int, finish: bool, validvector: bool, context: object) -> tuple[float, int, bool]:
  """ Plain Gaussian centred at the origin. Valid for any number of dimensions.  Minimum value is the number of dimensions."""
  finish = False;
  objective = 1e300 if not validvector else sum([x*x+1 for x in p])
  return objective, fcall+1, finish

def main():
  # 2D
  opts = diver.defaults(upperbounds=[2,2], lowerbounds=[-2,-2])
  s = diver.run(gauss, opts)
  print("Min in 2D: ", s[0])
  print("Found at: ", s[1])

if __name__=="__main__":
    main()


C++ version of DE code by A. Putze

=== Installation ===
I added autoconfigure and automake files in order to ensure the transportability of the code. To configure and compile the code follow these 3 steps:
1. autoreconf -ivf
2. ./configure (optional: --prefix==/Path/to/Installation)
3. make
Optional:
4. make install (requires that /Path/to/Installation is in your LD_LIBRARY_PATH for the executable to run)

=== Run the code ===
For now no input is needed for this version of the code and you can run the executable directly:
   /Path/to/Installation/bin/DE_example

=== Comments ===
5 Different likelihood functions are implemented in DE_C++/examples/include/TestFunctions.h. Plots for these functions are saved as .pdfs in DE_C++/plots. If you want to change the function you have to edit line 25 in DE_C++/examples/DE_example.cpp, i.e. you need to put the class name of the wanted likelihoodfunction, and recompile the code.

The code depends on no other library. This has one minor disadvantage: the random number generator is pretty bad.

I implemented one Mutation, Crossover and Selection class in order to have a simple DE/rand/1/bin. To implement other versions of DE, one has just to write the given class and use it as a template argument when creating the DE object.

=== TO DO ===
- add a better random generator
- comment the code
- add other mutation, crossover, and selection criteria
- add a convergence criterium
- ...

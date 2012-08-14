#include "Trial.h"
#include "DE.h"
#include "Mutation.h"
#include "Crossover.h"
#include "Selection.h"
#include "TestFunctions.h"
#include <algorithm>

#ifdef HAVE_MPI
	#include "mpi.h"
#endif

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
	MPI::Status status;
	// Initialise the MPI environment
	MPI::Init();

	// Get the number of processes
	int world_size = MPI::COMM_WORLD.Get_size(); 
	
	// Get the rank of the process
	int world_rank = MPI::COMM_WORLD.Get_rank();
#endif

	// Defining DE with a trial type, mutation scheme, crossover scheme, and selection scheme. The construction takes the number of generations and individuals per generation.
	DE<Trial, Mutation<Trial>, Crossover<Trial>, Selection> *randbest2bin = new DE<Trial, Mutation<Trial>, Crossover<Trial>, Selection>(20, 100);
	// Defining parameters of our problem: name, lower bound, upper bound
	randbest2bin -> GetTrial() -> AddParameter("x", 0., 1.);
	randbest2bin -> GetTrial() -> AddParameter("y", 0., 1.);
	// Refining mutation scheme
	randbest2bin -> GetMutation() -> SetNDifferenceVectors(2);
	randbest2bin -> GetMutation() -> SetScalingFactor(0, 1.);
	randbest2bin -> GetMutation() -> SetScalingFactor(1, 0.7);
	randbest2bin -> GetMutation() -> SetLambda(.5);
	// Refining crossover scheme
	randbest2bin -> GetCrossover() -> SetCrossoverRate(0.9);

	// Print out DE with all its characteristics
	std::cout << *randbest2bin << std::endl;
	
	// Define a likelihood function
	GaussianHill likelihoodFnc;
	
	// Run DE with the above define likelihood function
	randbest2bin -> Run(likelihoodFnc);
	
	// Freeing memory
	delete randbest2bin;
	randbest2bin = NULL;
	
#ifdef HAVE_MPI
	// Finalise the MPI environment.
	MPI::Finalize();
#endif
	return 0;
}

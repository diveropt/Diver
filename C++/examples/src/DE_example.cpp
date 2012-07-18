#include "Trial.h"
#include "DE.h"
#include "Mutation.h"
#include "Crossover.h"
#include "Selection.h"
#include "TestFunctions.h"
#include <algorithm>

int main(int argc, char** argv)
{	
	DE<Trial, Mutation, Crossover, Selection> *randbest2bin = new DE<Trial, Mutation, Crossover, Selection>(20, 100);
	randbest2bin -> GetTrial() -> AddParameter("x", 0., 1.);
	randbest2bin -> GetTrial() -> AddParameter("y", 0., 1.);
	randbest2bin -> GetMutation() -> SetNDifferenceVectors(2);
	randbest2bin -> GetMutation() -> SetScalingFactor(0, 1.);
	randbest2bin -> GetMutation() -> SetScalingFactor(1, 0.7);
	randbest2bin -> GetMutation() -> SetLambda(.5);
	randbest2bin -> GetCrossover() -> SetCrossoverRate(0.9);

	std::cout << *randbest2bin << std::endl;
	
	GaussianHill likelihoodFnc;
	
	randbest2bin -> Run(likelihoodFnc);
	
	delete randbest2bin;
	randbest2bin = NULL;

	return 0;
}

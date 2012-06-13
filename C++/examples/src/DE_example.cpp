#include "DE.h"
#include "Mutation.h"
#include "Crossover.h"
#include "Selection.h"
#include "TestFunctions.h"
#include <algorithm>

int main(int argc, char** argv)
{
	Model *MyModel = new Model();
	MyModel -> AddParameter("x", 0., 1.);
	MyModel -> AddParameter("y", 0., 1.);
	
	std::cout << *MyModel << std::endl;
		
	DE<Mutation, Crossover, Selection> *randbest2bin = new DE<Mutation, Crossover, Selection>(20, 100, MyModel);
	randbest2bin -> GetMutation() -> SetNDifferenceVectors(2);
	randbest2bin -> GetMutation() -> SetScalingFactor(0, 1.);
	randbest2bin -> GetMutation() -> SetScalingFactor(1, 0.7);
	randbest2bin -> GetMutation() -> SetLambda(.5);
	randbest2bin -> GetCrossover() -> SetCrossoverRate(0.9);

	GaussianHill likelihoodFnc;
	
	randbest2bin -> Run(likelihoodFnc);
	
	delete randbest2bin;
	randbest2bin = NULL;
	
	delete MyModel;
	MyModel = NULL;
	
	return 0;
}

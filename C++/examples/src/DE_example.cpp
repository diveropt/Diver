#include "DE.h"
#include "Mutation.h"
#include "Crossover.h"
#include "Selection.h"
#include "TestFunctions.h"

int main()
{
	Model *MyModel = new Model();
	MyModel -> AddParameter("x", 0., 1.);
	MyModel -> AddParameter("y", 0., 1.);
	
	std::cout << *MyModel << std::endl;
		
	DE<Mutation, Crossover, Selection> *rand1bin = new DE<Mutation, Crossover, Selection>(MyModel);
	rand1bin -> GetMutation() -> SetScalingFactor(1.);
	rand1bin -> GetCrossover() -> SetCrossoverRate(0.9);
				
	const int nPop = 20;
	const int nGen = 100;
	
	std::vector<std::vector<Trial*> > triallist(nGen + 1);
	
	StairCase likelihoodFnc;
	
	std::cout << " ==> First initialisation" << std::endl;
	for(int i = 0; i < nPop; i++)
	{
		triallist.at(0).push_back(rand1bin -> FirstTrial(likelihoodFnc));
	}

	std::cout << " ** The best individual ** \n" << *(rand1bin -> GetSelection() -> SelectBest(triallist.at(0))) << std::endl;
	
	for(int i = 0; i < nGen; i++)
	{
		std::cout << "+++ Generation " << i << std::endl;

		for(int j = 0; j < nPop; j++)
		{
			triallist.at(i+1).push_back(rand1bin -> NextTrial(triallist.at(i), j, likelihoodFnc));
		}
		
		std::cout << " ** The best individual ** \n" << *(rand1bin -> GetSelection() -> SelectBest(triallist.at(i+1))) << std::endl;
	}

	// Freeing memory
	std::vector<Trial*> cleanup;
	std::vector<Trial*>::iterator it;
	for(int i = 0; i <= nGen; i++)
	{
		for(int j = 0; j < nPop; j++)
			cleanup.push_back(triallist.at(i).at(j));
	}
	
	sort(cleanup.begin(), cleanup.end());
	it = unique(cleanup.begin(), cleanup.end());
	cleanup.resize(it - cleanup.begin());
	for(it = cleanup.begin(); it != cleanup.end(); it++)
		delete *it;
	
	cleanup.clear();
		
	delete rand1bin;
	rand1bin = NULL;
	
	delete MyModel;
	MyModel = NULL;
	
	return 0;
}
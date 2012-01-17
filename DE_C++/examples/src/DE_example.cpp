#include "DE.h"
#include "Model.h"
#include "Mutation.h"
#include "Crossover.h"
#include "Selection.h"
#include "TestFunctions.h"

int main()
{
	Model *MyModel = new Model();
	MyModel -> AddParameter("x", 0., 1.);
	MyModel -> AddParameter("y", 0., 1.);
	
	MyModel -> Show();
		
	DE<Mutation, Crossover, Selection> *rand1bin = new DE<Mutation, Crossover, Selection>(MyModel);
	rand1bin -> GetMutation() -> SetScalingFactor(1.);
	rand1bin -> GetCrossover() -> SetCrossoverRate(0.9);
				
	const int nPop = 20;
	const int nGen = 100;
	
	std::vector<std::vector<Trial*> > triallist = std::vector<std::vector<Trial*> >(nGen + 1);
	
	StairCase likelihoodFnc;
	
	std::cout << " ==> First initialisation" << std::endl;
	for(int i = 0; i < nPop; i++)
	{
		triallist.at(0).push_back(new Trial(*rand1bin -> FirstTrial(likelihoodFnc)));
//		triallist.at(0).at(i) -> Show();
	}

	std::cout << " ** The best individual ** " << std::endl;
	rand1bin -> GetSelection() -> SelectBest(triallist.at(0)) -> Show();
	
	for(int i = 0; i < nGen; i++)
	{
		std::cout << "+++ Generation " << i << std::endl;

		for(int j = 0; j < nPop; j++)
		{
			triallist.at(i+1).push_back(new Trial(*rand1bin -> NextTrial(triallist.at(i), j, likelihoodFnc)));
//			triallist.at(i+1).at(j) -> Show();
		}
		
		std::cout << " ** The best individual ** " << std::endl;
		rand1bin -> GetSelection() -> SelectBest(triallist.at(i+1)) -> Show();
	}

	for(int i = 0; i < nGen + 1; i++)
	{
		for(int j = 0; j < nPop; j++)
			delete triallist.at(i).at(j);

		triallist.at(i).clear();
	}
	
	triallist.clear();
	
	return 0;
}
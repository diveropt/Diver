#include "Selection.h"

Trial* Selection::Select(Trial *targetvector, Trial *trialvector)
{
	if(targetvector -> GetFitness() < trialvector -> GetFitness())
		return targetvector;
	else
		return trialvector;
}

Trial* Selection::SelectBest(std::vector<Trial*> population)
{
	int iBest = -1;
	double fBest = 1.e99;
	
	for(unsigned int i = 0; i < population.size(); i++)
	{
		if(population.at(i) -> GetFitness() < fBest)
		{
			fBest = population.at(i) -> GetFitness();
			iBest = i;
		}
	}

	return population.at(iBest);
}

void Selection::Show()
{
	std::cout << "=== Selection ===" << std::endl;
}
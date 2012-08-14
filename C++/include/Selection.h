#ifndef __SELECTION__H
#define __SELECTION__H

#include "global.h"
#include <iostream>

class Selection
{	
public:
	Selection() {}
	
	// Other member functions
	template<typename TrialT>
	TrialT* Select(TrialT *targetvector, TrialT *trialvector);
	
	template<typename TrialT>
	TrialT* SelectBest(std::vector<TrialT*> population);
	
	void Show(std::ostream &s = std::cout);
};

template<typename TrialT>
TrialT* Selection::Select(TrialT *targetvector, TrialT *trialvector)
{
	if(targetvector -> GetFitness() < trialvector -> GetFitness())
		return targetvector;
	else
		return trialvector;
}

template<typename TrialT>
TrialT* Selection::SelectBest(std::vector<TrialT*> population)
{
	typename std::vector<TrialT*>::iterator it_best = std::min_element(population.begin(), population.end(), CompareFitness<TrialT>);
	
	return *it_best;
}

std::ostream & operator<<(std::ostream &s, Selection &selection);

#endif // __SELECTION__H
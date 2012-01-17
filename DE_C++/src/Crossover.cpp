#include "Crossover.h"
#include <cstdlib>
#include <ctime>

Crossover& Crossover::operator=(Crossover crossover)
{
	std::swap(fCrossoverRate, crossover.fCrossoverRate);
	std::swap(fTrialVector, crossover.fTrialVector);

	srand(time(0));
}

void Crossover::CrossOver(const std::vector<double> &targetvector, const std::vector<double> &donorvector)
{
	if(targetvector.size() != donorvector.size())
	{
		std::cout << "Target and donor vector don't have the same size (" << targetvector.size() << " vs. " << donorvector.size() << ")" << std::endl;
		exit(1);
	}

	fTrialVector.resize(targetvector.size());

	double randomnumber = 0.;
	int randominteger = -1;
	
	for(unsigned int i = 0; i < targetvector.size(); i++)
	{
		randomnumber = (double)rand()/((double) RAND_MAX + 1.);
		randominteger = rand()%targetvector.size();
		
		if(randomnumber <= fCrossoverRate || i == randominteger)
			fTrialVector.at(i) = donorvector.at(i);
		else
			fTrialVector.at(i) = targetvector.at(i);
	}
}

void Crossover::Show()
{
	std::cout << "=== CrossOver ===" << std::endl;
	std::cout << "Cross over rate: " << fCrossoverRate << std::endl;
	std::cout << "Trial vector: ";
	if(fTrialVector.size() == 0)
		std::cout << "empty" << std::endl;
	else
	{
		for(unsigned int i = 0; i < fTrialVector.size(); i++)
			std::cout << fTrialVector.at(i) << "\t";
		std::cout << std::endl;
	}
}
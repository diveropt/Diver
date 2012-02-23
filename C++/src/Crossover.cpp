#include "Crossover.h"

#include "global.h"

Crossover& Crossover::operator=(Crossover crossover)
{
	std::swap(fCrossoverRate, crossover.fCrossoverRate);
	std::swap(fTrialVector, crossover.fTrialVector);

	srand(time(0));
	
	return *this;
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
	
	for(::vector_size_t i = 0; i < targetvector.size(); i++)
	{
		randomnumber = static_cast<double> (rand())/(RAND_MAX + 1.);
		randominteger = rand()%targetvector.size();
		
		if(randomnumber <= fCrossoverRate || i == static_cast< ::vector_size_t > (randominteger))
			fTrialVector.at(i) = donorvector.at(i);
		else
			fTrialVector.at(i) = targetvector.at(i);
	}
}

void Crossover::Show(std::ostream &s)
{
	s << "=== CrossOver ===" << std::endl;
	s << "Cross over rate: " << fCrossoverRate << std::endl;
	s << "Trial vector: ";
	if(fTrialVector.size() == 0)
		s << "empty" << std::endl;
	else
	{
		for(::vector_size_t i = 0; i < fTrialVector.size(); i++)
			s << fTrialVector.at(i) << "\t";
		s << std::endl;
	}
}

std::ostream & operator<<(std::ostream &s, Crossover &crossover)
{
	crossover.Show(s);
	
	return s;
}
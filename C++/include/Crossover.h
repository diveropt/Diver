#ifndef __CROSSOVER__H
#define __CROSSOVER__H

#include <iostream>
#include <vector>
#include <tr1/random>
#include "global.h"

template <typename TrialT>
class Crossover
{
	double fCrossoverRate;
	typename TrialT::Vector_t fTrialVector;
	
	//! random number generator
	std::tr1::ranlux64_base_01 generator;

public:
	// Constructors
	Crossover() : fCrossoverRate(1), fTrialVector() {generator.seed(time(NULL));}
	Crossover(const double &crossoverrate) : fCrossoverRate(crossoverrate), fTrialVector() {generator.seed(time(NULL));}
	//! Copy constructor
	Crossover(const Crossover& crossover) : fCrossoverRate(crossover.fCrossoverRate), fTrialVector(crossover.fTrialVector) {generator.seed(time(NULL));}
	//! Destructor
	~Crossover() {}
	//! Copy assignment operator
	Crossover& operator=(Crossover crossover);
	
	// Getters
	double GetCrossoverRate() const					{return fCrossoverRate;}
	std::vector<double> GetTrialVector() const		{return fTrialVector;}
	
	// Setters
	void SetCrossoverRate(const double &crossoverrate)	{fCrossoverRate = crossoverrate;}
	
	// other member functions
	void CrossOver(const typename TrialT::Vector_t &targetvector, const typename TrialT::Vector_t &donorvector);
	void Show(std::ostream &s = std::cout);
};


template<typename TrialT>
void Crossover<TrialT>::CrossOver(const typename TrialT::Vector_t &targetvector, const typename TrialT::Vector_t &donorvector)
{
	if(targetvector.size() != donorvector.size())
	{
		std::cout << "Target and donor vector don't have the same size (" << targetvector.size() << " vs. " << donorvector.size() << ")" << std::endl;
		exit(1);
	}
	
	fTrialVector.resize(targetvector.size());
	
	std::tr1::uniform_real<double> unifReal(0., 1.);
	std::tr1::uniform_int<int> unifInt(0, targetvector.size()-1);
	
	for(::vector_size_t i = 0; i < fTrialVector.size(); i++)
	{
		if(unifReal(generator) <= fCrossoverRate || i == static_cast< ::vector_size_t > (unifInt(generator)))
			fTrialVector.at(i) = donorvector.at(i);
		else
			fTrialVector.at(i) = targetvector.at(i);
	}
}

template<typename TrialT>
Crossover<TrialT>& Crossover<TrialT>::operator=(Crossover<TrialT> crossover)
{
	std::swap(fCrossoverRate, crossover.fCrossoverRate);
	std::swap(fTrialVector, crossover.fTrialVector);
	
	generator.seed(time(NULL));
	
	return *this;
}

template<typename TrialT>
void Crossover<TrialT>::Show(std::ostream &s)
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

template <typename TrialT>
std::ostream & operator<<(std::ostream &s, Crossover<TrialT> &crossover)
{
	crossover.Show(s);
	
	return s;
}
#endif // __CROSSOVER__H

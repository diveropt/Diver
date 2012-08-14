#ifndef __MUTATION__H
#define __MUTATION__H

#include <ctime>
#include <algorithm>
#include <tr1/random>
#include <iostream>

#include "global.h"

template <typename TrialT>
class Mutation
{
	int fNDifferenceVectors;				//!< Number of difference vectors. Default = 1
	std::vector<double> fScalingFactors;	//!< Mutation scaling factor somewhat related to the variance between points in each generation. [0, 2]
	typename TrialT::Vector_t fDonorVector;		//!< Mutant donor vector used to create a potential member of the populaton for the next generation.
	int fIBaseVector;						//!< Index of the chosen base vector in the population.
	double fLambda;							//!< Mixing parameter between base and best vector. [0, 1]. Default = 0 (Base instead of best vector)

	std::tr1::mt19937 generator;			//!< random number generator
	
public:
	//! Constructors
	Mutation() : fNDifferenceVectors(1), fScalingFactors(fNDifferenceVectors, 1.), fDonorVector(), fIBaseVector(-1), fLambda(0.) {generator.seed(time(NULL));}
	//! Copy constructor
	Mutation(const Mutation& mutation) : fNDifferenceVectors(mutation.fNDifferenceVectors), fScalingFactors(mutation.fScalingFactors), fDonorVector(mutation.fDonorVector), fIBaseVector(mutation.fIBaseVector), fLambda(mutation.fLambda) {generator.seed(time(NULL));}
	//! Destructor
	~Mutation() {}
	//! Copy assignment operator
	Mutation& operator=(Mutation mutation);
		
	// Getters
	double GetScalingFactor(const ::vector_size_t& index) const		{return fScalingFactors.at(index);}
	std::vector<double> GetScalingFactors() const		{return fScalingFactors;}
	typename TrialT::Vector_t GetDonorVector() const			{return fDonorVector;}
	double GetLambda() const							{return fLambda;}
	
	TrialT* GetBaseVector(std::vector<TrialT*> Population, const int& targetvector);	
	TrialT* GetBestVector(std::vector<TrialT*> Population);
	
	// Setters
	void SetNDifferenceVectors(const int &ndiffvectors)			{fNDifferenceVectors = ndiffvectors; fScalingFactors.resize(fNDifferenceVectors, 1.);}
	void SetScalingFactors(const double &scalingfactor)			{fScalingFactors = std::vector<double>(fNDifferenceVectors, scalingfactor);}
	void SetScalingFactors(std::vector<double> scalingfactors)	{fScalingFactors = scalingfactors;}
	void SetScalingFactor(const ::vector_size_t& index, const double &scalingfactor)	{fScalingFactors.at(index) = scalingfactor;}
	void SetLambda(const double &lambda)						{fLambda = lambda;}
	
	// Other member functions
	void Mutate(std::vector<TrialT*> Population, const ::vector_size_t& targetvector);
	
	void Show(std::ostream &s = std::cout);
};

template<typename TrialT>
Mutation<TrialT>& Mutation<TrialT>::operator=(Mutation<TrialT> mutation)
{
	std::swap(fScalingFactors, mutation.fScalingFactors);
	std::swap(fDonorVector, mutation.fDonorVector);
	std::swap(fIBaseVector, mutation.fIBaseVector);
	
	generator.seed(time(NULL));
	
	return *this;
}

template<typename TrialT>
TrialT* Mutation<TrialT>::GetBaseVector(std::vector<TrialT*> Population, const int& targetvector)
{
	::vector_size_t nP = Population.size();
	
	::vector_size_t pos;
	
	std::tr1::uniform_int<int> unifInt(0, nP-1);
	
	do
		pos = unifInt(generator);
	while (pos == targetvector);
	
	return *(Population.begin() + pos);
}

template<typename TrialT>
TrialT* Mutation<TrialT>::GetBestVector(std::vector<TrialT*> Population)
{
	typename std::vector<TrialT*>::iterator it_best = std::min_element(Population.begin(), Population.end(), CompareFitness<TrialT>);
	
	return *it_best;
}

template<typename TrialT>
void Mutation<TrialT>::Mutate(std::vector<TrialT*> Population, const ::vector_size_t& targetvector)
{
	std::vector<double> basevector = GetBaseVector(Population, targetvector) -> GetPoint();
	std::vector<double> bestvector = GetBestVector(Population) -> GetPoint();
	
	::vector_size_t nP = Population.size();
	fDonorVector.resize(basevector.size());
	
	std::vector< ::vector_size_t> iRandomVector1(fNDifferenceVectors, -1);
	std::vector< ::vector_size_t> iRandomVector2(fNDifferenceVectors, -1);
	
	std::vector<std::vector<double> > RandomVector1(fNDifferenceVectors);
	std::vector<std::vector<double> > RandomVector2(fNDifferenceVectors);
	
	std::tr1::uniform_int<int> unifInt(0, nP-1);

	for(::vector_size_t i = 0; i < (::vector_size_t)fNDifferenceVectors; i++)
	{
		do
		{
			iRandomVector1.at(i) = unifInt(generator);
			iRandomVector2.at(i) = unifInt(generator);
		}
		while(iRandomVector1.at(i) == targetvector || iRandomVector1.at(i) == iRandomVector2.at(i) || fIBaseVector == iRandomVector1.at(i));
		
		RandomVector1.at(i) = Population.at(iRandomVector1.at(i)) -> GetPoint();
		RandomVector2.at(i) = Population.at(iRandomVector2.at(i)) -> GetPoint();
	}
	
	for(::vector_size_t i = 0; i < basevector.size(); i++)
	{
		fDonorVector.at(i) = fLambda*bestvector.at(i) + (1. - fLambda)*basevector.at(i);
		
		for(::vector_size_t j = 0; j < (::vector_size_t)fNDifferenceVectors; j++)
		{
			fDonorVector.at(i) += fScalingFactors.at(j)*(RandomVector1.at(j).at(i) - RandomVector2.at(j).at(i));
		}
	}
}

template <typename TrialT>
void Mutation<TrialT>::Show(std::ostream &s)
{
	s << "=== Mutation ===" << std::endl;
	s << "Number of difference vectors: " << fNDifferenceVectors << std::endl;
	s << "Scaling factors: ";
	if(fScalingFactors.size() == 0)
		s << "empty" << std::endl;
	else
	{
		for(::vector_size_t i = 0; i < fScalingFactors.size(); i++)
			s << fScalingFactors.at(i) << "\t";
		s << std::endl;
	}
	s << "Donor vector: ";
	if(fDonorVector.size() == 0)
		s << "empty" << std::endl;
	else
	{
		for(::vector_size_t i = 0; i < fDonorVector.size(); i++)
			s << fDonorVector.at(i) << "\t";
		s << std::endl;
	}
}

template <typename TrialT>
std::ostream & operator<<(std::ostream &s, Mutation<TrialT> &mutation)
{
	mutation.Show(s);
	
	return s;
}

#endif // __MUTATION__H


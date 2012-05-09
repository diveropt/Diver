#ifndef __MUTATION__H
#define __MUTATION__H

#include <cstdlib>

#include "Trial.h"

class Mutation
{
	int fNDifferenceVectors;				//!< Number of difference vectors. Default = 1
	std::vector<double> fScalingFactors;	//!< Mutation scaling factor somewhat related to the variance between points in each generation. [0, 2]
	std::vector<double> fDonorVector;		//!< Mutant donor vector used to create a potential member of the populaton for the next generation.
	int fIBaseVector;						//!< Index of the chosen base vector in the population.
	int fIBestVector;						//!< Index of the best vector in the population.
	double fLambda;							//!< Mixing parameter between base and best vector. [0, 1]. Default = 0 (Base instead of best vector)
	
public:
	//! Constructors
	Mutation() : fNDifferenceVectors(1), fScalingFactors(fNDifferenceVectors, 1.), fDonorVector(), fIBaseVector(-1), fIBestVector(-1), fLambda(0.) {std::srand(time(0));}
	//! Copy constructor
	Mutation(const Mutation& mutation) : fNDifferenceVectors(mutation.fNDifferenceVectors), fScalingFactors(mutation.fScalingFactors), fDonorVector(mutation.fDonorVector), fIBaseVector(mutation.fIBaseVector), fIBestVector(mutation.fIBestVector), fLambda(mutation.fLambda) {std::srand(time(0));}
	//! Destructor
	~Mutation() {}
	//! Copy assignment operator
	Mutation& operator=(Mutation mutation);
	
	// Getters
	double GetScalingFactor(const ::vector_size_t& index) const		{return fScalingFactors.at(index);}
	std::vector<double> GetScalingFactors() const		{return fScalingFactors;}
	std::vector<double> GetDonorVector() const			{return fDonorVector;}
	double GetLambda() const							{return fLambda;}
	Trial* GetBaseVector(std::vector<Trial*> Population, const int& targetvector);
	Trial* GetBestVector(std::vector<Trial*> Population);
	
	// Setters
	void SetNDifferenceVectors(const int &ndiffvectors)			{fNDifferenceVectors = ndiffvectors; fScalingFactors.resize(fNDifferenceVectors, 1.);}
	void SetScalingFactors(const double &scalingfactor)			{fScalingFactors = std::vector<double>(fNDifferenceVectors, scalingfactor);}
	void SetScalingFactors(std::vector<double> scalingfactors)	{fScalingFactors = scalingfactors;}
	void SetScalingFactor(const ::vector_size_t& index, const double &scalingfactor)	{fScalingFactors.at(index) = scalingfactor;}
	void SetLambda(const double &lambda)						{fLambda = lambda;}
	
	// Other member functions
	void Mutate(std::vector<Trial*> Population, const ::vector_size_t& targetvector);
	void Show(std::ostream &s = std::cout);
};

std::ostream & operator<<(std::ostream &s, Mutation &mutation);

#endif // __MUTATION__H


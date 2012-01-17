#ifndef __MUTATION__H
#define __MUTATION__H

#include <iostream>
#include <vector>
#include "Trial.h"

class Mutation
{
	double fScalingFactor;				//!< Mutation scaling factor somewhat related to the variance between points in each generation. [0, 2]
	std::vector<double> fDonorVector;	//!< Mutant donor vector used to create a potential member of the populaton for the next generation.
	int fIBaseVector;					//!< Index of the chosen base vector in the population.
	
public:
	//! Constructors
	Mutation() : fScalingFactor(1.), fDonorVector(std::vector<double>()), fIBaseVector(-1) {srand(time(0));}
	explicit Mutation(const double& scalingfactor) : fScalingFactor(scalingfactor), fDonorVector(std::vector<double>()), fIBaseVector(-1) {srand(time(0));}
	//! Copy constructor
	Mutation(const Mutation& mutation) : fScalingFactor(mutation.fScalingFactor), fDonorVector(mutation.fScalingFactor), fIBaseVector(mutation.fIBaseVector) {srand(time(0));}
	//! Destructor
	~Mutation() {fDonorVector.clear();}
	//! Copy assignment operator
	Mutation& operator=(Mutation mutation);
	
	// Getters
	double GetScalingFactor() const				{return fScalingFactor;}
	std::vector<double> GetDonorVector() const	{return fDonorVector;}
	Trial* GetBaseVector(std::vector<Trial*> Population, const int& targetvector);
	
	// Setters
	void SetScalingFactor(const double &scalingfactor)		{fScalingFactor = scalingfactor;}
	
	// Other member functions
	void Mutate(std::vector<Trial*> Population, const int& targetvector);
	void Show();
};

#endif // __MUTATION__H


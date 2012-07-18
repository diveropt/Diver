#ifndef __PARAMETERSPACE__H
#define __PARAMETERSPACE__H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#include "global.h"

class ParameterSpace
{
	std::vector<std::string> fParameters;
	std::vector<double> fLowerBounds;
	std::vector<double> fUpperBounds;
	
	void CheckBounds();
	
public:	
	// Constructors
	ParameterSpace() : fParameters(), fLowerBounds(), fUpperBounds() {}
	ParameterSpace(std::vector<std::string> parameters, std::vector<double> lowerBounds, std::vector<double> upperBounds) : fParameters(parameters), fLowerBounds(lowerBounds), fUpperBounds(upperBounds) {CheckBounds();}
	//! Copy Constructor
	ParameterSpace(const ParameterSpace& parspace) : fParameters(parspace.fParameters), fLowerBounds(parspace.fLowerBounds), fUpperBounds(parspace.fUpperBounds) {}
	//! Destructor
	~ParameterSpace() {}
	//! Copy assignment operator
	ParameterSpace& operator=(ParameterSpace parspace);
	
	// Getters
	::vector_size_t GetNParameters() const								{return fParameters.size();}
	std::string GetParameterName(const ::vector_size_t& index)	const	{return fParameters.at(index);}
	double GetParameterLowerBound(const ::vector_size_t& index) const	{return fLowerBounds.at(index);}
	double GetParameterUpperBound(const ::vector_size_t& index) const	{return fUpperBounds.at(index);}
	std::vector<std::string> GetParameterNames() const					{return fParameters;}
	std::vector<double> GetParameterLowerBounds() const					{return fLowerBounds;}
	std::vector<double> GetParameterUpperBounds() const					{return fUpperBounds;}
		
	// other member functions
	void AddParameter(std::string parameter, double lowerBound, double upperBound);
	void Show(std::ostream &s = std::cout);
};

std::ostream & operator<<(std::ostream &s, ParameterSpace &parspace);

#endif // __PARAMETERSPACE__H

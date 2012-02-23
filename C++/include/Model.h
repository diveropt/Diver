/*
 *  Model.h
 *  DE
 *
 *  Created by Antje Putze on 06/11/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __MODEL__H
#define __MODEL__H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#include "global.h"

class Model
{
	std::vector<std::string> fParameters;
	std::vector<double> fLowerBounds;
	std::vector<double> fUpperBounds;
	
	void CheckBounds();
	
public:	
	// Constructors
	Model() : fParameters(), fLowerBounds(), fUpperBounds() {}
	Model(std::vector<std::string> parameters, std::vector<double> lowerBounds, std::vector<double> upperBounds) : fParameters(parameters), fLowerBounds(lowerBounds), fUpperBounds(upperBounds) {CheckBounds();}
	//! Copy Constructor
	Model(const Model& model) : fParameters(model.fParameters), fLowerBounds(model.fLowerBounds), fUpperBounds(model.fUpperBounds) {}
	//! Destructor
	~Model() {}
	//! Copy assignment operator
	Model& operator=(Model model);
	
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

std::ostream & operator<<(std::ostream &s, Model &model);

#endif // __MODEL__H

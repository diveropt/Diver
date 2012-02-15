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

#include "cstdlib"
#include "iostream"
#include "vector"
#include "string"

class Model
{
	std::vector<std::string> fParameters;
	std::vector<double> fUpperBounds;
	std::vector<double> fLowerBounds;
	
	void CheckBounds();
	
public:
	// Constructors
	Model() : fParameters(std::vector<std::string>()), fUpperBounds(std::vector<double>()), fLowerBounds(std::vector<double>()) {}
	Model(std::vector<std::string> parameters, std::vector<double> lowerBounds, std::vector<double> upperBounds) : fParameters(parameters), fUpperBounds(lowerBounds), fLowerBounds(upperBounds) {CheckBounds();}
	//! Copy Constructor
	Model(const Model& model) : fParameters(model.fParameters), fLowerBounds(model.fLowerBounds), fUpperBounds(model.fUpperBounds) {}
	//! Destructor
	~Model();
	//! Copy assignment operator
	Model& operator=(Model model);
	
	// Getters
	unsigned int GetNParameters() const						{return fParameters.size();}
	std::string GetParameterName(int index)	const			{return fParameters.at(index);}
	double GetParameterLowerBound(int index) const			{return fLowerBounds.at(index);}
	double GetParameterUpperBound(int index) const			{return fUpperBounds.at(index);}
	
	// Setters
	void SetNParameters(const int& size)								{fParameters.resize(size); fLowerBounds.resize(size); fUpperBounds.resize(size);}
	void SetParameterName(const int& index, const std::string& name)	{fParameters.at(index) = name;}
	void SetParameterLowerBound(const int& index, double lowerBound)	{fLowerBounds.at(index) = lowerBound;}
	void SetParameterUpperBound(const int& index, double upperBound)	{fUpperBounds.at(index) = upperBound;}
	
	// other member functions
	void AddParameter(std::string parameter, double lowerBound, double upperBound);
	void Show();
};

#endif // __MODEL__H

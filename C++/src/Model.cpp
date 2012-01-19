/*
 *  Model.cpp
 *  DE
 *
 *  Created by Antje Putze on 06/11/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "Model.h"

//! Destructor
Model::~Model()
{
	fParameters.clear();
	fLowerBounds.clear();
	fUpperBounds.clear();
}

//! Copy assignment operator
Model& Model::operator=(Model model)
{
	std::swap(fParameters, model.fParameters);
	std::swap(fLowerBounds, model.fLowerBounds);
	std::swap(fUpperBounds, model.fUpperBounds);
}
	
void Model::CheckBounds()
{
	if(fParameters.size() != fLowerBounds.size() && fParameters.size() != fUpperBounds.size())
	{
		std::cout << "The numbers of parameter names and ranges are not the same. Please check!" << std::endl;
		std::cout << "fParameters.size:\t" << fParameters.size() << std::endl;
		std::cout << "fLowerBounds.size:\t" << fLowerBounds.size() << std::endl;
		std::cout << "fUpperBounds.size:\t" << fUpperBounds.size() << std::endl;
		exit(1);
	}
	
	for(unsigned int i = 0; i < fParameters.size(); i++)
	{
		if(fLowerBounds.at(i) > fUpperBounds.at(i))
		{
			std::cout << "The lower bound of parameter " << fParameters.at(i) << " is greater than its upper bound " << fLowerBounds.at(i) << " > " << fUpperBounds.at(i) << std::endl;
			std::cout << "Swapping values" << std::endl;
			std::swap(fLowerBounds.at(i), fUpperBounds.at(i));
		}
	}
}
	
void Model::AddParameter(std::string parameter, double lowerBound, double upperBound)
{
	fParameters.push_back(parameter);
	fLowerBounds.push_back(lowerBound);
	fUpperBounds.push_back(upperBound);
	
	CheckBounds();
}

void Model::Show()
{
	std::cout << "=== Model ===" << std::endl;
	std::cout << "Number of parameters: " << fParameters.size() << std::endl;
	for(unsigned int i = 0; i < fParameters.size(); i++)
		std::cout << "Parameter " << fParameters.at(i) << " is defined between " << fLowerBounds.at(i) << " and " << fUpperBounds.at(i) << std::endl;
}


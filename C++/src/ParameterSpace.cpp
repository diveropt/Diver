#include "ParameterSpace.h"

//! Copy assignment operator
ParameterSpace& ParameterSpace::operator=(ParameterSpace parspace)
{
	std::swap(fParameters, parspace.fParameters);
	std::swap(fLowerBounds, parspace.fLowerBounds);
	std::swap(fUpperBounds, parspace.fUpperBounds);
	
	return *this;
}
	
void ParameterSpace::CheckBounds()
{
	// Check wether the vectors have the same number of entries. If not exit.
	if(fParameters.size() != fLowerBounds.size() && fParameters.size() != fUpperBounds.size())
	{
		std::cout << "The numbers of parameter names and ranges are not the same. Please check!" << std::endl;
		std::cout << "fParameters.size:\t" << fParameters.size() << std::endl;
		std::cout << "fLowerBounds.size:\t" << fLowerBounds.size() << std::endl;
		std::cout << "fUpperBounds.size:\t" << fUpperBounds.size() << std::endl;
		std::exit(1);
	}
	
	for(::vector_size_t i = 0; i < fParameters.size(); i++)
	{
		if(fLowerBounds.at(i) > fUpperBounds.at(i))
		{
			std::cout << "The lower bound of parameter " << fParameters.at(i) << " is greater than its upper bound " << fLowerBounds.at(i) << " > " << fUpperBounds.at(i) << std::endl;
			std::cout << "Swapping values" << std::endl;
			std::swap(fLowerBounds.at(i), fUpperBounds.at(i));
		}
	}
}
	
void ParameterSpace::AddParameter(std::string parameter, double lowerBound, double upperBound)
{
	// Check wether lowerBound < upperBound. If not, swap values.
	if(lowerBound > upperBound)
		std::swap(lowerBound, upperBound);
	
	fParameters.push_back(parameter);
	fLowerBounds.push_back(lowerBound);
	fUpperBounds.push_back(upperBound);
}

void ParameterSpace::Show(std::ostream &s)
{
	s << "=== ParameterSpace ===" << std::endl;
	s << "Number of parameters: " << fParameters.size() << std::endl;
	for(::vector_size_t i = 0; i < fParameters.size(); i++)
		s << "Parameter " << fParameters.at(i) << " is defined between " << fLowerBounds.at(i) << " and " << fUpperBounds.at(i) << std::endl;
}

std::ostream & operator<<(std::ostream &s, ParameterSpace &parspace)
{
	parspace.Show(s);
	
	return s;
}
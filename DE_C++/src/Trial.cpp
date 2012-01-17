#include "Trial.h"
#include <iomanip>

//! Copy assignment operator
Trial& Trial::operator=(Trial trial)
{
	std::swap(fPoint, trial.fPoint);
	std::swap(fFitness, trial.fFitness);
}

//! Shows trial's characteristics
void Trial::Show()
{
	std::cout << "Point: (";
	for(unsigned int i = 0; i < fPoint.size(); i++)
	{
		std::cout << std::setw(8) << fPoint.at(i);
		if(i < fPoint.size() - 1)
			std::cout << ", ";
		else
			std::cout << ")\t";
	}
	std::cout << "Fitness: " << std::setw(8) << fFitness << std::endl;
}
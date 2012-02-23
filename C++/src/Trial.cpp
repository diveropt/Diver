#include "Trial.h"
#include <iomanip>

//! Copy assignment operator
Trial& Trial::operator=(Trial trial)
{
	std::swap(fPoint, trial.fPoint);
	std::swap(fFitness, trial.fFitness);
	
	return *this;
}

//! Shows trial's characteristics
void Trial::Show(std::ostream &s)
{
	std::cout << "Point: (";
	for(::vector_size_t i = 0; i < fPoint.size(); i++)
	{
		s << std::setw(8) << fPoint.at(i);
		if(i < fPoint.size() - 1)
			s << ", ";
		else
			s << ")\t";
	}
	s << "Fitness: " << std::setw(8) << fFitness << std::endl;
}

std::ostream & operator<<(std::ostream &s, Trial &trial)
{
	trial.Show(s);
	
	return s;
}
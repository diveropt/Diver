#include <iomanip>
#include "Trial.h"

//! Static member initialisation
std::vector<std::string> Trial::fParameters;
Trial::Vector_t Trial::fLowerBounds;
Trial::Vector_t Trial::fUpperBounds;
std::tr1::ranlux64_base_01 Trial::generator;

//! Copy assignment operator
Trial& Trial::operator=(Trial trial)
{
	std::swap(fPoint, trial.fPoint);
	std::swap(fFitness, trial.fFitness);
	std::swap(fParameters, trial.fParameters);
	std::swap(fLowerBounds, trial.fLowerBounds);
	std::swap(fUpperBounds, trial.fUpperBounds);
	
	return *this;
}

void Trial::Init()
{
	// Set seed if not done already
	static bool IsSeed = false;
	if(!IsSeed)
	{
		generator.seed(time(NULL));
		IsSeed = true;
	}
	
	// Initialise each point with a random value drawn uniformly between its lower and upper bound
	for(::vector_size_t i = 0; i < fPoint.size(); i++)
	{
		std::tr1::uniform_real<double> unifReal(fLowerBounds.at(i), fUpperBounds.at(i));
		fPoint.at(i) = unifReal(generator);
	}
}

//! Add a new parameter with its name and validity range to the parameter space.
void Trial::AddParameter(std::string parameter, Vector_type_t lowerBound, Vector_type_t upperBound)
{
	// Check wether lowerBound < upperBound. If not, swap values.
	if(lowerBound > upperBound)
		std::swap(lowerBound, upperBound);
	
	// Add parameter name and bounds
	fParameters.push_back(parameter);
	fLowerBounds.push_back(lowerBound);
	fUpperBounds.push_back(upperBound);
	
	// Resize point vector
	fPoint.resize(fParameters.size());
}

//! Show the parameter space's characteristics
void Trial::ShowParameters(std::ostream &s)
{
	s << "=== Parameter Space ===" << std::endl;
	s << "Number of parameters: " << fParameters.size() << std::endl;
	for(::vector_size_t i = 0; i < fParameters.size(); i++)
		s << "Parameter " << fParameters.at(i) << " has the value " << fPoint.at(i) << " and is defined between " << fLowerBounds.at(i) << " and " << fUpperBounds.at(i) << std::endl;
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

//! Define the << operator for the Trial class
std::ostream & operator<<(std::ostream &s, Trial &trial)
{
	trial.Show(s);
	
	return s;
}
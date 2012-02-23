#ifndef __CROSSOVER__H
#define __CROSSOVER__H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

class Crossover
{
	double fCrossoverRate;
	std::vector<double> fTrialVector;
	
public:
	// Constructors
	Crossover() : fCrossoverRate(1), fTrialVector() {std::srand(time(0));}
	Crossover(const double &crossoverrate) : fCrossoverRate(crossoverrate), fTrialVector() {std::srand(time(0));}
	//! Copy constructor
	Crossover(const Crossover& crossover) : fCrossoverRate(crossover.fCrossoverRate), fTrialVector(crossover.fTrialVector) {std::srand(time(0));}
	//! Destructor
	~Crossover() {}
	//! Copy assignment operator
	Crossover& operator=(Crossover crossover);
	
	// Getters
	double GetCrossoverRate() const					{return fCrossoverRate;}
	std::vector<double> GetTrialVector() const		{return fTrialVector;}
	
	// Setters
	void SetCrossoverRate(const double &crossoverrate)	{fCrossoverRate = crossoverrate;}
	
	// other member functions
	void CrossOver(const std::vector<double> &targetvector, const std::vector<double> &donorvector);
	void Show(std::ostream &s = std::cout);
};

std::ostream & operator<<(std::ostream &s, Crossover &crossover);

#endif // __CROSSOVER__H

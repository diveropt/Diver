#ifndef __CROSSOVER__H
#define __CROSSOVER__H

#include <iostream>
#include <vector>

class Crossover
{
	double fCrossoverRate;
	std::vector<double> fTrialVector;
	
public:
	// Constructors
	Crossover() : fCrossoverRate(1), fTrialVector(std::vector<double>()) {srand(time(0));}
	explicit Crossover(const double &crossoverrate) : fCrossoverRate(crossoverrate), fTrialVector(std::vector<double>()) {srand(time(0));}
	//! Copy constructor
	Crossover(const Crossover& crossover) : fCrossoverRate(crossover.fCrossoverRate), fTrialVector(crossover.fTrialVector) {srand(time(0));}
	//! Destructor
	~Crossover() {fTrialVector.clear();}
	//! Copy assignment operator
	Crossover& operator=(Crossover crossover);
	
	// Getters
	double GetCrossoverRate() const					{return fCrossoverRate;}
	std::vector<double> GetTrialVector() const		{return fTrialVector;}
	
	// Setters
	void SetCrossoverRate(const double &crossoverrate)	{fCrossoverRate = crossoverrate;}
	
	// other member functions
	void CrossOver(const std::vector<double> &targetvector, const std::vector<double> &donorvector);
	void Show();
};

#endif // __CROSSOVER__H
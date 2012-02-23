#ifndef __TRIAL__H
#define __TRIAL__H

#include <iostream>
#include <vector>

#include "global.h"

class Trial
{
	std::vector<double> fPoint;		//!< Parameter space point
	double fFitness;				//!< Fitness of the parameter space point
	
public:
	// Constructors
	Trial() : fPoint(), fFitness(0.) {}
	explicit Trial(const std::vector<double> &point) : fPoint(point), fFitness(0.) {}
	Trial(const std::vector<double> &point, const double &fitness) : fPoint(point), fFitness(fitness) {}
	//! Copy constructors
	Trial(const Trial& trial) : fPoint(trial.fPoint), fFitness(trial.fFitness) {};
	//! Destructor
	~Trial() {}
	//! Copy assignment operator
	Trial& operator=(Trial trial);
	
	// Getters
	::vector_size_t GetN() const						{return fPoint.size();}
	std::vector<double> GetPoint() const				{return fPoint;}
	double GetFitness()	const							{return fFitness;}
	double GetValue(const ::vector_size_t &index) const {return fPoint.at(index);}
	
	// Setters
	void SetN(const ::vector_size_t &n)								{fPoint.resize(n);}
	void SetPoint(const std::vector<double> &point)							{fPoint = point;}
	void SetFitness(const double &fitness)									{fFitness = fitness;}
	void SetValue(const ::vector_size_t &index, const double &value)	{fPoint.at(index) = value;}
	
	// Other member functions
	void Show(std::ostream &s = std::cout);
};

std::ostream & operator<<(std::ostream &s, Trial &trial);
#endif // __TRIAL__H
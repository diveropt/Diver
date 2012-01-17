#ifndef __TRIAL__H
#define __TRIAL__H

#include <iostream>
#include <vector>

class Trial
{
	std::vector<double> fPoint;		//!< Parameter space point
	double fFitness;				//!< Fitness of the parameter space point
	
public:
	// Constructors
	Trial() : fPoint(std::vector<double>()), fFitness(0.) {}
	explicit Trial(const std::vector<double> &point) : fPoint(point), fFitness(0.) {}
	Trial(const std::vector<double> &point, const double &fitness) : fPoint(point), fFitness(fitness) {}
	//! Copy constructors
	Trial(const Trial& trial) : fPoint(trial.fPoint), fFitness(trial.fFitness) {};
	//! Destructor
	~Trial() {fPoint.clear();}
	//! Copy assignment operator
	Trial& operator=(Trial trial);
	
	// Getters
	unsigned int GetN() const				{return fPoint.size();}
	std::vector<double> GetPoint() const	{return fPoint;}
	double GetFitness()	const				{return fFitness;}
	double GetValue(const int &index) const {return fPoint.at(index);}
	
	// Setters
	void SetN(const unsigned int &n)					{fPoint.resize(n);}
	void SetPoint(const std::vector<double> &point)		{fPoint = point;}
	void SetFitness(const double &fitness)				{fFitness = fitness;}
	void SetValue(const int &index, const double &value){fPoint.at(index) = value;}
	
	// Other member functions
	void Show();
};

#endif // __TRIAL__H
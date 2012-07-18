#ifndef __TRIAL__H
#define __TRIAL__H

#include <iostream>
#include <vector>

#include <tr1/random>
#include "global.h"

class Trial
{
	std::vector<double> fPoint;		//!< Parameter space point
	double fFitness;				//!< Fitness of the parameter space point
	
	// Static member variables
	static std::vector<std::string> fParameters;	//! Parameter names
	static std::vector<double> fLowerBounds;		//! Parameter lower bounds
	static std::vector<double> fUpperBounds;		//! Parameter upper bounds
	
public:
	// Constructors
	Trial() : fPoint(fParameters.size()), fFitness(0.) {}
	explicit Trial(const std::vector<double> &point) : fPoint(point), fFitness(0.) {}
	Trial(const std::vector<double> &point, const double &fitness) : fPoint(point), fFitness(fitness) {}
	//! Copy constructors
	Trial(const Trial& trial) : fPoint(trial.fPoint), fFitness(trial.fFitness) {}
	//! Destructor
	~Trial() {}
	//! Copy assignment operator
	Trial& operator=(Trial trial);
	
	// Getters
	::vector_size_t GetN() const						{return fPoint.size();}
	std::vector<double> GetPoint() const				{return fPoint;}
	double GetFitness()	const							{return fFitness;}
	double GetValue(const ::vector_size_t &index) const {return fPoint.at(index);}
	std::string GetParameterName(const ::vector_size_t& index)	const	{return fParameters.at(index);}
	double GetParameterLowerBound(const ::vector_size_t& index) const	{return fLowerBounds.at(index);}
	double GetParameterUpperBound(const ::vector_size_t& index) const	{return fUpperBounds.at(index);}
	std::vector<std::string> GetParameterNames() const					{return fParameters;}
	std::vector<double> GetParameterLowerBounds() const					{return fLowerBounds;}
	std::vector<double> GetParameterUpperBounds() const					{return fUpperBounds;}
	
	// Setters
	void SetN(const ::vector_size_t &n)									{fPoint.resize(n);}
	void SetPoint(const std::vector<double> &point)						{fPoint = point;}
	void SetFitness(const double &fitness)								{fFitness = fitness;}
	void SetValue(const ::vector_size_t &index, const double &value)	{fPoint.at(index) = value;}
	
	// Other member functions
	void AddParameter(std::string parameter, double lowerBound, double upperBound);
	void ShowParameters(std::ostream &s = std::cout);
	void Show(std::ostream &s = std::cout);
};

std::ostream & operator<<(std::ostream &s, Trial &trial);
#endif // __TRIAL__H
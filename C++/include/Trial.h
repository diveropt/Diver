#ifndef __TRIAL__H
#define __TRIAL__H

#include <iostream>
#include <vector>
#include <tr1/random>

#include "global.h"

class Trial
{
public:
	typedef double Fitness_type_t;
	typedef double Vector_type_t;
	typedef std::vector<Vector_type_t> Vector_t;

private:
	Vector_t fPoint;				//!< Parameter space point
	Fitness_type_t fFitness;		//!< Fitness of the parameter space point
	
	// Static member variables
	static std::vector<std::string> fParameters;	//! Parameter names
	static Vector_t fLowerBounds;		//! Parameter lower bounds
	static Vector_t fUpperBounds;		//! Parameter upper bounds
	
	//! Random number generator
	static std::tr1::ranlux64_base_01 generator;
	
public:
	// Constructors
	Trial() : fPoint(fParameters.size()), fFitness(0.) {}
	explicit Trial(const Vector_t &point) : fPoint(point), fFitness(0.) {}
	Trial(const Vector_t &point, const Fitness_type_t &fitness) : fPoint(point), fFitness(fitness) {}
	//! Copy constructors
	Trial(const Trial& trial) : fPoint(trial.fPoint), fFitness(trial.fFitness) {}
	//! Destructor
	~Trial() {}
	//! Copy assignment operator
	Trial& operator=(Trial trial);
	//! operator
	
	// Getters
	::vector_size_t GetN() const												{return fPoint.size();}
	Vector_t GetPoint() const													{return fPoint;}
	Fitness_type_t GetFitness()	const											{return fFitness;}
	Vector_type_t GetValue(const ::vector_size_t &index) const					{return fPoint.at(index);}
	std::string GetParameterName(const ::vector_size_t& index)	const			{return fParameters.at(index);}
	Vector_type_t GetParameterLowerBound(const ::vector_size_t& index) const	{return fLowerBounds.at(index);}
	Vector_type_t GetParameterUpperBound(const ::vector_size_t& index) const	{return fUpperBounds.at(index);}
	std::vector<std::string> GetParameterNames() const							{return fParameters;}
	Vector_t GetParameterLowerBounds() const									{return fLowerBounds;}
	Vector_t GetParameterUpperBounds() const									{return fUpperBounds;}
	
	// Setters
	void SetN(const ::vector_size_t &n)										{fPoint.resize(n);}
	void SetPoint(const Vector_t &point)									{fPoint = point;}
	void SetFitness(const Fitness_type_t &fitness)							{fFitness = fitness;}
	void SetValue(const ::vector_size_t &index, const Vector_type_t &value)	{fPoint.at(index) = value;}
	
	// Other member functions
	void Init(); //! Initialises the parameter space point
	void AddParameter(std::string parameter, Vector_type_t lowerBound, Vector_type_t upperBound);
	void ShowParameters(std::ostream &s = std::cout);
	void Show(std::ostream &s = std::cout);
};

std::ostream & operator<<(std::ostream &s, Trial &trial);
#endif // __TRIAL__H
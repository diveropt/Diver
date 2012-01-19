#ifndef __TESTFUNCTIONS__H
#define __TESTFUNCTIONS__H

#include "vector"
#include "math.h"

class TwoDGaussian
{
public:
	double operator() (std::vector<double> parameters) const
	{
		double r  = 0.;
		double f1 = 0.;
		
		bool flag = false;
		
		for(unsigned int i = 0; i < parameters.size(); i++)
		{
			if(parameters.at(i) > 0. && parameters.at(i) < 1.)
				r += pow(parameters.at(i) - 0.5, 2.);
			else
				flag = true;
		}
		
		f1 = exp(r/-0.15);
		
		if(flag)
			return 1.e99;
		else
			return 1./f1;
	}
};

class StairCase
{
public:
	double operator() (std::vector<double> parameters) const
	{
		double r  = 0.;
		double f2 = 0.;
		
		bool flag = false;
		
		for(unsigned int i = 0; i < parameters.size(); i++)
		{
			if(parameters.at(i) > 0. && parameters.at(i) < 1.)
				r += int(10*parameters.at(i) - 1.e-6)/9.;
			else
				flag = true;
		}
		
		f2 = r/(double)parameters.size();
		
		if(flag)
			return 1.e99;
		else
			return 1./f2;
	}
};

class CocentricRings
{
public:
	double operator() (std::vector<double> parameters) const
	{
		double r = 0.;
		double f3 = 0.;
		
		bool flag = false;
		
		for(unsigned int i = 0; i < parameters.size(); i++)
		{
			if(parameters.at(i) > 0. && parameters.at(i) < 1.)
				r += pow(parameters.at(i) - 0.5, 2.);
			else
				flag = true;
		}
		
		f3 = pow(cos(9*M_PI*sqrt(r)), 2.)*exp(-r/0.15);
		
		if(flag)
			return 1.e99;
		else
			return 1./f3;
	}
};

class BroadNarrowGaussians
{
public:
	double operator() (std::vector<double> parameters) const
	{
		double r1 = 0.;
		double r2 = 0.;
		double f4 = 0.;
		
		bool flag = false;
		
		for(unsigned int i = 0; i < parameters.size(); i++)
		{
			if(parameters.at(i) > 0. && parameters.at(i) < 1.)
			{
				r1 += pow(parameters.at(i) - 0.5, 2.);
				r2 += pow(parameters.at(i) - 0.2, 2.);
			}
			else
				flag = true;
		}
		
		double A1 = 0.7;
		double A2 = 1 - 0.7*exp(-r1/0.15);
		
		f4 = A1*exp(-r1/0.15) + A2*exp(-r2/0.005);
		if(flag)
			return 1.e99;
		else
			return 1./f4;
	}
};

class GaussianHill
{
public:
	double operator() (std::vector<double> parameters) const
	{
		double r = 1.;
		double f5 = 0.;
		const int n = 9.;
		
		bool flag = false;
		
		for(unsigned int i = 0; i < parameters.size(); i++)
		{
			if(parameters.at(i) > 0. && parameters.at(i) < 1.)
				r *= pow(parameters.at(i)*(1 - parameters.at(i))*sin(n*M_PI*parameters.at(i)), 2.);
			else
				flag = true;
		}
		
		f5 = 16*16*r;
		
		if(flag)
			return 1.e99;
		else
			return 1./f5;
	}
};

#endif // __TESTFUNCTIONS__H
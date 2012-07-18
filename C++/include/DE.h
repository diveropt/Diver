#ifndef __DE__H
#define __DE__H

#include <vector>
#include <iostream>
#include "global.h"

template <typename TrialT, typename MutationSchemeT, typename CrossoverSchemeT, typename SelectionSchemeT>
class DE
{
	// Pointers to template parameters
	TrialT *fTrial;
	MutationSchemeT *fMutation;
	CrossoverSchemeT *fCrossover;
	SelectionSchemeT *fSelection;	
	
	//! number of individuals in the population
	const int fnInd;
	//! number of generations
	const int fnGen;
	
	// other member functions
	template <typename Functor>
	TrialT* FirstTrial(Functor const &function);
	template <typename Functor>
	TrialT* NextTrial(std::vector<TrialT*> Population, const int& targetvector, Functor const &function);
	
public:
	//! Constructors
	DE(const int& nInd = 20, const int& nGen = 100);
	//! Copy constructor
	DE(const DE &de);
	//! Destructor
	~DE();
	//! Copy assignment operator
	DE& operator=(DE de);
	
	// Getter
	TrialT* GetTrial() const				{return fTrial;}
	MutationSchemeT* GetMutation() const	{return fMutation;}
	CrossoverSchemeT* GetCrossover() const	{return fCrossover;}
	SelectionSchemeT* GetSelection() const	{return fSelection;}
	
	const int GetNIndividuals() const	{return fnInd;}
	const int GetNGenerations() const	{return fnGen;}
		
	//! Run the DE
	template <typename Functor>
	void Run(Functor const &function);	
	
	//! Show the DE configuration
	void Show(std::ostream &s = std::cout);
};

//! Constructor
template <typename TrialT, typename MutationSchemeT, typename CrossoverSchemeT, typename SelectionSchemeT>
DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT>::DE(const int& nInd, const int& nGen) : fTrial(new TrialT), fMutation(new MutationSchemeT), fCrossover(new CrossoverSchemeT), fSelection(new SelectionSchemeT), fnInd(nInd), fnGen(nGen)
{
	srand(time(0));	
}

//! Copy constructor
template <typename TrialT, typename MutationSchemeT, typename CrossoverSchemeT, typename SelectionSchemeT>
DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT>::DE(const DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT> &de)
{
	if(fTrial)
		delete fTrial;
	fTrial = new TrialT(*de.fTrial);
	
	if(fMutation)
		delete fMutation;
	fMutation = new MutationSchemeT(*de.fMutation);
	
	if(fCrossover)
		delete fCrossover;
	fCrossover = new CrossoverSchemeT(*de.fCrossover);
	
	if(fSelection)
		delete fSelection;
	fSelection = new SelectionSchemeT(*de.fSelection);
	
	fnInd = de.fnInd;
	fnGen = de.fnGen;
	
	srand(time(0));	
}	

//! Destructor
template <typename TrialT, typename MutationSchemeT, typename CrossoverSchemeT, typename SelectionSchemeT>
DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT>::~DE()
{
	delete fTrial;
	fTrial = NULL;
	
	delete fMutation;
	fMutation = NULL;
	
	delete fCrossover;
	fCrossover = NULL;
	
	delete fSelection;
	fSelection = NULL;	
}

//! Copy assignment operator
template <typename TrialT, typename MutationSchemeT, typename CrossoverSchemeT, typename SelectionSchemeT>
DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT>& DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT>::operator=(DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT> de)
{
	if(fTrial)
		delete fTrial;
	fTrial = new TrialT(de.fTrial);
	
	if(fMutation)
		delete fMutation;
	fMutation = new MutationSchemeT(de.fMutation);
	
	if(fCrossover)
		delete fCrossover;
	fCrossover = new CrossoverSchemeT(de.fCrossover);

	if(fSelection)
		delete fSelection;
	fSelection = new SelectionSchemeT(de.fSelection);

	srand(time(0));
	
	return *this;
}

template <typename TrialT, typename MutationSchemeT, typename CrossoverSchemeT, typename SelectionSchemeT>
template <typename Functor>
TrialT* DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT>::FirstTrial(Functor const &function)
{
	TrialT *first = new TrialT();

	double randomnumber = 0;
	for(::vector_size_t i = 0; i < first -> GetN(); i++)
	{
		randomnumber = first -> GetParameterLowerBound(i) + (static_cast<double>(rand())/(RAND_MAX + 1.0))*(first -> GetParameterUpperBound(i) - first -> GetParameterLowerBound(i));
		first -> SetValue(i, randomnumber);
	}

	first -> SetFitness(function(first -> GetPoint()));

	return first;
}

template <typename TrialT, typename MutationSchemeT, typename CrossoverSchemeT, typename SelectionSchemeT>
template <typename Functor>
TrialT* DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT>::NextTrial(std::vector<TrialT*> Population, const int& targetvector, Functor const &function)
{
	// Mutate 
	fMutation -> Mutate(Population, targetvector);

	// Cross over
	fCrossover -> CrossOver(Population.at(targetvector) -> GetPoint(), fMutation -> GetDonorVector());
	
	TrialT *trial = new TrialT(fCrossover -> GetTrialVector());
	trial -> SetFitness(function(trial -> GetPoint()));
	
	// Select
	TrialT *result = fSelection -> Select(Population.at(targetvector), trial);
	
	// Free memory
	if(result != trial)
	{
		delete trial;
		trial = NULL;
	}
	
	return result;
}

template <typename TrialT, typename MutationSchemeT, typename CrossoverSchemeT, typename SelectionSchemeT>
template <typename Functor>
void DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT>::Run(Functor const &function)
{
	std::vector<std::vector<TrialT*> > triallist(fnGen + 1);
	
	std::cout << " ==> First initialisation" << std::endl;
	for(int i = 0; i < fnInd; i++)
	{
		triallist.at(0).push_back(FirstTrial(function));
	}
	
	std::cout << " ** The best individual ** \n" << *(GetSelection() -> SelectBest(triallist.at(0))) << std::endl;
	
	for(int i = 0; i < fnGen; i++)
	{
		std::cout << "+++ Generation " << i << std::endl;
		
		for(int j = 0; j < fnInd; j++)
		{
			triallist.at(i+1).push_back(NextTrial(triallist.at(i), j, function));
		}
		
		std::cout << " ** The best individual ** \n" << *(GetSelection() -> SelectBest(triallist.at(i+1))) << std::endl;
	}
	
	// Freeing memory
	std::vector<TrialT*> cleanup;
	typename std::vector<TrialT*>::iterator it;
	for(int i = 0; i <= fnGen; i++)
	{
		for(int j = 0; j < fnInd; j++)
			cleanup.push_back(triallist.at(i).at(j));
	}
	
	sort(cleanup.begin(), cleanup.end());
	it = unique(cleanup.begin(), cleanup.end());
	cleanup.resize(it - cleanup.begin());
	for(it = cleanup.begin(); it != cleanup.end(); it++)
		delete *it;
	
	cleanup.clear();

}
	
template <typename TrialT, typename MutationSchemeT, typename CrossoverSchemeT, typename SelectionSchemeT>
void DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT>::Show(std::ostream &s)
{
	fTrial -> ShowParameters(s);
	fMutation -> Show(s);
	fCrossover -> Show(s);
	fSelection -> Show(s);
}

template <typename TrialT, typename MutationSchemeT, typename CrossoverSchemeT, typename SelectionSchemeT>
std::ostream & operator<<(std::ostream &s, DE<TrialT, MutationSchemeT, CrossoverSchemeT, SelectionSchemeT> &de)
{
	de.Show(s);
	
	return s;
}

#endif // __DE__H
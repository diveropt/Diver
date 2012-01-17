#ifndef __DE__H
#define __DE__H

#include "Model.h"
#include "Trial.h"
#include "cstdlib"
#include "ctime"

template <typename MutationT, typename CrossoverT, typename SelectionT>
class DE
{
	MutationT *fMutation;
	CrossoverT *fCrossover;
	SelectionT *fSelection;	
	Model *fModel;
	
	bool IsModel;
	
public:
	//! Constructors
	DE();
	explicit DE(Model *model);
	//! Copy constructor
	DE(const DE &de);
	//! Destructor
	~DE();
	//! Copy assignment operator
	DE& operator=(DE de);
	
	// Getter
	MutationT* GetMutation() const		{return fMutation;}
	CrossoverT* GetCrossover() const	{return fCrossover;}
	SelectionT* GetSelection() const	{return fSelection;}
	Model* GetModel() const				{return fModel;}
	
	// Setters
	void SetModel(Model *model);
	
	// other member functions
	template <typename Functor>
	Trial* FirstTrial(Functor const &function);
	template <typename Functor>
	Trial* NextTrial(std::vector<Trial*> Population, const int& targetvector, Functor const &function);
	
	void Show();
};

//! Default constructor
template <typename MutationT, typename CrossoverT, typename SelectionT>
DE<MutationT, CrossoverT, SelectionT>::DE()
{
	fMutation = new MutationT();
	fCrossover = new CrossoverT();
	fSelection = new SelectionT();
	SetModel(0);
	
	srand((unsigned) time(0));	
}

//! Constructor
template <typename MutationT, typename CrossoverT, typename SelectionT>
DE<MutationT, CrossoverT, SelectionT>::DE(Model *model)
{
	fMutation = new MutationT();
	fCrossover = new CrossoverT();
	fSelection = new SelectionT();
	SetModel(model);
	
	srand((unsigned) time(0));	
}

//! Copy constructor
template <typename MutationT, typename CrossoverT, typename SelectionT>
DE<MutationT, CrossoverT, SelectionT>::DE(const DE<MutationT, CrossoverT, SelectionT> &de)
{
	fMutation = new MutationT(*de.fMutation);
	fCrossover = new CrossoverT(*de.fCrossover);
	fSelection = new SelectionT(*de.fSelection);
	SetModel(de.fModel);
	
	srand((unsigned) time(0));	
}	

//! Destructor
template <typename MutationT, typename CrossoverT, typename SelectionT>
DE<MutationT, CrossoverT, SelectionT>::~DE()
{
	delete fMutation;
	fMutation = 0;
	
	delete fCrossover;
	fCrossover = 0;
	
	delete fSelection;
	fSelection = 0;
	
	delete fModel;
	fModel = 0;
}

//! Copy assignment operator
template <typename MutationT, typename CrossoverT, typename SelectionT>
DE<MutationT, CrossoverT, SelectionT>& DE<MutationT, CrossoverT, SelectionT>::operator=(DE<MutationT, CrossoverT, SelectionT> de)
{
	fMutation = new MutationT(de.fMutation);
	fCrossover = new CrossoverT(de.fCrossover);
	fSelection = new SelectionT(de.fSelection);
	SetModel(*de.fModel);
	
	srand((unsigned) time(0));	
}

//! Sets the model
template <typename MutationT, typename CrossoverT, typename SelectionT>
void DE<MutationT, CrossoverT, SelectionT>::SetModel(Model *model)
{
	if(model)
	{
		fModel = new Model(*model);
		IsModel = true;
	}
	else
	{
		fModel = 0;
		IsModel = false;
	}
}

template <typename MutationT, typename CrossoverT, typename SelectionT>
template <typename Functor>
Trial* DE<MutationT, CrossoverT, SelectionT>::FirstTrial(Functor const &function)
{
	Trial *first = new Trial();
	first -> SetN(fModel -> GetNParameters());

	double randomnumber = 0;
	for(unsigned int i = 0; i < fModel -> GetNParameters(); i++)
	{
		randomnumber = fModel -> GetParameterLowerBound(i) + ((double)rand()/(RAND_MAX + 1.0))*(fModel -> GetParameterUpperBound(i) - fModel -> GetParameterLowerBound(i));
		first -> SetValue(i, randomnumber);
	}
	
	first -> SetFitness(function(first -> GetPoint()));

	return first;
}

template <typename MutationT, typename CrossoverT, typename SelectionT>
template <typename Functor>
Trial* DE<MutationT, CrossoverT, SelectionT>::NextTrial(std::vector<Trial*> Population, const int& targetvector, Functor const &function)
{
	// Mutate 
	fMutation -> Mutate(Population, targetvector);

	// Cross over
	fCrossover -> CrossOver(Population.at(targetvector) -> GetPoint(), fMutation -> GetDonorVector());
	
	Trial *trial = new Trial(fCrossover -> GetTrialVector());
	trial -> SetFitness(function(trial -> GetPoint()));
	
	// Select
	return fSelection -> Select(Population.at(targetvector), trial);
}

template <typename MutationT, typename CrossoverT, typename SelectionT>
void DE<MutationT, CrossoverT, SelectionT>::Show()
{
	fMutation -> Show();
	fCrossover -> Show();
	fSelection -> Show();
}

#endif // __DE__H
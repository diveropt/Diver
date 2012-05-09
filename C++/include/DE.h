#ifndef __DE__H
#define __DE__H

#include "Model.h"
#include "Trial.h"

template <typename MutationT, typename CrossoverT, typename SelectionT>
class DE
{
	MutationT *fMutation;
	CrossoverT *fCrossover;
	SelectionT *fSelection;	
	Model *fModel;
	
public:
	//! Constructors
	explicit DE(Model *model = NULL);
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
	
	void Show(std::ostream &s = std::cout);
};

//! Constructor
template <typename MutationT, typename CrossoverT, typename SelectionT>
DE<MutationT, CrossoverT, SelectionT>::DE(Model *model) : fMutation(new MutationT), fCrossover(new CrossoverT), fSelection(new SelectionT), fModel(new Model(*model))
{
	srand(time(0));	
}

//! Copy constructor
template <typename MutationT, typename CrossoverT, typename SelectionT>
DE<MutationT, CrossoverT, SelectionT>::DE(const DE<MutationT, CrossoverT, SelectionT> &de)
{
	if(fMutation)
		delete fMutation;
	fMutation = new MutationT(*de.fMutation);
	
	if(fCrossover)
		delete fCrossover;
	fCrossover = new CrossoverT(*de.fCrossover);
	
	if(fSelection)
		delete fSelection;
	fSelection = new SelectionT(*de.fSelection);
	
	SetModel(de.fModel);
	
	srand(time(0));	
}	

//! Destructor
template <typename MutationT, typename CrossoverT, typename SelectionT>
DE<MutationT, CrossoverT, SelectionT>::~DE()
{
	delete fMutation;
	fMutation = NULL;
	
	delete fCrossover;
	fCrossover = NULL;
	
	delete fSelection;
	fSelection = NULL;
	
	delete fModel;
	fModel = NULL;
}

//! Copy assignment operator
template <typename MutationT, typename CrossoverT, typename SelectionT>
DE<MutationT, CrossoverT, SelectionT>& DE<MutationT, CrossoverT, SelectionT>::operator=(DE<MutationT, CrossoverT, SelectionT> de)
{
	if(fMutation)
		delete fMutation;
	fMutation = new MutationT(de.fMutation);
	
	if(fCrossover)
		delete fCrossover;
	fCrossover = new CrossoverT(de.fCrossover);

	if(fSelection)
		delete fSelection;
	fSelection = new SelectionT(de.fSelection);

	SetModel(*de.fModel);
	
	srand(time(0));
	
	return *this;
}

//! Sets the model
template <typename MutationT, typename CrossoverT, typename SelectionT>
void DE<MutationT, CrossoverT, SelectionT>::SetModel(Model *model)
{
	if(model)
		fModel = new Model(*model);
	else
		fModel = NULL;
}

template <typename MutationT, typename CrossoverT, typename SelectionT>
template <typename Functor>
Trial* DE<MutationT, CrossoverT, SelectionT>::FirstTrial(Functor const &function)
{
	Trial *first = new Trial();
	first -> SetN(fModel -> GetNParameters());

	double randomnumber = 0;
	for(::vector_size_t i = 0; i < fModel -> GetNParameters(); i++)
	{
		randomnumber = fModel -> GetParameterLowerBound(i) + (static_cast<double>(rand())/(RAND_MAX + 1.0))*(fModel -> GetParameterUpperBound(i) - fModel -> GetParameterLowerBound(i));
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
	Trial *result = fSelection -> Select(Population.at(targetvector), trial);
	
	// Free memory
	if(result != trial)
	{
		delete trial;
		trial = NULL;
	}
	
	return result;
}

template <typename MutationT, typename CrossoverT, typename SelectionT>
void DE<MutationT, CrossoverT, SelectionT>::Show(std::ostream &s)
{
	fMutation -> Show(s);
	fCrossover -> Show(s);
	fSelection -> Show(s);
}

template <typename MutationT, typename CrossoverT, typename SelectionT>
std::ostream & operator<<(DE<MutationT, CrossoverT, SelectionT> &de, std::ostream &s)
{
	de.Show(s);
	
	return s;
}

#endif // __DE__H
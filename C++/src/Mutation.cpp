#include "Mutation.h"

Mutation& Mutation::operator=(Mutation mutation)
{
	std::swap(fScalingFactors, mutation.fScalingFactors);
	std::swap(fDonorVector, mutation.fDonorVector);
	std::swap(fIBaseVector, mutation.fIBaseVector);
	
	srand(time(0));
	
	return *this;
}

Trial* Mutation::GetBaseVector(std::vector<Trial*> Population, const int& targetvector)
{
	::vector_size_t nP = Population.size();
	
	fIBaseVector = -1;
	
	do
		fIBaseVector = rand()%nP;
	while (fIBaseVector == targetvector);
	
	return Population.at(fIBaseVector);
}

Trial* Mutation::GetBestVector(std::vector<Trial*> Population)
{
	::vector_size_t nP = Population.size();
	
	fIBestVector = -1;
	
	double fitness = 1.e99;
	for(::vector_size_t i = 0; i < Population.size(); i++)
	{
		if(Population.at(i) -> GetFitness() < fitness)
		{
			fitness = Population.at(i) -> GetFitness();
			fIBestVector = i;
		}
	}
	
	return Population.at(fIBestVector);
}

void Mutation::Mutate(std::vector<Trial*> Population, const ::vector_size_t& targetvector)
{
	std::vector<double> basevector = GetBaseVector(Population, targetvector) -> GetPoint();
	std::vector<double> bestvector = GetBestVector(Population) -> GetPoint();

	::vector_size_t nP = Population.size();
	fDonorVector.resize(basevector.size());
	
	std::vector< ::vector_size_t> iRandomVector1(fNDifferenceVectors, -1);
	std::vector< ::vector_size_t> iRandomVector2(fNDifferenceVectors, -1);

	std::vector<std::vector<double> > RandomVector1(fNDifferenceVectors);
	std::vector<std::vector<double> > RandomVector2(fNDifferenceVectors);
	
	for(::vector_size_t i = 0; i < (::vector_size_t)fNDifferenceVectors; i++)
	{
		do
		{
			iRandomVector1.at(i) = rand()%nP;
			iRandomVector2.at(i) = rand()%nP;
		}
		while(iRandomVector1.at(i) == targetvector || iRandomVector1.at(i) == iRandomVector2.at(i) || fIBaseVector == iRandomVector1.at(i));
	
		RandomVector1.at(i) = Population.at(iRandomVector1.at(i)) -> GetPoint();
		RandomVector2.at(i) = Population.at(iRandomVector2.at(i)) -> GetPoint();
	}
	
	for(::vector_size_t i = 0; i < basevector.size(); i++)
	{
		fDonorVector.at(i) = fLambda*bestvector.at(i) + (1. - fLambda)*basevector.at(i);
		
		for(::vector_size_t j = 0; j < (::vector_size_t)fNDifferenceVectors; j++)
		{
			fDonorVector.at(i) += fScalingFactors.at(j)*(RandomVector1.at(j).at(i) - RandomVector2.at(j).at(i));
		}
	}
}

void Mutation::Show(std::ostream &s)
{
	s << "=== Mutation ===" << std::endl;
	s << "Number of difference vectors: " << fNDifferenceVectors << std::endl;
	s << "Scaling factors: ";
	if(fScalingFactors.size() == 0)
		s << "empty" << std::endl;
	else
	{
		for(::vector_size_t i = 0; i < fScalingFactors.size(); i++)
			s << fScalingFactors.at(i) << "\t";
		s << std::endl;
	}
	s << "Donor vector: ";
	if(fDonorVector.size() == 0)
		s << "empty" << std::endl;
	else
	{
		for(::vector_size_t i = 0; i < fDonorVector.size(); i++)
			s << fDonorVector.at(i) << "\t";
		s << std::endl;
	}
	s << "Index of chosen base vector: " << fIBaseVector << std::endl;
	s << "Index of best vector: " << fIBestVector << std::endl;
}

std::ostream & operator<<(std::ostream &s, Mutation &mutation)
{
	mutation.Show(s);
	
	return s;
}
#include "Mutation.h"

Mutation& Mutation::operator=(Mutation mutation)
{
	std::swap(fScalingFactor, mutation.fScalingFactor);
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

void Mutation::Mutate(std::vector<Trial*> Population, const int& targetvector)
{
	std::vector<double> basevector = GetBaseVector(Population, targetvector) -> GetPoint();
	::vector_size_t nP = Population.size();
	fDonorVector.resize(basevector.size());
	
	int iRandomVector1 = -1;
	int iRandomVector2 = -1;

	do
	{
		iRandomVector1 = rand()%nP;
		iRandomVector2 = rand()%nP;
	}
	while(iRandomVector1 == targetvector || iRandomVector1 == iRandomVector2 || fIBaseVector == iRandomVector1);
	
	std::vector<double> RandomVector1 = Population.at(iRandomVector1) -> GetPoint();
	std::vector<double> RandomVector2 = Population.at(iRandomVector2) -> GetPoint();
	
	for(::vector_size_t i = 0; i < basevector.size(); i++)
		fDonorVector.at(i) = basevector.at(i) + fScalingFactor*(RandomVector1.at(i) - RandomVector2.at(i));
}

void Mutation::Show(std::ostream &s)
{
	s << "=== Mutation ===" << std::endl;
	s << "Scaling factor: " << fScalingFactor << std::endl;
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
}

std::ostream & operator<<(std::ostream &s, Mutation &mutation)
{
	mutation.Show(s);
	
	return s;
}
#include "Mutation.h"
#include "cstdlib"
#include "ctime"

Mutation& Mutation::operator=(Mutation mutation)
{
	std::swap(fScalingFactor, mutation.fScalingFactor);
	std::swap(fDonorVector, mutation.fDonorVector);
	std::swap(fIBaseVector, mutation.fIBaseVector);
	
	srand(time(0));
}

Trial* Mutation::GetBaseVector(std::vector<Trial*> Population, const int& targetvector)
{
	unsigned int nP = Population.size();
	
	fIBaseVector = -1;
	
	do
		fIBaseVector = rand()%nP;
	while (fIBaseVector == targetvector);
	
	return Population.at(fIBaseVector);
}

void Mutation::Mutate(std::vector<Trial*> Population, const int& targetvector)
{
	std::vector<double> basevector = GetBaseVector(Population, targetvector) -> GetPoint();
	unsigned int nP = Population.size();
	fDonorVector.resize(basevector.size());
	
	int iRandomVector1 = -1;
	int iRandomVector2 = -1;

	do
	{
		iRandomVector1 = rand()%nP;
		iRandomVector2 = rand()%nP;
	}
	while (iRandomVector1 == targetvector || iRandomVector1 == iRandomVector2 || fIBaseVector == iRandomVector1);
	
	std::vector<double> RandomVector1 = Population.at(iRandomVector1) -> GetPoint();
	std::vector<double> RandomVector2 = Population.at(iRandomVector2) -> GetPoint();
	
	for(unsigned int i = 0; i < basevector.size(); i++)
		fDonorVector.at(i) = basevector.at(i) + fScalingFactor*(RandomVector1.at(i) - RandomVector2.at(i));
}

void Mutation::Show()
{
	std::cout << "=== Mutation ===" << std::endl;
	std::cout << "Scaling factor: " << fScalingFactor << std::endl;
	std::cout << "Donor vector: ";
	if(fDonorVector.size() == 0)
		std::cout << "empty" << std::endl;
	else
	{
		for(unsigned int i = 0; i < fDonorVector.size(); i++)
			std::cout << fDonorVector.at(i) << "\t";
		std::cout << std::endl;
	}
	std::cout << "Index of chosen base vector: " << fIBaseVector << std::endl;
}
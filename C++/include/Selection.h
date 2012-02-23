#ifndef __SELECTION__H
#define __SELECTION__H

#include "Trial.h"

class Selection
{	
public:
	Selection() {}
	
	// Other member functions
	Trial* Select(Trial *targetvector, Trial *trialvector);
	Trial* SelectBest(std::vector<Trial*> population);
	
	void Show(std::ostream &s = std::cout);
};

std::ostream & operator<<(std::ostream &s, Selection &selection);

#endif // __SELECTION__H
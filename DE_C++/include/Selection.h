#ifndef __SELECTION__H
#define __SELECTION__H

#include <iostream>
#include "Trial.h"

class Selection
{	
public:
	Selection() {}
	
	// Other member functions
	Trial* Select(Trial *targetvector, Trial *trialvector);
	Trial* SelectBest(std::vector<Trial*> population);
	void Show();
};

#endif // __SELECTION__H
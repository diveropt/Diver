#ifndef __GLOBAL__H
#define __GLOBAL__H

#include <vector>

typedef std::vector<double>::size_type vector_size_t;

template<typename TrialT>
bool CompareFitness(TrialT* trial1, TrialT* trial2) {return (trial1 -> GetFitness() < trial2 -> GetFitness());}

#endif // __GLOBAL__H
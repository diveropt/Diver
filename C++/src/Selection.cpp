#include "Selection.h"

void Selection::Show(std::ostream &s)
{
	s << "=== Selection ===" << std::endl;
}

std::ostream & operator<<(std::ostream &s, Selection &selection)
{
	selection.Show(s);
	
	return s;
}
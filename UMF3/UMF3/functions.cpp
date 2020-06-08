#include "Structs.h" 


Func getTettaS(int i)
{
	return [](double x, double y, double t)
	{
		return 1;
	};
}

Func getTettaO(int i)
{
	if (i == 0)
		return [](double x, double y, double t)
	{
		return x * x - 1;
	};
	else
		return [](double x, double y, double t)
	{
		return 2 * y * x;
	};
}


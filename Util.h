#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <string>
#include <fstream>

namespace Util {

	std::string   getLineFromFile(std::string fileName, int lineNum);
			int   flip(int x);
	   double**   allocateMatrix(int r, int c);
		   void   deleteMatrix(double** m);
	  double***   allocateCube(int a, int b, int c);
		   void   deleteCube(double*** m, int a);

}

#endif
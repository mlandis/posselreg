#include "Util.h"





std::string Util::getLineFromFile(std::string fileName, int lineNum) {

	/* open the file */
	std::ifstream fileStream(fileName.c_str());
	if (!fileStream) 
		{
		std::cerr << "Cannot open file \"" + fileName + "\"" << std::endl;
		exit(1);
		}

	std::string linestring = "";
	int line = 0;
	while( getline(fileStream, linestring).good() )
		{
		line++;
		if (line == lineNum)
			break;
		}

	/* close the file */
	fileStream.close();
	
	if (line != lineNum)
		{
		std::cerr << "The file \"" + fileName + "\" has " << line << " lines. Could not find line " << lineNum << std::endl;
		exit(1);
		}
	
	return linestring;
}

int Util::flip(int x) {

	if (x == 0)
		return 1;
	return 0;
}

double** Util::allocateMatrix(int r, int c) {

	double **m = new double*[r];
	m[0] = new double[r * c];
	for (int i=1; i<r; i++)
		m[i] = m[i-1] + c;
	for (int i=0; i<r; i++)
		for (int j=0; j<c; j++)
			m[i][j] = 0.0;
	return m;
}

void Util::deleteMatrix(double** m) {

	delete [] m[0];
	delete [] m;
}

double*** Util::allocateCube(int a, int b, int c) {

	double*** m = new double**[a];
	for (int i=0; i<a; i++)
		{
		m[i] = new double*[b];
		m[i][0] = new double[b * c];
		for (int j=1; j<b; j++)
			m[i][j] = m[i][j-1] + c;
		}
	for (int i=0; i<a; i++)
		for (int j=0; j<b; j++)
			for (int k=0; k<c; k++)
				m[i][j][k] = 0.0;
	return m;
}

void Util::deleteCube(double*** m, int a) {

	for (int i=0; i<a; i++)
		{
		delete [] m[i][0];
		delete [] m[i];
		}
	delete [] m;
}



#include <cstdlib> 
#include <istream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "Alignment.h"
#include "code.h"
#include "iomanager.h"
#include "MbBitfield.h"



Alignment::Alignment(IoManager *iomngr) {

	if ( readAlignment(iomngr) == false )
		{
		std::cerr << "ERROR: Problem reading alignment file" << std::endl;
		exit(1);
		}
	isTranslated = false;
	geneticCode = new Code;
}

Alignment::~Alignment(void) {

	delete [] originalMatrix[0];
	delete [] originalMatrix;
	delete geneticCode;
	if (isTranslated == true)
		{
		for (int i=0; i<numTaxa; i++)
			for (int j=0; j<numCodonSites; j++)
				delete codonMatrix[i][j];
		delete [] codonMatrix[0];
		delete [] codonMatrix;
		}
}

int Alignment::getIndexForTaxon(std::string s) {

	for (int i=0; i<numTaxa; i++)
		{
		if (s == taxonNames[i])
			return i;
		}
	std::cerr << "WARNING: Could not find taxon \"" << s << "\" in alignment" << std::endl;
	return -1;
}

int Alignment::getNumSenseCodons(void) { 

	return geneticCode->getNumSenseCodons(); 
}

/*-------------------------------------------------------------------
|
|   GetPossibleNucs: 
|
|   This function initializes a vector, nuc[MAX_NUM_STATES]. The four elements
|   of nuc correspond to the four nucleotides in alphabetical order.
|   We are assuming that the nucCode is a binary representation of
|   the nucleotides that are consistent with the observation. For
|   example, if we observe an A, then the nucCode is 1 and the 
|   function initalizes nuc[0] = 1 and the other elements of nuc
|   to be 0.
|
|   Observation    nucCode        nuc
|        A            1           1000
|        C            2           0100
|        G            4           0010
|        T            8           0001
|        R            5           1010
|        Y           10           0101
|        M            3           1100
|        K           12           0011
|        S            6           0110
|        W            9           1001
|        H           11           1101
|        B           14           0111
|        V            7           1110
|        D           13           1011
|        N - ?       15           1111
|
-------------------------------------------------------------------*/
void Alignment::getPossibleNucs(int nucCode, int nuc[]) {

	if (nucCode == 1)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 2)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 3)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 4)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 5)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 6)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 7)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 8)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 9)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 10)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 11)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 12)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 13)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 14)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 15)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 16)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
}

void Alignment::print(void) {

	if (isTranslated == false)
		{
		for (int s=0; s<numSites; s++)
			{
			std::cout << std::setw(4) << s << " --";
			for (int t=0; t<numTaxa; t++)
				std::cout << std::setw(3) << originalMatrix[t][s];
			std::cout << std::endl;
			}
		}
	else
		{
		for (int s=0; s<numCodonSites; s++)
			{
			std::cout << std::setw(4) << s << " --";
			for (int t=0; t<numTaxa; t++)
				{
				for (int i=0; i<codonMatrix[t][s]->dim(); i++)
					{
					if (codonMatrix[t][s]->isBitSet(i) == true)
						std::cout << "1";
					else
						std::cout << "0";
					}
				std::cout << " ";
				//std::cout << *codonMatrix[t][s] << " ";
				}
			std::cout << std::endl;
			}
		}
}

bool Alignment::readAlignment(IoManager *iomngr) {

	/* open the file */
	ifstream seqStream;
	if ( iomngr->openFile(seqStream) == false )
		{
		std::cerr << "Cannot open file \"" + iomngr->getFileName() + "\"" << std::endl;
		exit(1);
		}
		
	string linestring = "";
	int line = 0;
	string theSequence = "";
	int taxonNum = 0;
	originalMatrix = NULL;
	numTaxa = numSites = 0;
	bool excludeLine = false, charSetLine = false;
	while( getline(seqStream, linestring).good() )
		{
		istringstream linestream(linestring);
		int ch;
		string word = "";
		int wordNum = 0;
		int siteNum = 0;
		excludeLine = false;
		charSetLine = false;
		string cmdString = "";
		do
			{
			word = "";
			linestream >> word;
			wordNum++;
			//cout << "word(" << wordNum << ") = " << word << endl;
			if (line == 0)
				{
				/* read the number of taxa/chars from the first line */
				int x;
				istringstream buf(word);
				buf >> x;
				if (wordNum == 1)
					numTaxa = x;
				else
					numSites = x;
				if (numTaxa > 0 && numSites > 0 && originalMatrix == NULL)
					{	
					originalMatrix = new int*[numTaxa];
					originalMatrix[0] = new int[numTaxa * numSites];
					for (int i=1; i<numTaxa; i++)
						originalMatrix[i] = originalMatrix[i-1] + numSites;
					for (int i=0; i<numTaxa; i++)
						for (int j=0; j<numSites; j++)
							originalMatrix[i][j] = 0;
					}
				}
			else
				{
				if (wordNum == 1)
					{
					taxonNames.push_back(word);
					taxonNum++;
					}
				else
					{
					for (int i=0; i<word.length(); i++)
						{
						char site = word.at(i);
						originalMatrix[taxonNum-1][siteNum++] = nucID(site);
						}
					}
				}
			} while ( (ch=linestream.get()) != EOF );
			
		//cout << linestring << endl;
		line++;
		}	
	
	/* close the file */
	iomngr->closeFile(seqStream);

	return true;
}

/*-------------------------------------------------------------------
|
|   NucID: 
|
|   Take a character, nuc, and return an integer:
|
|       nuc        returns
|        A            1 
|        C            2     
|        G            4      
|        T U          8     
|        R            5      
|        Y           10       
|        M            3      
|        K           12   
|        S            6     
|        W            9      
|        H           11      
|        B           14     
|        V            7      
|        D           13  
|        N - ?       15       
|
-------------------------------------------------------------------*/
int Alignment::nucID(char nuc) {

	char		n;
	
	if (nuc == 'U' || nuc == 'u')
		n = 'T';
	else
		n = nuc;

	if (n == 'A' || n == 'a')
		{
		return 1;
		}
	else if (n == 'C' || n == 'c')
		{
		return 2;
		}
	else if (n == 'G' || n == 'g')
		{
		return 4;
		}
	else if (n == 'T' || n == 't')
		{
		return 8;
		}
	else if (n == 'R' || n == 'r')
		{
		return 5;
		}
	else if (n == 'Y' || n == 'y')
		{
		return 10;
		}
	else if (n == 'M' || n == 'm')
		{
		return 3;
		}
	else if (n == 'K' || n == 'k')
		{
		return 12;
		}
	else if (n == 'S' || n == 's')
		{
		return 6;
		}
	else if (n == 'W' || n == 'w')
		{
		return 9;
		}
	else if (n == 'H' || n == 'h')
		{
		return 11;
		}
	else if (n == 'B' || n == 'b')
		{
		return 14;
		}
	else if (n == 'V' || n == 'v')
		{
		return 7;
		}
	else if (n == 'D' || n == 'd')
		{
		return 13;
		}
	else if (n == 'N' || n == 'n')
		{
		return 15;
		}
	else if (n == '-')
		{
		return 15;
		}
	else if (n == '?')
		{
		return 15;
		}
	else
		return -1;
}

bool Alignment::translate(void) {

	if (isTranslated == false)
		{
		// check that the sequence length is evenly divisible by 3
		if (numSites % 3 != 0)
			{
			std::cerr << "ERROR: Alignment is not evenly divisible by three" << std::endl;
			exit(1);
			}
		
		// translate
		numCodonSites = numSites / 3;
		codonMatrix = new MbBitfield**[numTaxa];
		codonMatrix[0] = new MbBitfield*[numTaxa * numCodonSites];
		for (int i=1; i<numTaxa; i++)
			codonMatrix[i] = codonMatrix[i-1] + numCodonSites;
		for (int i=0; i<numTaxa; i++)	
			for (int j=0; j<numCodonSites; j++)
				codonMatrix[i][j] = new MbBitfield( geneticCode->getNumSenseCodons() );
				
		for (int i=0; i<numTaxa; i++)
			{
			for (int j=0; j<numSites; j += 3)
				{
				MbBitfield *bf = codonMatrix[i][j/3];
				int nucCode1 = originalMatrix[i][j];
				int nucCode2 = originalMatrix[i][j+1];
				int nucCode3 = originalMatrix[i][j+2];
				int nuc1[4], nuc2[4], nuc3[4];
				getPossibleNucs(nucCode1, nuc1);
				getPossibleNucs(nucCode2, nuc2);
				getPossibleNucs(nucCode3, nuc3);
				int s = 0;
				for (int s1=0; s1<4; s1++)
					{
					for (int s2=0; s2<4; s2++)
						{
						for (int s3=0; s3<4; s3++)
							{
							if ( nuc1[s1] == 1 && nuc2[s2] == 1 && nuc3[s3] == 1 )
								bf->setBit(s);
							if (geneticCode->isStopCodon(s1, s2, s3) == false)	
								s++;
							}
						}
					}
				
				}
			}
		
		}
	isTranslated = true;
	return true;
}
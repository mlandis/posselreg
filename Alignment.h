#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <string>
#include <vector>

class Code;
class IoManager;
class MbBitfield;
class Alignment {

	public:
                            Alignment(IoManager *iomngr);
							~Alignment(void);
					  int   getNumTaxa(void) { return numTaxa; }
			  std::string   getNameForTaxon(int i) { return taxonNames[i]; }
			          int   getNumSenseCodons(void);
					  int   getNumCodonSites(void) { return numCodonSites; }
					  int   getIndexForTaxon(std::string s);
			  MbBitfield*   getCodonMatrixEntry(int t, int s) { return codonMatrix[t][s]; }
					  int   getOrigMatrixEntry(int t, int s) { return originalMatrix[t][s]; }
					 void   print(void);
					 bool   translate(void);
					Code*   getGeneticCode(void) { return geneticCode; }

	private:
	                 void   getPossibleNucs(int nucCode, int nuc[]);
	                  int   nucID(char nuc);
	                 bool   readAlignment(IoManager *iomngr);
	                  int   numTaxa;
					  int   numSites;
					  int   numCodonSites;
			          int   **originalMatrix;
			   MbBitfield   ***codonMatrix;
					 bool   isTranslated;
 std::vector<std::string>   taxonNames;
                     Code   *geneticCode;
};


#endif
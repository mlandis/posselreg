#ifndef MODEL_H
#define MODEL_H

#include "MbMatrix.h"
#include <string>
#include <vector>

class Alignment;
class CondLikes;
class Dnds;
class Freqs;
class MbRandom;
class MbTransitionMatrix;
class Node;
class OmegaPrior;
class Parm;
class Settings;
class SwRate;
class TiProbs;
class TiTv;
class Topology;
class Tree;

class Model {

	public:
                            Model(Alignment *ap, MbRandom *rp, Settings* sp, std::string treeStr);
							~Model(void);
					 void   calcProbs(double **dnds, double **posPr);
					 void   fillInRateMatrix(void);
				    Dnds*   getActiveDnds(void);
                   Freqs*   getActiveFreqs(void);
              OmegaPrior*   getActiveOmegaPrior(void);
                  SwRate*   getActiveSwRate(void);
                    TiTv*   getActiveTiTv(void);
				Topology*   getActiveTopology(void);
				    Tree*   getActiveTree(void);
					  int   getNumParameters(void) { return parameters.size(); }
					Parm*   getParameter(int i) { return parameters[i]; }
					  int   getNumSites(void) { return numSites; }
					 void   printAcceptanceInfo(void);
				    Parm*   pickParmAtRandom(void);
				   double   lnLikelihood(void);
				     void   printQ(void);
					 void   updateAllFlags(void);
					 void   updateQ(Parm* p);
					 void   restoreQ(Parm* p);
					 void   keep(Parm* p);
					 void   restore(Parm* p);
				   double   calcDwellTime(int i, int j, double v, int k);
				   double   calcExpectedOmega(int i, int j, double v, int k);

	private:
	                 void   initializeCondLikes(void);
	                 void   initializeRateMatrix(int nc, int no);
					 bool   isTransversion(int a, int b);
					 void   fillInStationaryProbs(void);
	            Alignment   *alignmentPtr;
	             MbRandom   *ranPtr;
	   std::vector<Parm*>   parameters;
			   CondLikes*   condLikes;
	  MbTransitionMatrix*   tiMatrix;
	             TiProbs*   tiProbs;
		 MbMatrix<double>   q;
		    MbMatrix<int>   qTemplate;
		 std::vector<int>   catOffsets;
	  std::vector<double>   proposalProbs;
	              double*   lnScalers;
	              double*   stationaryProbs;
		              int   numStates;
					  int   numSenseCodons;
					  int   numOmegaCategories;
					  int   numSites;
};


#endif
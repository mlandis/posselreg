#ifndef MCMC_H
#define MCMC_H

#include <fstream>
#include <string>

class Alignment;
class IoManager;
class MbRandom;
class Model;
class Parm;
class Settings;
class TreeSum;

class Mcmc {

	public:
                            Mcmc(Alignment* ap, Settings* sp, MbRandom* rp, Model* mp, IoManager* om);
							~Mcmc(void);

	private:
	                 void   printChainState(int n, double lnL);
	                 void   runChain(void);
				   double   safeExp(double lnX);
			  std::string   trueFalse(bool tf);
			         void   openFiles(std::string fn);
				   Model*   modelPtr;
				MbRandom*   ranPtr;
				Settings*   settingsPtr;
			   Alignment*   alignmentPtr;
					  int   numCycles;
					  int   printFrequency;
					  int   sampleFrequency;
					  int   summarizeFrequency;
			std::ofstream   treeFileStrm;
			std::ofstream   parmFileStrm;
			     TreeSum*   treeSummary;
};


#endif
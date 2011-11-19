#ifndef PARM_FREQS_H
#define PARM_FREQS_H

#include <string>
#include <vector>
#include "Parm.h"

class Model;

class Freqs {

	public:
                            Freqs(MbRandom *rp, Alignment *ap, Model* mp);
							Freqs(Freqs &c);
							~Freqs(void);
	                Freqs   &operator=(const Freqs &A);
		           double   change(void);
	 std::vector<double>&   getFreqs(void) { return f; }
		           double   lnProbability(void);
				     void   print(void);

	private:
				MbRandom*   ranPtr;
				   Model*   modelPtr;
					 void   clone(const Freqs &c);
					 void   normalizeFreqs(std::vector<double>& v, double threshhold);
					  int   numCodons;
	  std::vector<double>   a;
	  std::vector<double>   f;
};

class Alignment;
class MbRandom;

class CodonFrequencies : public Parm {

	public:
                            CodonFrequencies(MbRandom *rp, std::string pn, Model* mp, Alignment *ap);
							~CodonFrequencies(void);
				   Freqs*   getActiveFreqs(void) { return freqs[activeState]; }
		           double   lnPriorRatio(void);
		           double   change(void);
				     void   print(void);
					 void   keep(void);
					 void   restore(void);
	          std::string   getParameterStr(void);
	          std::string   getParameterHeader(void);

	private:
				    Freqs   *freqs[2];
};


#endif
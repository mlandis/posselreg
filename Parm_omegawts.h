#ifndef PARM_OMEGAWTS_H
#define PARM_OMEGAWTS_H

#include <string>
#include <vector>
#include "Parm.h"

class Model;

class OmegaPrior {

	public:
                            OmegaPrior(MbRandom *rp, int k, Model* mp);
							OmegaPrior(OmegaPrior &d);
			   OmegaPrior   &operator=(const OmegaPrior &d);
		           double   change(void);
	 std::vector<double>&   getDndsWts(void) { return dndsWts; }
		           double   lnProbability(void);
				     void   print(void);

	private:
				MbRandom*   ranPtr;
				   Model*   modelPtr;
					 void   clone(const OmegaPrior &d);
					  int   numCategories;
	  std::vector<double>   a;
	  std::vector<double>   dndsWts;
};

class Alignment;
class MbRandom;

class OmegaWts : public Parm {

	public:
                            OmegaWts(MbRandom *rp, std::string pn, Model* mp, int k);
							~OmegaWts(void);
			  OmegaPrior*   getActiveOmegaPrior(void) { return omegaPrior[activeState]; }
		           double   lnPriorRatio(void);
		           double   change(void);
				     void   print(void);
					 void   keep(void);
					 void   restore(void);
	          std::string   getParameterStr(void);
	          std::string   getParameterHeader(void);

	private:
			   OmegaPrior   *omegaPrior[2];
};

#endif
#ifndef PARM_OMEGA_H
#define PARM_OMEGA_H

#include <string>
#include <vector>
#include "Parm.h"

class Model;

class Dnds {

	public:
                            Dnds(MbRandom *rp, int k, Model* mp);
							Dnds(Dnds &d);
	                 Dnds   &operator=(const Dnds &d);
		           double   change(void);
	 std::vector<double>&   getDnds(void) { return dndsVals; }
		           double   lnProbability(void);
				     void   print(void);
	private:
				MbRandom*   ranPtr;
				   Model*   modelPtr;
					 void   clone(const Dnds &d);
					  int   numCategories;
	  std::vector<double>   dndsVals;
};

class Alignment;
class MbRandom;

class Omega : public Parm {

	public:
                            Omega(MbRandom *rp, std::string pn, Model* mp, int k);
							~Omega(void);
				    Dnds*   getActiveDnds(void) { return dnds[activeState]; }
		           double   lnPriorRatio(void);
		           double   change(void);
				     void   print(void);
					 void   keep(void);
					 void   restore(void);
	          std::string   getParameterStr(void);
	          std::string   getParameterHeader(void);

	private:
				     Dnds   *dnds[2];
};

#endif
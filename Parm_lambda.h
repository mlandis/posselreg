#ifndef PARM_LAMBDA_H
#define PARM_LAMBDA_H

#include <string>
#include "Parm.h"

class MbRandom;
class Model;

class SwRate {

	public:
                            SwRate(MbRandom *rp, Model* mp, double swLambda);
							SwRate(SwRate &r);
				   SwRate   &operator=(const SwRate &r);
		           double   change(void);
				   double   getRate(void) { return lam; }
		           double   lnProbability(void);
				     void   print(void);

	private:
					 void   clone(const SwRate &r);
				MbRandom*   ranPtr;
				   Model*   modelPtr;
				   double   lam;
				   double   exponentialPriorParm;
};

class Lambda : public Parm {

	public:
                            Lambda(MbRandom *rp, std::string pn, Model* mp, double swLambda);
							~Lambda(void);
				  SwRate*   getActiveSwRate(void) { return swrate[activeState]; }
		           double   lnPriorRatio(void);
		           double   change(void);
				     void   print(void);
					 void   keep(void);
					 void   restore(void);
	          std::string   getParameterStr(void);
	          std::string   getParameterHeader(void);

	private:
				   SwRate   *swrate[2];
};


#endif
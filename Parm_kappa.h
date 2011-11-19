#ifndef PARM_KAPPA_H
#define PARM_KAPPA_H

#include <string>
#include "Parm.h"

class MbRandom;
class Model;

class TiTv {

	public:
                            TiTv(MbRandom *rp, Model* mp);
							TiTv(TiTv &t);
	                 TiTv   &operator=(const TiTv &k);
		           double   change(void);
				   double   getRate(void) { return k; }
		           double   lnProbability(void);
				     void   print(void);

	private:
					 void   clone(const TiTv &t);
				MbRandom*   ranPtr;
				   Model*   modelPtr;
				   double   k;
};

class Kappa : public Parm {

	public:
                            Kappa(MbRandom *rp, std::string pn, Model* mp);
							~Kappa(void);
					TiTv*   getActiveTiTv(void) { return titvs[activeState]; }
		           double   lnPriorRatio(void);
		           double   change(void);
				     void   print(void);
					 void   keep(void);
					 void   restore(void);
	          std::string   getParameterStr(void);
	          std::string   getParameterHeader(void);

	private:
				    TiTv   *titvs[2];
};


#endif
#ifndef PARM_H
#define PARM_H

#include <string>

class MbRandom;
class Model;

class Parm {

	public:
                            Parm(MbRandom *rp, std::string pn, Model* mp);
		      std::string   getName(void) { return parmName; }
		   virtual double   lnPriorRatio(void) = 0;
		   virtual double   change(void) = 0;
		     virtual void   print(void) = 0;
			 virtual void   keep(void) = 0;
			 virtual void   restore(void) = 0;
			         void   incrementNumAccepted(void) { numAcceptedChanges++; }
					  int   getNumAcceptances(void) { return numAcceptedChanges; }
					  int   getNumAttempts(void) { return numAttemptedChanges; }
	  virtual std::string   getParameterStr(void) = 0;
	  virtual std::string   getParameterHeader(void) = 0;
		   
	protected:
				MbRandom*   ranPtr;
				   Model*   modelPtr;
              std::string   parmName;
			          int   activeState;
					  int   numAttemptedChanges;
					  int   numAcceptedChanges;
};


#endif
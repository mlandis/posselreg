#include "Parm.h"



Parm::Parm(MbRandom* rp, std::string pn, Model* mp) {

	ranPtr              = rp;
	modelPtr            = mp;
	parmName            = pn;
	activeState         = 0;
	numAttemptedChanges = 0;
	numAcceptedChanges  = 0;
}

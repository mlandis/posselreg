#ifndef TIPROBS_H
#define TIPROBS_H

#include "MbMatrix.h"

class MbTransitionMatrix;
class TiProbs {

	public:
                            TiProbs(int nn, MbTransitionMatrix* tm, Model* mp);
							~TiProbs(void);
		MbMatrix<double>&   getTiMatrix(int space, int node) { return tis[space][node]; }
					 void   setTransitionProbs(void);
					 void   print(void);

	private:
					 bool   isMatrixGood(int space, int node, double thresshold);
					  int   numNodes;
					  int   numStates;
	  MbTransitionMatrix*   tiMatrixPtr;
	               Model*   modelPtr;
	     MbMatrix<double>   *tis[2];
};


#endif
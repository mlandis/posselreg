#ifndef NODEVALS_H
#define NODEVALS_H

#include <vector>

class NodeVals {

	public:
                            NodeVals(void);
							~NodeVals(void);
					 void   addSample(int s, double v);
					  int   getNumSamples(void) { return count; }
					 void   print(void);
				   double   averageBrlen(int discardSamplesUpTo);
				      int   numSamplesWithThisBipartition(int discardSamplesUpTo);

	private:
	                  int   count;
		 std::vector<int>   sampleNum;
	  std::vector<double>   brlens;
};


#endif
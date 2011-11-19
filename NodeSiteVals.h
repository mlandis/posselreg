#ifndef NODESITEVALS_H
#define NODESITEVALS_H

#include <vector>

class NodeSiteVals {

	public:
                            NodeSiteVals(void);
							~NodeSiteVals(void);
					 void   addSample(int s, double p, double d);
					  int   getNumSamples(void) { return probPos.size(); }
					 void   print(void);
				   double   averageDnds(int discardSamplesUpTo);
				   double   averagePosProb(int discardSamplesUpTo);
	private:
		 std::vector<int>   sampleNum;
	  std::vector<double>   probPos;
	  std::vector<double>   aveDnds;
};


#endif
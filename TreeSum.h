#ifndef TREESUM_H
#define TREESUM_H

#include "MbBitfield.h"
#include "NodeKey.h"
#include "NodeSiteKey.h"
#include <map>
#include <sstream>
#include <string>

class Alignment;
class NodeVals;
class NodeSiteVals;
class Model;

class TreeSum {

	public:
                            TreeSum(Alignment* ap, Model* mp, std::string fn, int sf);
							~TreeSum(void);
					 void   addState(int sampleNum);
					 void   print(void);
                     void   addNewKeyValPair(MbBitfield& b, int s, NodeSiteVals* v);
			NodeSiteVals*   isKeyPresent(MbBitfield& b, int s);
			         void   printSummary(void);
				   double   getVal(MbBitfield& b, int site, bool calcDnds, double burnFraction);

	private:
	                  int   numSamples;
					  int   sampleFrequency;
	               Model*   modelPtr;
			   Alignment*   alignmentPtr;
			  std::string   sumFileName;
std::map<NodeSiteKey,NodeSiteVals*>   vals;
    std::map<NodeKey,NodeVals*>   uniqueParts;
};


#endif
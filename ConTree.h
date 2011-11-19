#ifndef CONTREE_H
#define CONTREE_H

#include "NodeKey.h"
#include "MbBitfield.h"
#include <map>
#include <string>
#include <vector>

class ConNode {

	public:
                            ConNode(void);
							~ConNode(void);
					 void   allocatePartition(int n);
				 ConNode*   getLft(void) { return lft; }
				 ConNode*   getSis(void) { return sis; }
				 ConNode*   getAnc(void) { return anc; }
				   double   getV(void) { return v; }
				   double   getProb(void) { return prob; }
			  std::string   getName(void) { return name; }
			          int   getIndex(void) { return index; }
					 bool   getFlag(void ) { return flag; }
					 void   setLft(ConNode *p) { lft = p; }
					 void   setSis(ConNode	*p) { sis = p; }
					 void   setAnc(ConNode *p) { anc = p; }
					 void   setV(double x) { v = x; }
					 void   setProb(double x) { prob = x; }
					 void   setName(std::string s) { name = s; }
					 void   setIndex(int x) { index = x; }
					 void   setFlag(bool tf) { flag = tf; }
			  MbBitfield*   getPartition(void) { return part; }

	private:
				 ConNode*   lft;
				 ConNode*   sis;
				 ConNode*   anc;
					  int   index;
				   double   v;
				   double   prob;
			  std::string   name;
					 bool   flag;
			  MbBitfield*   part;
};

class Alignment;
class MbBitfield;
class NodeVals;
class TreeSum;

class ConTree {

	public:
                            ConTree(std::map<NodeKey,NodeVals*>& parts, int ns, int dis, Alignment *ap, TreeSum* ts);
							~ConTree(void);
					 void   getDownPassSequence(void);
				 ConNode*   getDownPassNode(int i) { return downPassSequence[i]; }
				      int   getNumDownPassNodes(void) { return downPassSequence.size(); }
					 void   initializePartitions(void);
					 void   setAllNodeFlags(bool tf);
					 void   listNodes(void);
					  int   getNumTaxa(void) { return numTaxa; }
					 void   writeTree(ConNode *p, std::stringstream &ss, int site, bool calcDnds, int discardSamplesBefore);
				 ConNode*   getRoot(void) { return root; }

	private:
	             TreeSum*   treeSumPtr;
					  int   dex(ConNode* p);
	                 void   buildTreeFromPartitions(std::map<NodeKey,NodeVals*>& parts, Alignment *ap);
					 void   passDn(ConNode* p);
					 bool   allSibsFlagsAreTrue(ConNode* p);
	std::vector<ConNode*>   nodes;
	std::vector<ConNode*>   downPassSequence;
				  ConNode   *root;
				      int   numSamples;
					  int   numTaxa;
					  int   discardSamplesBefore;
};

#endif
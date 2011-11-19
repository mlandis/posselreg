#ifndef PARM_TREE_H
#define PARM_TREE_H

#include <sstream>
#include <string>
#include "Parm.h"
#include "MbBitfield.h"

class Node {

	public:
                            Node(void);
							~Node(void);
					 void   flipActiveCl(void) { (activeCl == 0 ? activeCl = 1 : activeCl = 0); }
					 void   flipActiveTi(void) { (activeTi == 0 ? activeTi = 1 : activeTi = 0); }
					Node*   getLft(void) { return lft; }
					Node*   getRht(void) { return rht; }
					Node*   getAnc(void) { return anc; }
				   double   getV(void) { return v; }
			  std::string   getName(void) { return name; }
			          int   getIndex(void) { return index; }
					  int   getActiveCl(void) { return activeCl; }
					  int   getActiveTi(void) { return activeTi; }
					 bool   getUpdateCl(void) { return updateCl; }
					 bool   getUpdateTi(void) { return updateTi; }
					 bool   getFlag(void) { return flag; }
					 void   setLft(Node *p) { lft = p; }
					 void   setRht(Node	*p) { rht = p; }
					 void   setAnc(Node *p) { anc = p; }
					 void   setV(double x) { v = x; }
					 void   setName(std::string s) { name = s; }
					 void   setIndex(int x) { index = x; }
					 void   setActiveCl(int x) { activeCl = x; }
					 void   setActiveTi(int x) { activeTi = x; }
					 void   setUpdateCl(bool tf) { updateCl = tf; }
					 void   setUpdateTi(bool tf) { updateTi = tf; }
					 void   setFlag(bool tf) { flag = tf; }
			  MbBitfield&   getPartition(void) { return *part; }

	private:
				    Node*   lft;
				    Node*   rht;
				    Node*   anc;
					  int   index;
				   double   v;
			  std::string   name;
			          int   activeCl;
					  int   activeTi;
					 bool   updateCl;
					 bool   updateTi;
					 bool   flag;
			   MbBitfield   *part;
};

class Model;

class Topology {

	public:
                            Topology(MbRandom *rp, Alignment *ap, Model* mp, double ep);
							Topology(MbRandom *rp, Alignment *ap, Model* mp, std::string treeStr, double ep);
							Topology(Topology &t);
							~Topology(void);
				 Topology   &operator=(const Topology &t);
		           double   change(void);
				    Node*   getDownPassNode(int i) { return downPassSequence[i]; }
					Node*   getRoot(void) { return root; }
					 void   getDownPassSequence(void);
					  int   getNumNodes(void) { return numNodes; }
			  std::string   getNewick(void);
		           double   lnProbability(void);
				     void   print(void);
					 void   updateAllCls(bool tf);
					 void   updateAllTis(bool tf);
					 void   flipAllActiveCls(void);
					 void   flipAllActiveTis(void);
					 void   initializeTaxonBipartitions(void);
					 void   getTaxonBipartitions(void);
					 void   printTaxonBipartitions(void);
					 void   markPathDownFromNode(Node* p);
	private:
	                 void   buildRandomTree(Alignment *ap);
					 void   buildTreeFromNewickString(Alignment *ap, std::string &ts);
					 void   clone(const Topology &t);
					  int   dex(Node* p);
					 void   passDn(Node *p, int *x);
					 void   showNodes(Node *p, int indent);
					 void   writeTree(Node *p, std::stringstream &ss);
					  int   numTaxa;
					  int   numNodes;
					 Node   *nodes;
				 MbRandom   *ranPtr;
				   Model*   modelPtr;
				     Node   *root;
					 Node   **downPassSequence;
				   double   brlenLambda;
};

class Alignment;
class MbRandom;

class Tree : public Parm {

	public:
                            Tree(MbRandom *rp, std::string pn, Model* mp, Alignment *ap, double ep);
                            Tree(MbRandom *rp, std::string pn, Model* mp, Alignment *ap, std::string treeStr, double ep);
							~Tree(void);
				Topology*   getActiveTopology(void) { return trees[activeState]; }
		           double   lnPriorRatio(void);
		           double   change(void);
				     void   print(void);
					 void   keep(void);
					 void   restore(void);
	          std::string   getParameterStr(void);
	          std::string   getParameterHeader(void);

	private:
				 Topology   *trees[2];
};


#endif
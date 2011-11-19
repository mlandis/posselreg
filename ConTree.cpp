#include <iomanip>
#include <iostream>
#include <sstream>
#include "Alignment.h"
#include "ConTree.h"
#include "MbBitfield.h"
#include "NodeVal.h"
#include "TreeSum.h"



ConNode::ConNode(void) {

	lft   = NULL;
	sis   = NULL;
	anc   = NULL;
	index = 0;
	v     = 0.0;
	prob  = 0.0;
	name  = "";
	flag  = false;
	part  = NULL;
}

ConNode::~ConNode(void) {

	if (part != NULL)
		delete part;
}

void ConNode::allocatePartition(int n) {

	if (part != NULL)
		delete part;
	part = new MbBitfield(n);
}

ConTree::ConTree(std::map<NodeKey,NodeVals*>& parts, int ns, int dis, Alignment *ap, TreeSum* ts) {

	treeSumPtr = ts;
	numTaxa = ap->getNumTaxa();
	numSamples = ns;
	discardSamplesBefore = dis;
	buildTreeFromPartitions(parts, ap);
	listNodes();
}

ConTree::~ConTree(void) {

	for (int i=0; i<nodes.size(); i++)
		delete nodes[i];
}

void ConTree::buildTreeFromPartitions(std::map<NodeKey,NodeVals*>& parts, Alignment *ap) {


	int nextNodeIndex = 0;
	// add the tip nodes
	for (int i=0; i<ap->getNumTaxa(); i++)
		{
		ConNode* p = new ConNode;
		p->setIndex(nextNodeIndex++);
		p->setName( ap->getNameForTaxon(i) );
		nodes.push_back( p );
		}
		
	// set up a star tree
	ConNode* r = new ConNode;
	r->setIndex(nextNodeIndex++);
	root = r;
	ConNode* p = nodes[0];
	r->setLft( p );
	p->setAnc( r );
	for (int i=1; i<nodes.size(); i++)
		{
		ConNode* q = nodes[i];
		p->setSis( q );
		q->setAnc( r );
		p = q;
		}
	nodes.push_back( r );
	
	// set up the consensus tree using the partitions
	for (std::map<NodeKey, NodeVals*>::iterator it = parts.begin(); it != parts.end(); it++)
		{
		NodeKey k = it->first;
		NodeVals* v = it->second;
		int num = v->numSamplesWithThisBipartition(discardSamplesBefore);
		MbBitfield* part = k.getPartition();
		double partProb = (double)num / numSamples;
		if ( partProb >= 0.5 )
			{
			//std::cout << *part << " " << partProb << " " << std::endl;
			
			// initialize the pre-order traversal sequence for the tree
			getDownPassSequence();
			
			// flag the nodes that define the tips to be grouped
			setAllNodeFlags(false);
			for (int i=0; i<part->dim(); i++)
				if (part->isBitSet(i) == true)
					nodes[i]->setFlag(true);
			for (int i=0; i<downPassSequence.size(); i++)
				{
				p = downPassSequence[i];
				if (p->getLft() != NULL && p != root)
					{
					if (allSibsFlagsAreTrue(p) == true)
						{
						p->setFlag(true);
						}
					}
				}

			//std::cout << "Before:" << std::endl;
			//listNodes();
				
			// initialize a vector containing the sister taxa
			std::vector<ConNode*> sisterList;
			for (int i=0; i<downPassSequence.size(); i++)
				{
				p = downPassSequence[i];
				if (p->getAnc() != NULL)
					{
					if (p->getFlag() == true && p->getAnc()->getFlag() == false)
						sisterList.push_back(p);
					}
				}
			/*std::cout << "Sister List: ";
			for (int i=0; i<sisterList.size(); i++)
				std::cout << sisterList[i]->getIndex() << " ";
			std::cout << std::endl;*/
				
			if (sisterList.size() > 1)
				{
				// get anc for the sister list
				ConNode* sr = sisterList[0]->getAnc();
				for (int i=1; i<sisterList.size(); i++)
					{
					if (sisterList[i]->getAnc() != sr)
						{
						std::cerr << "ERROR: Incompatible taxon bipartition detected!" << std::endl;
						exit(1);
						}
					}
					
				// remove the tips in the sister list from the tree
				for (int i=0; i<sisterList.size(); i++)	
					{
					ConNode* nde = sisterList[i];
					ConNode* ndeSis = nde->getSis();
					ConNode* ndeLftSis = NULL;
					ConNode* q = sr->getLft();
					while (q != NULL)
						{
						if (q->getSis() == nde)
							ndeLftSis = q;
						q = q->getSis();
						}
					
					if (sr->getLft() == nde)
						{
						// node to be removed is left-most above sr
						if (ndeSis == NULL)
							{
							std::cerr << "ERROR: Really confusing situation in the consensus tree!" << std::endl;
							exit(1);
							}
						sr->setLft(ndeSis);
						}
					else 
						{
						// node to be removed is not the left-most above sr
						ndeLftSis->setSis(ndeSis);
						}
					nde->setAnc(NULL);
					nde->setSis(NULL);
					}
					
				// make the sister taxa share a common ancestral node
				ConNode* ndeAnc = new ConNode;
				ndeAnc->setIndex(nextNodeIndex++);
				nodes.push_back(ndeAnc);
				ndeAnc->setLft(sisterList[0]);
				sisterList[0]->setAnc(ndeAnc);
				p = sisterList[0];
				for (int i=1; i<sisterList.size(); i++)
					{
					p->setSis(sisterList[i]);
					sisterList[i]->setAnc(ndeAnc);
					p = sisterList[i];
					}
					
				// add the sister taxa, with a common ancestor at ndeAnc, back into the tree
				ConNode* oldLft = sr->getLft();
				sr->setLft(ndeAnc);
				ndeAnc->setAnc(sr);
				ndeAnc->setSis(oldLft);
				}

			//std::cout << "After:" << std::endl;
			//listNodes();
			
			}
		
		}
		
#	if 0
	std::cout << "Partition List:" << std::endl;
	int cnt = 0;
	for (std::map<NodeKey, NodeVals*>::iterator it = parts.begin(); it != parts.end(); it++)
		{
		NodeKey k = it->first;
		NodeVals* v = it->second;
		int num = v->getNumSamples();
		MbBitfield* part = k.getPartition();
		double partProb = (double)num / numSamples;
		std::cout << cnt << " -- " << *part << " " << num << std::endl;
		cnt++;
		}
#	endif
	// scan through all of the nodes in the consensus tree and fill in the
	// average branch length and the probability of the clade
	initializePartitions();
	for (int i=0; i<nodes.size(); i++)
		{
		ConNode* p = nodes[i];
		MbBitfield* prt = p->getPartition();
		NodeKey queryKey(*prt);
		std::map<NodeKey, NodeVals*>::iterator it = parts.find( queryKey );
		if ( it != parts.end() )
			{
			NodeVals* v = it->second;
			p->setV( v->averageBrlen(0.1) );
			p->setProb( (double)v->getNumSamples()/numSamples );
			}
		else 
			{
			//std::cout << "ERROR: Could not find split " << *prt << " in the list of bipartitions" << std::endl;
			}
		}
	
}

bool ConTree::allSibsFlagsAreTrue(ConNode* p) {

	ConNode* q = p->getLft();
	bool allThrown = true;
	while(q != NULL)
		{
		if (q->getFlag() == false)
			allThrown = false;
		q = q->getSis();
		}
	return allThrown;
}

void ConTree::setAllNodeFlags(bool tf) {

	for (int i=0; i<nodes.size(); i++)
		nodes[i]->setFlag(tf);
}

void ConTree::getDownPassSequence(void) {

	downPassSequence.clear();
	passDn(root);
}

void ConTree::passDn(ConNode* p) {

	if (p != NULL)
		{
		passDn(p->getLft());
		downPassSequence.push_back(p);
		passDn(p->getSis());
		}
}

void ConTree::initializePartitions(void) {

	// allocate a partition for each node
	for (int i=0; i<nodes.size(); i++)
		{
		ConNode* p = nodes[i];
		p->allocatePartition(numTaxa);
		}
		
	// get the pre-order traversal sequence of nodes for the tree
	getDownPassSequence();
		
	// initialize the partitions for each node
	for (int i=0; i<getNumDownPassNodes(); i++)
		{
		ConNode* p = getDownPassNode(i);
		MbBitfield* prt = p->getPartition();
		if (p->getLft() == NULL)
			{
			prt->setBit(p->getIndex());
			}
		else 
			{
			ConNode* q = p->getLft();
			while (q != NULL)
				{
				MbBitfield* qPrt = q->getPartition();
				(*prt) |= (*qPrt);
				q = q->getSis();
				}
			}
		}
}

int ConTree::dex(ConNode* p) {

	if (p == NULL)
		return -1;
	return p->getIndex();
}

void ConTree::listNodes(void) {

	for (int i=0; i<nodes.size(); i++)
		{
		ConNode* p = nodes[i];
		std::cout << std::setw(4) << dex(p) << " (" << std::setw(4) << dex(p->getLft()) << std::setw(4) << dex(p->getSis()) << std::setw(4) << dex(p->getAnc()) << ") ";
		if (p->getAnc() != NULL)
			std::cout << p->getV() << " " << p->getProb();
		std::cout << std::endl;
		}
}

void ConTree::writeTree(ConNode *p, std::stringstream &ss, int site, bool calcDnds, int discardSamplesBefore) {

	ss << std::fixed << std::setprecision(5);
	
	if (p->getAnc() != NULL)
		if (p->getAnc()->getLft() == p)
			ss << "(";

	for (ConNode* q = p->getLft(); q != NULL; q = q->getSis())
		{
		if (q->getAnc() != NULL)
			if (q->getAnc()->getLft() != q)
				ss << ",";
		writeTree(q, ss, site, calcDnds, discardSamplesBefore);
		}
		
	if (p->getLft() == NULL)
		{
		double x = treeSumPtr->getVal( *(p->getPartition()), site, calcDnds, discardSamplesBefore);
		ss << p->getName() << ":" << p->getV() << ":" << x;
		}
		
	if (p->getSis() == NULL && p->getAnc() != NULL)
		{
		if (p->getAnc()->getAnc() == NULL)
			ss << ")";
		else 
			{
			double x = treeSumPtr->getVal( *(p->getPartition()), site, calcDnds, discardSamplesBefore);
			ss << ")" << ":" << p->getAnc()->getV() << ":" << x;
			}
		}
}

#include <cstdlib> 
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "Alignment.h"
#include "MbRandom.h"
#include "Parm.h"
#include "Parm_tree.h"
#include "Util.h"



Node::Node(void) {

	lft      = NULL;
	rht      = NULL;
	anc      = NULL;
	index    = 0;
	v        = 0.0;
	name     = "";
	activeCl = 0;
	activeTi = 0;
	updateCl = false;
	updateTi = false;
	flag     = false;
	part     = new MbBitfield;
}

Node::~Node(void) {

	delete part;
}

Topology::Topology(MbRandom *rp, Alignment *ap, Model* mp, double ep) {

	// set base state for some pointers
	nodes = NULL;
	downPassSequence = NULL;

	// initialize some basic tree variables
	ranPtr = rp;
	modelPtr = mp;
	numTaxa = ap->getNumTaxa();
	numNodes = 2 * numTaxa - 2;
	brlenLambda = ep;
	
	// build a random tree
	buildRandomTree(ap);
	initializeTaxonBipartitions();
	getTaxonBipartitions();
	print();
}

Topology::Topology(MbRandom *rp, Alignment *ap, Model* mp, std::string treeStr, double ep) {

	// set base state for some pointers
	nodes = NULL;
	downPassSequence = NULL;

	// initialize some basic tree variables
	ranPtr = rp;
	modelPtr = mp;
	numTaxa = ap->getNumTaxa();
	numNodes = 2 * numTaxa - 2;
	brlenLambda = ep;
	
	// build tree from string
	buildTreeFromNewickString(ap, treeStr);
	initializeTaxonBipartitions();
	getTaxonBipartitions();
	print();
}

Topology::Topology(Topology &t) {

	numTaxa = 0;
	numNodes = 0;
	nodes = NULL;
	downPassSequence = NULL;
	clone(t);
}

Topology::~Topology(void) {

	delete [] nodes;
	delete [] downPassSequence;
}

Topology& Topology::operator=(const Topology &t) {

	if (this != &t)
		clone(t);
	return *this;
}

void Topology::buildRandomTree(Alignment *ap) {

	if (numTaxa < 3)
		{
		std::cerr << "ERROR: Too few taxa!" << std::endl;
		exit(1);
		}
	
	// allocate the nodes & down pass sequence vector
	nodes = new Node[numNodes];
	downPassSequence = new Node*[numNodes];
	
	// set the node indices
	for (int i=0; i<numNodes; i++)
		nodes[i].setIndex(i);
		
	// set the taxon names
	for (int i=0; i<numTaxa; i++)
		nodes[i].setName( ap->getNameForTaxon(i) );
		
	// build the three species tree
	Node **availableNodes = new Node*[numNodes];
	int numAvailableNodes = 0;
	int nextIntNode = numTaxa;
	Node *p = &nodes[0];
	root = p;
	Node *q = &nodes[nextIntNode++];
	availableNodes[numAvailableNodes++] = q;
	p->setLft(q);
	q->setAnc(p);
	p = &nodes[1];
	availableNodes[numAvailableNodes++] = p;
	q->setLft(p);
	p->setAnc(q);
	p = &nodes[2];
	availableNodes[numAvailableNodes++] = p;
	q->setRht(p);
	p->setAnc(q);

	// add the remaining taxa to the tree
	int nextTipNode = 3;
	for (int i=3; i<numTaxa; i++)
		{
		// pick a branch to attach the new branch to
		int whichNode = (int)(ranPtr->uniformRv()*numAvailableNodes);
		p = availableNodes[whichNode];
		Node *pAnc = p->getAnc();
		Node *newTip = &nodes[nextTipNode++];
		Node *newInt = &nodes[nextIntNode++];
		if (p == pAnc->getLft())
			{
			// p is to the left of pAnc
			pAnc->setLft(newInt);
			newInt->setAnc(pAnc);
			newInt->setLft(p);
			p->setAnc(newInt);
			newInt->setRht(newTip);
			newTip->setAnc(newInt);
			}
		else 
			{
			// p is to the right of pAnc
			pAnc->setRht(newInt);
			newInt->setAnc(pAnc);
			newInt->setRht(p);
			p->setAnc(newInt);
			newInt->setLft(newTip);
			newTip->setAnc(newInt);
			}
		}
	
	delete [] availableNodes;
	
	// initialize the down pass sequence for the tree
	getDownPassSequence();
	
	// set the branch lengths to 0.1 (for now)
	for (int n=0; n<numNodes; n++)
		{
		p = &nodes[n];
		if (p->getAnc() != NULL)
			p->setV( ranPtr->exponentialRv(brlenLambda) );
		}
}

void Topology::buildTreeFromNewickString(Alignment *ap, std::string &ts) {

	/* parse the tree string, and put each token into a vector of strings */
	std::vector<std::string> parsedNewick;
	std::string temp = "";
	bool readingBrlen = false;
	int nt = 0;
	for (int i=0; i<ts.size(); i++)
		{
		char c = ts[i];
		if ( c == ' ' )
			continue;
		if ( c == '(' || c == ')' || c == ',' || c == ':' || c == ';' )
			{
			temp = c;
			parsedNewick.push_back( temp );
			if ( c == ':' )
				readingBrlen = true;
			else
				readingBrlen = false;
			}
		else
			{
			/* the character is part of a taxon name */
			int j = i;
			std::string taxonName = "";
			while ( ts[j] != '(' && ts[j] != ')' && ts[j] != ',' && ts[j] != ':' && ts[j] != ';' )
				{
				taxonName += ts[j];
				j++;
				}
			parsedNewick.push_back( taxonName );
			i = j - 1;
			if ( readingBrlen == false )
				nt++;
			readingBrlen = false;
			}
		if ( c == ';' )
			break;
		}
		
	/* check that the number of taxa in the tree description is the same as the
	   number of taxa in the alignment */
	if ( nt != numTaxa )
		{
		std::cerr << "ERROR: The tree file is not the right size" << std::endl;
		std::cout << "nt = " << nt << " numTaxa = " << numTaxa << std::endl;
		exit(1);
		}
		
	// show the tokens
	/*for (std::vector<std::string>::iterator p=parsedNewick.begin(); p != parsedNewick.end(); p++)
		std::cout << "token = \"" << (*p) << "\"" << std::endl;*/

	// allocate the nodes
	nodes = new Node[numNodes];
	downPassSequence = new Node*[numNodes];
	
	// set the node indices
	for (int i=0; i<numNodes; i++)
		nodes[i].setIndex(i);

	// initialize the branch lengths
	for (int i=0; i<numNodes; i++)
		nodes[i].setV(ranPtr->exponentialRv(brlenLambda));
		
	// set the taxon names
	for (int i=0; i<numTaxa; i++)
		nodes[i].setName( ap->getNameForTaxon(i) );
		
	/* build up the tree, using the information stored in the parsed vector */
	int nextInteriorNode = numTaxa;
	Node *p = NULL;
	for (std::vector<std::string>::iterator t=parsedNewick.begin(); t != parsedNewick.end(); t++)
		{
		//std::cout << "Token: " << (*t) << std::endl;
		if ( (*t) == "(" )
			{
			/* add a new interior node */
			if ( p == NULL )
				{
				p = &nodes[nextInteriorNode++];
				}
			else
				{
				Node *q = &nodes[nextInteriorNode++];
				if (p->getLft() == NULL)
					{
					p->setLft( q );
					q->setAnc( p );
					p = q;
					}
				else if (p->getRht() == NULL)
					{
					p->setRht( q );
					q->setAnc( p );
					p = q;
					}
				else if (p->getAnc() == NULL)
					{
					p->setAnc( q );
					q->setLft( p );
					p = q;
					p->setFlag( true );
					}
				else
					{
					std::cout << "ERROR: Problem reading the Newick-formatted tree (1)" << std::endl;
					exit(1);
					}
				}
			readingBrlen = false;
			}
		else if ( (*t) == ")" )
			{
			/* we hit a right parantheses, so we should go down the tree...unless
			   we should go up! */
			if (p->getFlag() == false && p->getAnc() != NULL)
				p = p->getAnc();
			else if (p->getFlag() == true && p->getLft() != NULL)
				p = p->getLft();
			else
				{
				std::cout << "ERROR: Problem reading the Newick-formatted tree (2)" << std::endl;
				exit(1);
				}
			readingBrlen = false;
			}
		else if ( (*t) == "," )
			{
			/* we hit a comma, so we should go down the tree...unless
			   we should go up! */
			if (p->getFlag() == false && p->getAnc() != NULL)
				p = p->getAnc();
			else if (p->getFlag() == true && p->getLft() != NULL)
				p = p->getLft();
			else
				{
				std::cout << "ERROR: Problem reading the Newick-formatted tree (3)" << std::endl;
				exit(1);
				}
			readingBrlen = false;
			}
		else if ( (*t) == ":" )
			{
			readingBrlen = true;
			}
		else if ( (*t) == ";" )
			{
			/* We are at the end of the tree description. I guess we don't have
			   to do anything. */
			}
		else
			{
			if (readingBrlen == false)
				{
				/* read in a taxon name, and add the node to the tree */
				std::string theName = (*t);
				int theIndex = ap->getIndexForTaxon( (*t) );
				if (theIndex == -1)
					{
					std::cerr << "ERROR: Could not find taxon " << (*t) << " in the alignment" << std::endl;
					exit(1);
					}

				Node *q = &nodes[theIndex];
				q->setName( theName );
				
				if (p->getLft() == NULL)
					{
					p->setLft( q );
					q->setAnc( p );
					p = q;
					}
				else if (p->getRht() == NULL)
					{
					p->setRht( q );
					q->setAnc( p );
					p = q;
					}
				else if (p->getAnc() == NULL)
					{
					p->setAnc( q );
					q->setLft( p );
					p = q;
					p->setFlag( true );
					}
				else
					{
					std::cout << "ERROR: Problem reading the Newick-formatted tree (4)" << std::endl;
					exit(1);
					}
				}
			else
				{
				/* reading a branch length */
				double x = 0.0;
				std::istringstream buf(*t);
				buf >> x;
				if (x < 0.00001)
					x = 0.0001;
				if (p->getFlag() == false)
					p->setV(x);
				else 
					p->getLft()->setV(x);
				readingBrlen = false;
				}
			}
		}

	/* set the pointer to the root */
	Node *q = p;
	while (q->getAnc() != NULL)
		q = q->getAnc();
	root = q;

	// initialize the down pass sequence for the tree
	getDownPassSequence();
	
	// set the root branch length to zero
	root->setV( 0.0 );
			
}

void Topology::clone(const Topology &t) {

	numTaxa     = t.numTaxa;
	numNodes    = t.numNodes;
	brlenLambda = t.brlenLambda;
	ranPtr      = t.ranPtr;
	modelPtr    = t.modelPtr;
	
	if (nodes == NULL)
		nodes = new Node[t.numNodes];
	if (downPassSequence == NULL)
		downPassSequence = new Node*[t.numNodes];
		
	for (int n=0; n<numNodes; n++)
		{
		Node *p = &t.nodes[n];
		Node *q = &nodes[n];
		
		if (p->getLft() == NULL)
			q->setLft(NULL);
		else 
			q->setLft( &nodes[p->getLft()->getIndex()] );

		if (p->getRht() == NULL)
			q->setRht(NULL);
		else 
			q->setRht( &nodes[p->getRht()->getIndex()] );

		if (p->getAnc() == NULL)
			q->setAnc(NULL);
		else 
			q->setAnc( &nodes[p->getAnc()->getIndex()] );

		q->setIndex( p->getIndex() );
		q->setV( p->getV() );
		q->setName( p->getName() );
		q->setActiveCl( p->getActiveCl() );
		q->setActiveTi( p->getActiveTi() );
		q->setUpdateCl( p->getUpdateCl() );
		q->setUpdateTi( p->getUpdateTi() );
		q->setFlag( p->getFlag() );
		
		if (p == t.root)
			root = q;
		
		p = t.downPassSequence[n];
		downPassSequence[n] = &nodes[p->getIndex()];
		}
}

double Topology::change(void) {

	// update conditional likelihood and transition probability flags
	updateAllCls(true);
	updateAllTis(true);
	flipAllActiveCls();
	flipAllActiveTis();

	// pick a branch at random
	Node* p = NULL;
	do
		{
		p = &nodes[(int)(ranPtr->uniformRv()*numNodes)];
		} while ( p->getAnc() == NULL );
	
	// update flags
	/*Node* q = p;
	while (q != NULL)
		{
		q->flipActiveTi();
		q->setUpdateTi(true);
		q->flipActiveCl();
		q->setUpdateCl(true);
		q = q->getAnc();
		}*/
		
	double oldV = p->getV();
	double tuning = log(4.0);
	double newV = oldV * exp( tuning*(ranPtr->uniformRv()-0.5) );
	p->setV(newV);
	
	return log(newV) - log(oldV);
}

int Topology::dex(Node *p) {

	if (p != NULL)
		return p->getIndex();
	return -1;
}

double Topology::lnProbability(void) {

	double lnP = 0.0;
	for (int n=0; n<numNodes; n++)
		{
		Node *p = &nodes[n];
		if (p->getAnc() != NULL)
			lnP += ranPtr->lnExponentialPdf( brlenLambda, p->getV() );
		}
	return lnP;
}

void Topology::updateAllCls(bool tf) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = &nodes[n];
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
			p->setUpdateCl(tf);
		}
}

void Topology::updateAllTis(bool tf) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = &nodes[n];
		if (p->getAnc() != NULL)
			p->setUpdateTi(tf);
		}
}

void Topology::flipAllActiveCls(void) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = &nodes[n];
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
			p->flipActiveCl();
		}
}

void Topology::flipAllActiveTis(void) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = &nodes[n];
		if (p->getAnc() != NULL)
			p->flipActiveTi();
		}
}

void Topology::getDownPassSequence(void) {

	int i = 0;
	passDn(root, &i);
}

void Topology::passDn(Node *p, int *x) {

	if (p != NULL)
		{
		passDn(p->getLft(), x);
		passDn(p->getRht(), x);
		downPassSequence[(*x)++] = p;
		}
}

void Topology::print(void) {

	showNodes(root, 0);
}

void Topology::showNodes(Node *p, int indent) {

	if (p != NULL)
		{
		for (int i=0; i<indent; i++)
			std::cout << " ";
		std::cout << std::setw(2) << dex(p) << " (" << std::setw(2) << dex(p->getLft()) << "," << std::setw(2) << dex(p->getRht()) << "," << std::setw(2) << dex(p->getAnc()) << ") ";
		std::cout << p->getPartition() << " ";
		std::cout << std::fixed << std::setprecision(5) << p->getV() << " ";
		std::cout << p->getName() << " ";
		if (p == root)
			std::cout << " <= Root ";
		std::cout << std::endl;
		showNodes(p->getLft(), indent+2);
		showNodes(p->getRht(), indent+2);
		}
}

std::string Topology::getNewick(void) {

	std::stringstream ss;
	writeTree(root->getLft(), ss);
	std::string newick = ss.str();
	return newick;
}

void Topology::markPathDownFromNode(Node* p) {

	for (int i=0; i<numNodes; i++)
		nodes[i].setFlag(false);
	Node* q = p;
	while (q != NULL)
		{
		q->setFlag(true);
		q = q->getAnc();
		}
}


void Topology::writeTree(Node *p, std::stringstream &ss) {

	if (p != NULL)
		{
		
		if (p->getLft() == NULL)
			{
			ss << p->getIndex() + 1 << ":" << std::fixed << std::setprecision(5) << p->getV();
			}
		else
			{
			if (p->getAnc() != NULL)
				ss << "(";
			writeTree (p->getLft(), ss);
			ss << ",";
			writeTree (p->getRht(), ss);	
			if (p->getAnc() != NULL)
				{
				if (p->getAnc()->getAnc() == NULL)
					ss << "," << p->getAnc()->getIndex() + 1 << ":" << std::fixed << std::setprecision(5) << p->getV();
				if (p->getAnc()->getAnc() != NULL)
					ss << "):" << std::fixed << std::setprecision(5) << p->getV();
				else
					ss << ")";
				}
			else
				{
				if (p->getAnc() == NULL)
					ss << ")";
				else
					ss << "):" << std::fixed << std::setprecision(5) << p->getV();
				}
			}
		}
}

void Topology::initializeTaxonBipartitions(void) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = &nodes[n];
		MbBitfield *tipBf = new MbBitfield(numTaxa);
		p->getPartition() = *tipBf;
		delete tipBf;
		if (p->getLft() == NULL || p->getRht() == NULL || p->getAnc() == NULL)
			{
			p->getPartition().setBit(p->getIndex());
			}
		}
}

void Topology::getTaxonBipartitions(void) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
			{
			MbBitfield &lp = p->getLft()->getPartition();
			MbBitfield &rp = p->getRht()->getPartition();
			MbBitfield &pp = p->getPartition();
			pp = (lp | rp);
			}
		}
	for (int n=0; n<numNodes; n++)
		{
		Node *p = getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
			{
			MbBitfield &pp = p->getPartition();
			if (pp.isBitSet(0) == true)
				{
				pp.flipBits();
				}
			}
		}
}

void Topology::printTaxonBipartitions(void) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = &nodes[n];
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
			{
			MbBitfield &pp = p->getPartition();
			std::cout << p->getIndex() << " -- " << pp << std::endl;
			}
		}
}

Tree::Tree(MbRandom *rp, std::string pn, Model* mp, Alignment *ap, double ep) : Parm(rp, pn, mp) {

	trees[0] = new Topology(rp, ap, modelPtr, ep);
	trees[1] = new Topology( *trees[0] );
}

Tree::Tree(MbRandom *rp, std::string pn, Model* mp, Alignment *ap, std::string treeStr, double ep) : Parm(rp, pn, mp) {

	trees[0] = new Topology(rp, ap, modelPtr, treeStr, ep);
	trees[1] = new Topology( *trees[0] );
}

Tree::~Tree(void) {

	delete trees[0];
	delete trees[1];
}

double Tree::change(void) {

	numAttemptedChanges++;
	return trees[activeState]->change();
}

double Tree::lnPriorRatio(void) {

	return trees[activeState]->lnProbability() - trees[Util::flip(activeState)]->lnProbability();
}

void Tree::print(void) {

	trees[activeState]->print();
}

void Tree::keep(void) {

	*trees[Util::flip(activeState)] = *trees[activeState];
}

void Tree::restore(void) {

	*trees[activeState] = *trees[Util::flip(activeState)];
}

std::string Tree::getParameterStr(void) {

	std::string pStr = trees[activeState]->getNewick();
	return pStr;
}

std::string Tree::getParameterHeader(void) {

	Topology* t = trees[activeState];
	
	std::string pStr = "";
	pStr += "#Nexus\n\n";
	pStr += "begin trees;\n";
	pStr += "   translate\n";
	int i = 1;
	for (int n=0; n<t->getNumNodes(); n++)
		{
		Node* p = t->getDownPassNode(n);
		if (p->getLft() == NULL && p->getRht() == NULL && p->getAnc() != NULL)
			{
			char temp[100];
			sprintf(temp, "   %d %s,\n", i++, p->getName().c_str() );
			pStr += temp;
			}
		else if (p->getAnc() == NULL)
			{
			char temp[100];
			sprintf(temp, "   %d %s;\n", i++, p->getName().c_str() );
			pStr += temp;
			}
		}

	return pStr;
}

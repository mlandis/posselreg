#include <string>
#include "ConTree.h"
#include "MbBitfield.h"
#include "Model.h"
#include "NodeSiteKey.h"
#include "NodeVal.h"
#include "NodeSiteVals.h"
#include "Parm_tree.h"
#include "TreeSum.h"
#include "Util.h"




TreeSum::TreeSum(Alignment* ap, Model* mp, std::string fn, int sf) {

	numSamples      = 0;
	alignmentPtr    = ap;
	modelPtr        = mp;
	sumFileName     = fn;
	sampleFrequency = sf;
}

TreeSum::~TreeSum(void) {

}

void TreeSum::addState(int sampleNum) {

	numSamples++;
	
	// get the taxon bipartitions for the tree
	Topology* t = modelPtr->getActiveTopology();
	t->getTaxonBipartitions();

	// calculate the probabilities on the tree
	double **aveDnds = Util::allocateMatrix( t->getNumNodes(), modelPtr->getNumSites() );
	double **posProb = Util::allocateMatrix( t->getNumNodes(), modelPtr->getNumSites() );
	modelPtr->calcProbs(aveDnds, posProb);
	
	// add these probabilities to the summary table
	for (int n=0; n<t->getNumNodes(); n++)
		{
		Node* p = t->getDownPassNode(n);
		if (p->getAnc() != NULL)
			{
			// get the partition for the branch
			MbBitfield &bf = p->getPartition();
			
			// potentially add the partition to the list of unique taxon bipartitions
			NodeKey queryKey(bf);
			std::map<NodeKey, NodeVals*>::iterator it = uniqueParts.find( queryKey );
			if ( it != uniqueParts.end() )
				{
				NodeVals* v = it->second;
				v->addSample( sampleNum, p->getV() );
				}
			else 
				{
				NodeVals* newNodeVal = new NodeVals;
				newNodeVal->addSample( sampleNum, p->getV() );
				uniqueParts.insert( std::make_pair(queryKey, newNodeVal) );
				}
			
			// find the partition in the list of partitions in the partitions (vals)
			for (int c=0; c<modelPtr->getNumSites(); c++)
				{
				NodeSiteVals* v = isKeyPresent(bf, c);
				if (v == NULL)
					{
					NodeSiteVals* newVal = new NodeSiteVals;
					addNewKeyValPair(bf, c, newVal);
					newVal->addSample( sampleNum, posProb[p->getIndex()][c], aveDnds[p->getIndex()][c] );
					}
				else 
					{
					v->addSample( sampleNum, posProb[p->getIndex()][c], aveDnds[p->getIndex()][c] );
					}
				}
			}
		}
		
	// free memory
	Util::deleteMatrix(aveDnds);
	Util::deleteMatrix(posProb);
}

NodeSiteVals* TreeSum::isKeyPresent(MbBitfield& b, int s) {

	NodeSiteKey queryKey(b, s);
	std::map<NodeSiteKey, NodeSiteVals*>::iterator it = vals.find( queryKey );
	if ( it != vals.end() )
		return it->second;
	return NULL;
}

void TreeSum::addNewKeyValPair(MbBitfield& b, int s, NodeSiteVals* v) {

	NodeSiteKey newKey(b, s);
	vals.insert( std::make_pair(newKey, v) );
}

void TreeSum::print(void) {

	for (std::map<NodeSiteKey, NodeSiteVals*>::iterator it = vals.begin(); it != vals.end(); it++)
		{
		NodeSiteKey k = it->first;
		NodeSiteVals *v = it->second;
		k.print();
		std::cout << std::endl;
		v->print();
		}
}

void TreeSum::printSummary(void) {

	for (int b=0; b<10; b++)
		{
		// calculate burn fraction
		double burnFraction = b * 0.1;
		int chainCycle = numSamples * sampleFrequency;
		int discardSamplesBefore = burnFraction * chainCycle;
		int totalNumSamples = (chainCycle - chainCycle * burnFraction) / sampleFrequency;
		if ( (int)(chainCycle - chainCycle * burnFraction) % sampleFrequency != 0 )
			totalNumSamples++;
		
		// construct the consensus tree
		ConTree* ct = new ConTree(uniqueParts, totalNumSamples, discardSamplesBefore, alignmentPtr, this);
	
		// open files for output
		char temp[20];
		sprintf(temp, "%d", b*10);
		std::string tempStr = temp;
		std::string fd = sumFileName + ".burn." + tempStr + ".dnds_sum";
		std::string fp = sumFileName + ".burn." + tempStr + ".prob_sum";
		std::ofstream   dStrm, pStrm;
		dStrm.open( fd.c_str(), std::ios::out );
		if ( !dStrm )
			{
			std::cerr << "Cannot open file \"" + fd + "\"" << std::endl;
			exit(1);
			}
		pStrm.open( fp.c_str(), std::ios::out );
		if ( !pStrm )
			{
			std::cerr << "Cannot open file \"" + fp + "\"" << std::endl;
			exit(1);
			}
		
		// print header info
		dStrm << "begin basstrees;" << std::endl;
		dStrm << "   dimensions=" << ct->getNumTaxa() << " ntrees=" << modelPtr->getNumSites() << ";" << std::endl;
		pStrm << "begin basstrees;" << std::endl;
		pStrm << "   dimensions=" << ct->getNumTaxa() << " ntrees=" << modelPtr->getNumSites() << ";" << std::endl;
		
		// format output
		for (int i=0; i<modelPtr->getNumSites(); i++)
			{
			std::stringstream ss;
			ct->writeTree(ct->getRoot(), ss, i, true, discardSamplesBefore);
			std::string newick = ss.str();
			dStrm << "   tree site_" << i+1 << " = " << newick << ";" << std::endl;
			}
		for (int i=0; i<modelPtr->getNumSites(); i++)
			{
			std::stringstream ss;
			ct->writeTree(ct->getRoot(), ss, i, false, discardSamplesBefore);
			std::string newick = ss.str();
			pStrm << "   tree site_" << i+1 << " = " << newick << ";" << std::endl;
			}
			
		// print end-of-file info
		dStrm << "end;" << std::endl;
		pStrm << "end;" << std::endl;
			
		// close files
		dStrm.close();
		pStrm.close();
		
		// delete the consensus tree
		delete ct;
		}
}

double TreeSum::getVal(MbBitfield& b, int site, bool calcDnds, double burnFraction) {

	NodeSiteVals* val = isKeyPresent( b, site );
	if (val == NULL)
		{
		std::cout << "ERROR: Could not find key in map" << std::endl;
		exit(1);
		}
	double theVal = 0.0;
	if (calcDnds == true)
		theVal = val->averageDnds(burnFraction);
	else 
		theVal = val->averagePosProb(burnFraction);
	return theVal;
}





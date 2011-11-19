#include <cmath>
#include <iomanip>
#include <iostream>
#include "MbRandom.h"
#include "Mcmc.h"
#include "Model.h"
#include "Parm.h"
#include "Parm_tree.h"
#include "Settings.h"
#include "TreeSum.h"



Mcmc::Mcmc(Alignment* ap, Settings* sp, MbRandom* rp, Model* mp, IoManager* om) {

	// remember important objects
	settingsPtr  = sp;
	ranPtr       = rp;
	modelPtr     = mp;
	alignmentPtr = ap;
	
	// initialize parameters of the chain
	numCycles          = sp->getChainLength();
	printFrequency     = sp->getPrintFrequency();
	sampleFrequency    = sp->getSampleFrequency();
	summarizeFrequency = sp->getSummarizeFrequency();
	
	// open the output files
	openFiles( settingsPtr->getOutputFileName() );
	
	// create object that keeps track of probabilities on the tree
	treeSummary = new TreeSum(alignmentPtr, modelPtr, settingsPtr->getOutputFileName(), sampleFrequency);
	
	// run the Markov chain
	runChain();
}

Mcmc::~Mcmc(void) {

	// close the output files
	treeFileStrm.close();
	parmFileStrm.close();
	delete treeSummary;
}

void Mcmc::runChain(void) {

	// get initial likelihood
	modelPtr->updateAllFlags();
	double oldLnL = modelPtr->lnLikelihood();
	
	for (int n=1; n<=numCycles; n++)
		{
		// pick a parameter to change
		Parm* parm = modelPtr->pickParmAtRandom();
		
		// propose a new state for the parameter
		double lnProposalRatio = parm->change();
		
		// calculate the prior ratio 
		double lnPriorRatio = parm->lnPriorRatio();
		
		// calculate the likelihood
		modelPtr->updateQ(parm);
		double lnL = modelPtr->lnLikelihood();
		double lnLikelihoodRatio = lnL - oldLnL;
		
		// calculate the acceptance probability of the proposed state
		double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
		
		// accept or reject the proposed state as the next state of the chain
		bool acceptState = false;
		if (ranPtr->uniformRv() < r)
			acceptState = true;
			
		// print information to the screen
		if ( n % printFrequency == 0 )
			{
			std::cout << std::setw(5) << n << " -- ";
			std::cout << std::fixed << std::setprecision(8) << oldLnL << " -> " << lnL;
			if (acceptState == true)
				std::cout << " -- Accepted change to " << parm->getName() << std::endl;
			else
				std::cout << " -- Rejected change to " << parm->getName() << std::endl;
			}
		
		// update the state of the chain
		if (acceptState == true)
			{
			parm->incrementNumAccepted();
			oldLnL = lnL;
			modelPtr->keep(parm);
			}
		else 
			{
			modelPtr->restore(parm);
			modelPtr->restoreQ(parm);
			}

		// print out current chain state
		if ( n == 1 || n % sampleFrequency == 0 || n == numCycles )
			printChainState(n, oldLnL);
			
		// print summary of chain
		if ( n % summarizeFrequency == 0 )
			treeSummary->printSummary();
		}
		
	treeSummary->print();
	modelPtr->printAcceptanceInfo();
}

void Mcmc::openFiles(std::string fn) {

	std::string tf = fn + ".t";
	std::string pf = fn + ".p";
	
	treeFileStrm.open( tf.c_str(), std::ios::out );
	if ( !treeFileStrm )
		{
		std::cerr << "ERROR: Problem opening tree output file" << std::endl;
		exit(1);
		}

	parmFileStrm.open( pf.c_str(), std::ios::out );
	if ( !parmFileStrm )
		{
		std::cerr << "ERROR: Problem opening parameter output file" << std::endl;
		exit(1);
		}
}

void Mcmc::printChainState(int n, double lnL) {

	if (n == 1)
		{
		std::string pHeaderStr = "";
		std::string tHeaderStr = "";
		for (int i=0; i<modelPtr->getNumParameters(); i++)
			{
			Parm* p = modelPtr->getParameter(i);
			Tree* derivedPtr = dynamic_cast<Tree*> (p);
			if (derivedPtr != 0)
				tHeaderStr += p->getParameterHeader();
			else 
				pHeaderStr += p->getParameterHeader();
			}

		treeFileStrm << tHeaderStr;
		parmFileStrm << "Cycle\tlnL\t" << pHeaderStr << std::endl;
		}
		
	std::string pStr = "";
	std::string tStr = "";
	for (int i=0; i<modelPtr->getNumParameters(); i++)
		{
		Parm* p = modelPtr->getParameter(i);
		Tree* derivedPtr = dynamic_cast<Tree*> (p);
		if (derivedPtr != 0)
			tStr += p->getParameterStr();
		else 
			pStr += p->getParameterStr();
		}
	treeFileStrm << "   tree_" << n << " = " << tStr << ";" << std::endl;
	parmFileStrm << n << '\t' << std::fixed << std::setprecision(2) << lnL << '\t' << pStr << std::endl;
	
	if (n == numCycles)
		{
		treeFileStrm << "end;" << std::endl;
		}
		
	if (n != 1)
		treeSummary->addState(n);
}

double Mcmc::safeExp(double lnX) {

	if (lnX < -300.0)
		return 0.0;
	else if (lnX > 0.0)
		return 1.0;
	else
		return exp(lnX);
}

std::string Mcmc::trueFalse(bool tf) {

	if (tf == true)
		return "True";
	return "False";
}
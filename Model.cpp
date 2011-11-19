#include <iomanip>
#include <iostream>
#include "Alignment.h"
#include "code.h"
#include "CondLikes.h"
#include "MbRandom.h"
#include "MbTransitionMatrix.h"
#include "Model.h"
#include "Parm_freqs.h"
#include "Parm_kappa.h"
#include "Parm_lambda.h"
#include "Parm_omega.h"
#include "Parm_omegawts.h"
#include "Parm_tree.h"
#include "Settings.h"
#include "TiProbs.h"
#include "Util.h"

#define ZERO_ENTRY		0
#define SYN_TI_ENTRY	1
#define SYN_TV_ENTRY	2
#define NONSYN_TI_ENTRY	3
#define NONSYN_TV_ENTRY	4

#define RESCALE_LIKES



Model::Model(Alignment *ap, MbRandom *rp, Settings* sp, std::string treeStr) {

	// keep track of pointers to important objects
	alignmentPtr = ap;
	ranPtr       = rp;
	numSites     = alignmentPtr->getNumCodonSites();
	
	// allocate a vector holding the stationary probabilities
	stationaryProbs = new double[ alignmentPtr->getNumSenseCodons() * sp->getNumOmegaCats() ];
	
	// add the parameters to the phylogenetic model
	if (treeStr == "")
		parameters.push_back( new Tree(ranPtr, "Tree", this, alignmentPtr, sp->getBrlenLambda()) );
	else
		parameters.push_back( new Tree(ranPtr, "Tree", this, alignmentPtr, treeStr, sp->getBrlenLambda()) );
	parameters.push_back( new CodonFrequencies(ranPtr, "Codon Frequencies", this, alignmentPtr) );
	parameters.push_back( new Kappa(ranPtr, "Kappa", this) );
	parameters.push_back( new Lambda(ranPtr, "Lambda", this, 2.0) );
	parameters.push_back( new Omega(ranPtr, "Omega", this, sp->getNumOmegaCats()) );
	parameters.push_back( new OmegaWts(ranPtr, "Omega Prior", this, sp->getNumOmegaCats()) );
	
	// initialize the proposal probabilities
	proposalProbs.push_back( 5.0 ); // probability of proposing a change to the tree
	proposalProbs.push_back( 1.0 ); // probability of proposing a change to the codon frequencies
	proposalProbs.push_back( 1.0 ); // probability of proposing a change to the transition/transversion rate ratio
	proposalProbs.push_back( 1.0 ); // probability of proposing a change to the switching rate
	proposalProbs.push_back( 1.0 ); // probability of proposing a change to the dN/dS rate ratio
	proposalProbs.push_back( 1.0 ); // probability of proposing a change to the prior probability of a category
	double sum = 0.0;
	for (int i=0; i<proposalProbs.size(); i++)
		sum += proposalProbs[i];
	for (int i=0; i<proposalProbs.size(); i++)
		proposalProbs[i] /= sum;

	// set up the rate matrix
	initializeRateMatrix( alignmentPtr->getNumSenseCodons(), sp->getNumOmegaCats() );

	// instantiate the conditional likelihoods
	condLikes = new CondLikes(alignmentPtr, sp->getNumOmegaCats());
	
	// instantiate the transition probabilities
	tiMatrix = new MbTransitionMatrix(q, false);
	tiProbs = new TiProbs( 2 * alignmentPtr->getNumTaxa() - 2, tiMatrix, this );
	
#	if defined(RESCALE_LIKES)
	// allocate a vector that holds (temporarily) the log scalers for log likelihood calculations
	lnScalers = new double[numSites];
#	endif

	// print the parameters
	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		(*p)->print();
}

Model::~Model(void) {

	delete condLikes;
	delete tiProbs;
	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		delete (*p);
	delete [] stationaryProbs;
#	if defined(RESCALE_LIKES)
	delete [] lnScalers;
#	endif
}

void Model::calcProbs(double **dnds, double **posPr) {

	// get a pointer to the current tree topology
	Topology* t = getActiveTopology();

	// get the active dNdS values
	std::vector<double> w = getActiveDnds()->getDnds();
	
	// update the transition probabilities, if necessary
	tiProbs->setTransitionProbs();
	MbMatrix<double> ti;
	
	// initialize the conditional likelihoods in each direction for each node on the tree
	initializeCondLikes();
	
	// allocate a square matrix to hold the joint probabilities
	double** jointProbs = Util::allocateMatrix(numOmegaCategories, numOmegaCategories);
	double*** eT = Util::allocateCube(numOmegaCategories, numOmegaCategories, numOmegaCategories);
	double*** eK = Util::allocateCube(numOmegaCategories, numOmegaCategories, numOmegaCategories);
	
	// calculate the stationary probabilities
	fillInStationaryProbs();
	
	for (int n=0; n<t->getNumNodes(); n++)
		{
		Node* p = t->getDownPassNode(n);
		if (p->getAnc() != NULL)
			{
			bool pToLeft = true;
			if (p->getAnc()->getLft() != p)
				pToLeft = false;
			ti = tiProbs->getTiMatrix( p->getActiveTi(), p->getIndex() );
			double* clP = condLikes->getClPtr( p->getActiveCl(), p->getIndex(), 0 );
			double* clA;
			if (pToLeft == true)
				clA = condLikes->getClPtrUpR( p->getAnc()->getActiveCl(), p->getAnc()->getIndex(), 0 );
			else 
				clA = condLikes->getClPtrUpL( p->getAnc()->getActiveCl(), p->getAnc()->getIndex(), 0 );

			for (int i=0; i<numOmegaCategories; i++)
				{
				for (int j=0; j<numOmegaCategories; j++)
					{
					for (int k=0; k<numOmegaCategories; k++)
						{
						eT[i][j][k] = calcDwellTime(i, j, p->getV(), k);
						eK[i][j][k] = calcExpectedOmega(i, j, p->getV(), k);
						}
					}
				}
			
			for (int c=0; c<numSites; c++)
				{
				
				for (int i=0; i<numOmegaCategories; i++)
					for (int j=0; j<numOmegaCategories; j++)
						jointProbs[i][j] = 0.0;
				for (int i=0; i<numStates; i++)
					{
					for (int j=0; j<numStates; j++)
						{
						double prob = stationaryProbs[i] * clA[i] * ti[i][j] * clP[j];
						//std::cout << std::fixed << std::setprecision(20);
						//std::cout << "prob = " << prob << std::endl;
						//std::cout << "clA[" << i << "] = " << clA[i] << std::endl;
						//std::cout << "clP[" << j << "] = " << clP[j] << std::endl;
						//std::cout << "ti[" << i << "][" << j << "] = " << ti[i][j] << std::endl;
						jointProbs[i/numSenseCodons][j/numSenseCodons] += prob;
						}
					}
				double sum = 0.0;
				for (int i=0; i<numOmegaCategories; i++)
					for (int j=0; j<numOmegaCategories; j++)
						sum += jointProbs[i][j];
				for (int i=0; i<numOmegaCategories; i++)
					for (int j=0; j<numOmegaCategories; j++)
						jointProbs[i][j] /= sum;
						
				dnds[p->getIndex()][c] = 0.0;
				posPr[p->getIndex()][c] = 0.0;
				for (int i=0; i<numOmegaCategories; i++)
					{
					for (int j=0; j<numOmegaCategories; j++)
						{
						for (int k=0; k<numOmegaCategories; k++)
							{
							dnds[p->getIndex()][c] += (eK[i][j][k] / p->getV()) * jointProbs[i][j];
							if ( w[k] > 1.0 )
								posPr[p->getIndex()][c] += (eT[i][j][k] / p->getV()) * jointProbs[i][j];
							}
						}
					}
					
				/*std::cout << "Branch/Site: " << p->getIndex() << "/" << c << ":" << std::endl;
				for (int i=0; i<numOmegaCategories; i++)
					{
					for (int j=0; j<numOmegaCategories; j++)
						std::cout << std::fixed << std::setprecision(30) << jointProbs[i][j] << " ";
					std::cout << std::endl;
					}
				std::cout << dnds[p->getIndex()][c] << " " << posPr[p->getIndex()][c] << std::endl;*/
					
				clP += numStates;
				clA += numStates;
				}
			
			}
		}
	
	
	
	// free the matrix holding the joint probabilities
	Util::deleteMatrix(jointProbs);
	Util::deleteCube(eT, numOmegaCategories);
	Util::deleteCube(eK, numOmegaCategories);
}

void Model::initializeCondLikes(void) {

	// get a pointer to the current tree topology
	Topology* t = getActiveTopology();
	
	// update the transition probabilities, if necessary
	tiProbs->setTransitionProbs();
	
	// allocate some matrices to temporarily hold the transition matrices
	MbMatrix<double> tiL;
	MbMatrix<double> tiR;
	MbMatrix<double> tiA;
	
	// calculate the conditional likelihoods down the tree
	for (int n=0; n<t->getNumNodes(); n++)
		{
		Node* p = t->getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
			{
#			if defined(RESCALE_LIKES)
			double* clS = condLikes->getClScalerPtr( p->getActiveCl(), p->getIndex(), 0 );
#			endif
			int lftIdx = p->getLft()->getIndex();
			int rhtIdx = p->getRht()->getIndex();
			tiL = tiProbs->getTiMatrix( p->getLft()->getActiveTi(), lftIdx );
			tiR = tiProbs->getTiMatrix( p->getRht()->getActiveTi(), rhtIdx );
			double* clL = condLikes->getClPtr( p->getLft()->getActiveCl(), lftIdx, 0 );
			double* clR = condLikes->getClPtr( p->getRht()->getActiveCl(), rhtIdx, 0 );
			double* clP = condLikes->getClPtr( p->getActiveCl(), p->getIndex(), 0 );
			for (int c=0; c<numSites; c++)
				{
#				if defined(RESCALE_LIKES)
				double maxCl = 0.0;
#				endif
				for (int i=0; i<numStates; i++)
					{
					double sumL = 0.0, sumR = 0.0;
					for (int j=0; j<numStates; j++)
						{
						sumL += tiL[i][j] * clL[j];
						sumR += tiR[i][j] * clR[j];
						}
					clP[i] = sumL * sumR;
#					if defined(RESCALE_LIKES)
					if (clP[i] > maxCl)
						maxCl = clP[i];
#					endif
					}
#				if defined(RESCALE_LIKES)
				double scaler = 1.0 / maxCl;
				for (int i=0; i<numStates; i++)
					clP[i] *= scaler;
				clS[c] = log(maxCl);
#				endif
				clL += numStates;
				clR += numStates;
				clP += numStates;
				}
			}
		}

	// calculate the conditional likelihoods up the tree
	for (int n=t->getNumNodes()-1; n>=0; n--)
		{
		Node* p = t->getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
			{
			bool pToLeft = true;
			if (p->getAnc()->getLft() != p)
				pToLeft = false;
#			if defined(RESCALE_LIKES)
			double* clSUpL = condLikes->getClScalerPtrUpL( p->getActiveCl(), p->getIndex(), 0 );
			double* clSUpR = condLikes->getClScalerPtrUpR( p->getActiveCl(), p->getIndex(), 0 );
#			endif
			int lftIdx = p->getLft()->getIndex();
			int rhtIdx = p->getRht()->getIndex();
			int ancIdx = p->getAnc()->getIndex();
			tiL = tiProbs->getTiMatrix( p->getLft()->getActiveTi(), lftIdx );
			tiR = tiProbs->getTiMatrix( p->getRht()->getActiveTi(), rhtIdx );
			tiA = tiProbs->getTiMatrix( p->getActiveTi(),           p->getIndex() );
			double* clL = condLikes->getClPtr( p->getLft()->getActiveCl(), lftIdx, 0 );
			double* clR = condLikes->getClPtr( p->getRht()->getActiveCl(), rhtIdx, 0 );
			double* clA;
			if (pToLeft == true)
				clA = condLikes->getClPtrUpR( p->getAnc()->getActiveCl(), ancIdx, 0 );
			else 
				clA = condLikes->getClPtrUpL( p->getAnc()->getActiveCl(), ancIdx, 0 );
			double* clPUpL = condLikes->getClPtrUpL( p->getActiveCl(), p->getIndex(), 0 );
			double* clPUpR = condLikes->getClPtrUpR( p->getActiveCl(), p->getIndex(), 0 );
			for (int c=0; c<numSites; c++)
				{
#				if defined(RESCALE_LIKES)
				double maxCl = 0.0;
#				endif
				for (int i=0; i<numStates; i++)
					{
					double sumL = 0.0, sumA = 0.0;
					for (int j=0; j<numStates; j++)
						{
						sumL += tiL[i][j] * clL[j];
						sumA += tiA[i][j] * clA[j];
						}
					clPUpL[i] = sumL * sumA;
#					if defined(RESCALE_LIKES)
					if (clPUpL[i] > maxCl)
						maxCl = clPUpL[i];
#					endif
					}
#				if defined(RESCALE_LIKES)
				double scaler = 1.0 / maxCl;
				for (int i=0; i<numStates; i++)
					clPUpL[i] *= scaler;
				clSUpL[c] = log(maxCl);
#				endif

#				if defined(RESCALE_LIKES)
				maxCl = 0.0;
#				endif
				for (int i=0; i<numStates; i++)
					{
					double sumR = 0.0, sumA = 0.0;
					for (int j=0; j<numStates; j++)
						{
						sumR += tiR[i][j] * clR[j];
						sumA += tiA[i][j] * clA[j];
						}
					clPUpR[i] = sumR * sumA;
#					if defined(RESCALE_LIKES)
					if (clPUpR[i] > maxCl)
						maxCl = clPUpR[i];
#					endif
					}
#				if defined(RESCALE_LIKES)
				scaler = 1.0 / maxCl;
				for (int i=0; i<numStates; i++)
					clPUpR[i] *= scaler;
				clSUpR[c] = log(maxCl);
#				endif

				clL += numStates;
				clR += numStates;
				clA += numStates;
				clPUpL += numStates;
				clPUpR += numStates;
				}
			}
		}



}

void Model::fillInStationaryProbs(void) {

	std::vector<double> f = getActiveFreqs()->getFreqs();
	std::vector<double> p = getActiveOmegaPrior()->getDndsWts();
	double *x = &stationaryProbs[0];
	for (int k=0; k<numOmegaCategories; k++)
		{
		for (int i=0; i<numSenseCodons; i++)
			{
			(*x) = f[i] * p[k];
			x++;
			}
		}
}

void Model::fillInRateMatrix(void) {

	// get the parameters of the model
	std::vector<double> f = getActiveFreqs()->getFreqs();
	std::vector<double> w = getActiveDnds()->getDnds();
	std::vector<double> p = getActiveOmegaPrior()->getDndsWts();
	double kappa          = getActiveTiTv()->getRate();
	double lambda         = getActiveSwRate()->getRate();
	
	// set diagonal elements of rate matrix to zero
	for (int i=0; i<numStates; i++)
		q[i][i] = 0.0;
		
	// fill in off diagonal rate matrices corresponding to switch events
	for (int k1=0; k1<numOmegaCategories; k1++)
		{
		for (int k2=0; k2<numOmegaCategories; k2++)
			{
			if (k1 != k2)
				{
				int k1Offset = k1 * numSenseCodons;
				int k2Offset = k2 * numSenseCodons;
				double r = lambda * p[k2];
				for (int s=0; s<numSenseCodons; s++)
					{
					int i = k1Offset + s;
					int j = k2Offset + s;
					q[i][j] = r;
					q[i][i] -= r;
					//std::cout << "i=" << i << " j= " << j << " r = " << r << " q[i][j] = " << q[i][j] << " q[i][i] = " << q[i][i] << std::endl;
					}
					
				}
			}
		}

	// fill in the diagonal elements of the matrix, corresponding to substitutions
#	if 1
	double averageRate = 0.0;
	for (int i=0; i<numSenseCodons; i++)
		{
		for (int j=0; j<numSenseCodons; j++)
			{
			if ( qTemplate[i][j] != ZERO_ENTRY )
				{
				if ( qTemplate[i][j] == SYN_TV_ENTRY )
					{
					double val = f[j];
					for (int k=0; k<numOmegaCategories; k++)
						{
						q[ i+catOffsets[k] ][ j+catOffsets[k] ] = val;
						q[ i+catOffsets[k] ][ i+catOffsets[k] ] -= val;
						averageRate += p[k] * f[i] * val;
						}
					}
				else if ( qTemplate[i][j] == SYN_TI_ENTRY )
					{
					double val = f[j] * kappa;
					for (int k=0; k<numOmegaCategories; k++)
						{
						q[ i+catOffsets[k] ][ j+catOffsets[k] ] = val;
						q[ i+catOffsets[k] ][ i+catOffsets[k] ] -= val;
						averageRate += p[k] * f[i] * val;
						}
					}
				else if ( qTemplate[i][j] == NONSYN_TV_ENTRY )
					{
					double val = f[j];
					for (int k=0; k<numOmegaCategories; k++)
						{
						double x = val * w[k];
						q[ i+catOffsets[k] ][ j+catOffsets[k] ] = x;
						q[ i+catOffsets[k] ][ i+catOffsets[k] ] -= x;
						averageRate += p[k] * f[i] * x;
						}
					}
				else if ( qTemplate[i][j] == NONSYN_TI_ENTRY )
					{
					double val = f[j] * kappa;
					for (int k=0; k<numOmegaCategories; k++)
						{
						double x = val * w[k];
						q[ i+catOffsets[k] ][ j+catOffsets[k] ] = x;
						q[ i+catOffsets[k] ][ i+catOffsets[k] ] -= x;
						averageRate += p[k] * f[i] * x;
						}
					}
				}
			}
		}
#	endif

	// scale the rate matrix
#	if 1
	double scaler = 1.0 / averageRate;
	for (int i=0; i<numSenseCodons; i++)
		{
		for (int j=0; j<numSenseCodons; j++)
			{
			if ( qTemplate[i][j] != ZERO_ENTRY )
				{
				for (int k=0; k<numOmegaCategories; k++)
					{
					int a = i + catOffsets[k];
					int b = j + catOffsets[k];
					double oldEntry = q[a][b];
					double newEntry = q[a][b] * scaler;
					q[a][b] = newEntry;
					q[a][a] += (oldEntry - newEntry);
					}
				}
			}
		}
#	endif

#	if 0
	double threshhold = 0.000001;
	bool isInError = false;
	for (int i=0; i<numStates; i++)
		{
		double rowSum = 0.0;
		for (int j=0; j<numStates; j++)
			rowSum += q[i][j];
		if (rowSum > threshhold || rowSum < -threshhold)
			{
			std::cout << "ERROR: Problem filling in rate matrix (" << rowSum << ", " << "i=" << i << ")" << std::endl;
			isInError = true;
			}
		}
	if (isInError == true)
		exit(1);
#	endif

#	if 0
	MbTransitionMatrix *tempMatrix = new MbTransitionMatrix(q, true);
	std::cout << "{";
	for (int i=0; i<numStates; i++)
		{
		std::cout << "{";
		for (int j=0; j<numStates; j++)
			{
			std::cout << std::fixed << std::setprecision(40) << q[i][j];
			if (j+1 < numStates)
				std::cout << ",";
			}
		std::cout << "}" << std::endl;
		if (i+1 < numStates)
			std::cout << "," << std::endl;
		}
	std::cout << "}" << std::endl;
	getchar();
	MbMatrix<double> tempP(numStates,numStates);
	tempP = tempMatrix->tiProbs(1.0, tempP);
	std::cout << tempP << std::endl;
	delete tempMatrix;
	exit(1);
#	endif
}

Freqs* Model::getActiveFreqs(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		{
		CodonFrequencies* derivedPtr = dynamic_cast<CodonFrequencies*> (*p);
		if (derivedPtr != 0)
			{
			return derivedPtr->getActiveFreqs();
			break;
			}
		}
	return NULL;
}

Dnds* Model::getActiveDnds(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		{
		Omega* derivedPtr = dynamic_cast<Omega*> (*p);
		if (derivedPtr != 0)
			{
			return derivedPtr->getActiveDnds();
			break;
			}
		}
	return NULL;
}

OmegaPrior* Model::getActiveOmegaPrior(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		{
		OmegaWts* derivedPtr = dynamic_cast<OmegaWts*> (*p);
		if (derivedPtr != 0)
			{
			return derivedPtr->getActiveOmegaPrior();
			break;
			}
		}
	return NULL;
}

SwRate* Model::getActiveSwRate(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		{
		Lambda* derivedPtr = dynamic_cast<Lambda*> (*p);
		if (derivedPtr != 0)
			{
			return derivedPtr->getActiveSwRate();
			break;
			}
		}
	return NULL;
}

TiTv* Model::getActiveTiTv(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		{
		Kappa* derivedPtr = dynamic_cast<Kappa*> (*p);
		if (derivedPtr != 0)
			{
			return derivedPtr->getActiveTiTv();
			break;
			}
		}
	return NULL;
}

Topology* Model::getActiveTopology(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		{
		Tree* derivedPtr = dynamic_cast<Tree*> (*p);
		if (derivedPtr != 0)
			{
			return derivedPtr->getActiveTopology();
			break;
			}
		}
	return NULL;
}

Tree* Model::getActiveTree(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		{
		Tree* derivedPtr = dynamic_cast<Tree*> (*p);
		if (derivedPtr != 0)
			{
			return derivedPtr;
			break;
			}
		}
	return NULL;
}

void Model::initializeRateMatrix(int nc, int no) {

	numSenseCodons     = nc;
	numOmegaCategories = no;
	numStates          = numSenseCodons * numOmegaCategories;
	
	q = MbMatrix<double>(numStates, numStates);
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			q[i][j] = 0.0;
			
	// fill in the template
	qTemplate = MbMatrix<int>(numSenseCodons, numSenseCodons);
	for (int i=0; i<numSenseCodons; i++)
		for (int j=0; j<numSenseCodons; j++)
			qTemplate[i][j] = ZERO_ENTRY;
	for (int i=0; i<numSenseCodons; i++)
		{
		int a1 = alignmentPtr->getGeneticCode()->getCodonNuc(i, 0);
		int a2 = alignmentPtr->getGeneticCode()->getCodonNuc(i, 1);
		int a3 = alignmentPtr->getGeneticCode()->getCodonNuc(i, 2);
		int aAa = alignmentPtr->getGeneticCode()->getTransCodon(i);
		for (int j=0; j<numSenseCodons; j++)
			{
			int b1 = alignmentPtr->getGeneticCode()->getCodonNuc(j, 0);
			int b2 = alignmentPtr->getGeneticCode()->getCodonNuc(j, 1);
			int b3 = alignmentPtr->getGeneticCode()->getCodonNuc(j, 2);
			int bAa = alignmentPtr->getGeneticCode()->getTransCodon(j);
			
			int numChanges = 0;
			if (a1 != b1)
				numChanges++;
			if (a2 != b2)
				numChanges++;
			if (a3 != b3)
				numChanges++;
			if (numChanges == 1)
				{
				bool isSynonymous = true;
				if (aAa != bAa)
					isSynonymous = false;
					
				bool isTransition = true;
				if ( isTransversion(a1, b1) == true )
					isTransition = false;
				
				if ( isSynonymous == true && isTransition == true )
					qTemplate[i][j] = SYN_TI_ENTRY;
				else if ( isSynonymous == true && isTransition == false )
					qTemplate[i][j] = SYN_TV_ENTRY;
				else if ( isSynonymous == false && isTransition == true )
					qTemplate[i][j] = NONSYN_TI_ENTRY;
				else if ( isSynonymous == false && isTransition == false )
					qTemplate[i][j] = NONSYN_TV_ENTRY;
				}
			}
		}

	// initialize the category offsets
	for (int k=0; k<numOmegaCategories; k++)
		catOffsets.push_back( k * numSenseCodons );

	// fill in the rate matrix
	fillInRateMatrix();
	printQ();
}

bool Model::isTransversion(int a, int b) {

	if (a == b)
		return false;
	int pA = (int)pow(2.0, (double)a);
	int pB = (int)pow(2.0, (double)b);
	if ( pA + pB == 5 || pA + pB == 10)
		return false;
	return true;
}

double Model::lnLikelihood(void) {

	// get a pointer to the current tree topology
	Topology* t = getActiveTopology();
	
	// update the transition probabilities, if necessary
	tiProbs->setTransitionProbs();
	
	// allocate some matrices to temporarily hold the transition matrices
	MbMatrix<double> tiL;
	MbMatrix<double> tiR;
	MbMatrix<double> tiA;
	
	// set the log scalers to zero
#	if defined(RESCALE_LIKES)
	for (int c=0; c<numSites; c++)
		lnScalers[c] = 0.0;
#	endif
	
	for (int n=0; n<t->getNumNodes(); n++)
		{
		Node* p = t->getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
			{
			
#			if defined(RESCALE_LIKES)
			double* clS = condLikes->getClScalerPtr( p->getActiveCl(), p->getIndex(), 0 );
#			endif

			if (p->getUpdateCl() == true)
				{
				// update the conditional likelihoods for a node
				if (p->getAnc()->getAnc() == NULL)
					{
					// three-way split at root node
					int lftIdx = p->getLft()->getIndex();
					int rhtIdx = p->getRht()->getIndex();
					int ancIdx = p->getAnc()->getIndex();
					tiL = tiProbs->getTiMatrix( p->getLft()->getActiveTi(), lftIdx );
					tiR = tiProbs->getTiMatrix( p->getRht()->getActiveTi(), rhtIdx );
					tiA = tiProbs->getTiMatrix( p->getActiveTi(), p->getIndex() );
					double* clL = condLikes->getClPtr( p->getLft()->getActiveCl(), lftIdx, 0 );
					double* clR = condLikes->getClPtr( p->getRht()->getActiveCl(), rhtIdx, 0 );
					double* clA = condLikes->getClPtr( p->getAnc()->getActiveCl(), ancIdx, 0 );
					double* clP = condLikes->getClPtr( p->getActiveCl(), p->getIndex(), 0 );
					for (int c=0; c<numSites; c++)
						{
#						if defined(RESCALE_LIKES)
						double maxCl = 0.0;
#						endif
						for (int i=0; i<numStates; i++)
							{
							double sumL = 0.0, sumR = 0.0, sumA = 0.0;
							for (int j=0; j<numStates; j++)
								{
								sumL += tiL[i][j] * clL[j];
								sumR += tiR[i][j] * clR[j];
								sumA += tiA[i][j] * clA[j];
								}
							clP[i] = sumL * sumR * sumA;
#							if defined(RESCALE_LIKES)
							if (clP[i] > maxCl)
								maxCl = clP[i];
#							endif
							}
#						if defined(RESCALE_LIKES)
						double scaler = 1.0 / maxCl;
						for (int i=0; i<numStates; i++)
							clP[i] *= scaler;
						clS[c] = log(maxCl);
#						endif
			
						clL += numStates;
						clR += numStates;
						clA += numStates;
						clP += numStates;
						}
					
					}
				else 
					{
					// two-way split
					int lftIdx = p->getLft()->getIndex();
					int rhtIdx = p->getRht()->getIndex();
					tiL = tiProbs->getTiMatrix( p->getLft()->getActiveTi(), lftIdx );
					tiR = tiProbs->getTiMatrix( p->getRht()->getActiveTi(), rhtIdx );
					double* clL = condLikes->getClPtr( p->getLft()->getActiveCl(), lftIdx, 0 );
					double* clR = condLikes->getClPtr( p->getRht()->getActiveCl(), rhtIdx, 0 );
					double* clP = condLikes->getClPtr( p->getActiveCl(), p->getIndex(), 0 );
					for (int c=0; c<numSites; c++)
						{
#						if defined(RESCALE_LIKES)
						double maxCl = 0.0;
#						endif
						for (int i=0; i<numStates; i++)
							{
							double sumL = 0.0, sumR = 0.0;
							for (int j=0; j<numStates; j++)
								{
								sumL += tiL[i][j] * clL[j];
								sumR += tiR[i][j] * clR[j];
								}
							clP[i] = sumL * sumR;
#							if defined(RESCALE_LIKES)
							if (clP[i] > maxCl)
								maxCl = clP[i];
#							endif
							}
#						if defined(RESCALE_LIKES)
						double scaler = 1.0 / maxCl;
						for (int i=0; i<numStates; i++)
							clP[i] *= scaler;
						clS[c] = log(maxCl);
#						endif

						clL += numStates;
						clR += numStates;
						clP += numStates;
						}
					}
				p->setUpdateCl(false);
				}
			
#			if defined(RESCALE_LIKES)
			// collect log scaler
			for (int c=0; c<numSites; c++)
				lnScalers[c] += clS[c];
#			endif
			
			}
		}
	
	// calculate the likelihood by averaging over the states at the base of the tree
	Node* p = t->getRoot()->getLft();
	fillInStationaryProbs();
	double* clP = condLikes->getClPtr( p->getActiveCl(), p->getIndex(), 0 );
	double lnL = 0.0;
	for (int c=0; c<numSites; c++)
		{
		double like = 0.0;
		for (int i=0; i<numStates; i++)
			like += clP[i] * stationaryProbs[i];
#		if defined(RESCALE_LIKES)
		lnL += log(like) + lnScalers[c];
#		else
		lnL += log(like);
#		endif
		clP += numStates;
		}
	
	return lnL;
}

Parm* Model::pickParmAtRandom(void) {

	double u = ranPtr->uniformRv();
	double sum = 0.0;
	for (int i=0; i<proposalProbs.size(); i++)
		{
		sum += proposalProbs[i];
		if (u < sum)
			return parameters[i];
		}
	return NULL;
}

void Model::printQ(void) {

	for (int i=0; i<numStates; i++)
		{
		for (int j=0; j<numStates; j++)
			{
			if (i == j)
				std::cout << std::fixed << std::setprecision(3) << q[i][j] << " ";
			else 
				std::cout << std::fixed << std::setprecision(3) << " " << q[i][j] << " ";

			}
		std::cout << std::endl;
		}
}

void Model::updateAllFlags(void) {

	Topology* t = getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();
	tiProbs->setTransitionProbs();
}

void Model::updateQ(Parm* p) {

	bool shouldUpdate = false;
	{
	CodonFrequencies* derivedPtr = dynamic_cast<CodonFrequencies*> (p);
	if (derivedPtr != 0)
		{
		shouldUpdate = true;
		goto funcEnd;
		}
	}

	{
	Kappa* derivedPtr = dynamic_cast<Kappa*> (p);
	if (derivedPtr != 0)
		{
		shouldUpdate = true;
		goto funcEnd;
		}
	}

	{
	Lambda* derivedPtr = dynamic_cast<Lambda*> (p);
	if (derivedPtr != 0)
		{
		shouldUpdate = true;
		goto funcEnd;
		}
	}

	{
	Omega* derivedPtr = dynamic_cast<Omega*> (p);
	if (derivedPtr != 0)
		{
		shouldUpdate = true;
		goto funcEnd;
		}
	}

	{
	OmegaWts* derivedPtr = dynamic_cast<OmegaWts*> (p);
	if (derivedPtr != 0)
		{
		shouldUpdate = true;
		goto funcEnd;
		}
	}

	funcEnd:
		if (shouldUpdate == true)
			{
			fillInRateMatrix();
			tiMatrix->updateQ(q);
			}
}

void Model::restoreQ(Parm* p) {

	bool shouldUpdate = false;
	{
	CodonFrequencies* derivedPtr = dynamic_cast<CodonFrequencies*> (p);
	if (derivedPtr != 0)
		{
		shouldUpdate = true;
		goto funcEnd;
		}
	}

	{
	Kappa* derivedPtr = dynamic_cast<Kappa*> (p);
	if (derivedPtr != 0)
		{
		shouldUpdate = true;
		goto funcEnd;
		}
	}

	{
	Lambda* derivedPtr = dynamic_cast<Lambda*> (p);
	if (derivedPtr != 0)
		{
		shouldUpdate = true;
		goto funcEnd;
		}
	}

	{
	Omega* derivedPtr = dynamic_cast<Omega*> (p);
	if (derivedPtr != 0)
		{
		shouldUpdate = true;
		goto funcEnd;
		}
	}

	{
	OmegaWts* derivedPtr = dynamic_cast<OmegaWts*> (p);
	if (derivedPtr != 0)
		{
		shouldUpdate = true;
		goto funcEnd;
		}
	}

	funcEnd:
		if (shouldUpdate == true)
			{
			tiMatrix->restoreQ();
			}
}

void Model::keep(Parm* p) {

	p->keep();

	Tree* derivedPtr = dynamic_cast<Tree*> (p);
	if (derivedPtr == 0)
		{
		getActiveTree()->keep();
		}
}

void Model::restore(Parm* p) {

	p->restore();

	Tree* derivedPtr = dynamic_cast<Tree*> (p);
	if (derivedPtr == 0)
		{
		getActiveTree()->restore();
		}

}

void Model::printAcceptanceInfo(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		{
		int ni = (*p)->getNumAcceptances();
		int n  = (*p)->getNumAttempts();
		double a = (double)ni / n;
		std::cout << std::setw(5) << ni << " " << std::setw(5) << n << " -- ";
		std::cout << std::fixed << std::setprecision(2) << a << " ";
		std::cout << (*p)->getName() << std::endl;
		}
}

double Model::calcDwellTime(int i, int j, double v, int k) {

	double lam = getActiveSwRate()->getRate();
	std::vector<double> p = getActiveOmegaPrior()->getDndsWts();
	double lamV = lam * v;
	double expLamV = exp(lamV);
	
	if (i == k && j == k)
		{
		// 111
		double p1 = p[i];
		return ((-1.0 + p1)*(2.0*p1 + lamV*(-1.0 + p1)) + expLamV*p1*(2.0 + p1*(-2.0 + lamV)))/(lam*(1.0 + (-1.0 + expLamV)* p1));
		}
	else if (i == j && i != k)
		{
		// 121
		double p1 = p[i];
		double p2 = p[k];
		return (p1*p2*(2.0 + lamV + expLamV*(-2.0 + lamV)))/(lam*(1.0 + (-1.0 + expLamV)*p1));
		}
	else if (i != j && i == k)
		{
		// 112
		double p1 = p[i];
		return ((-1.0 + 2.0*p1)*v)/(-1.0 + expLamV) + (1.0 - 2.0*p1 + lamV*p1)/lam;
		}
	else if (i != j && j == k)
		{
		// 122
		double p2 = p[j];
		return ((-1.0 + 2.0*p2)*v)/(-1.0 + expLamV) + (1.0 - 2.0*p2 + lamV*p2)/lam;
		}
	else if (i != j && j != k && i != k)
		{
		// 123
		double p2 = p[k];
		return (p2*(-2.0 + lamV*((expLamV+1.0)/(expLamV-1.0))))/lam;
		}
	return 0.0;
}

double Model::calcExpectedOmega(int i, int j, double v, int k) {

	double lam = getActiveSwRate()->getRate();
	std::vector<double> p = getActiveOmegaPrior()->getDndsWts();
	std::vector<double> w = getActiveDnds()->getDnds();
	double lamV = lam * v;
	double expLamV = exp(lamV);
	
	if (i == k && j == k)
		{
		// 111
		double p1 = p[i];
		double w1 = w[i];
		return ((-1.0 + p1)*(2.0*p1 + lamV*(-1.0 + p1))*w1 + expLamV*p1*(2.0 + p1*(-2.0 + lamV))*w1)/(lam*(1.0 + (-1.0 + expLamV)*p1));
		}
	else if (i == j && i != k)
		{
		// 121
		double p1 = p[i];
		double p2 = p[k];
		double w2 = w[k];
		return (p1*p2*(2.0 + lamV + expLamV*(-2.0 + lamV))*w2)/(lam*(1.0 + (-1.0 + expLamV)*p1));
		}
	else if (i != j && i == k)
		{
		// 112
		double p1 = p[i];
		double w1 = w[i];
		return ((-1.0 + 2.0*p1)*v*w1)/(-1.0 + expLamV) + ((1.0 + p1*(-2.0 + lamV))*w1)/lam;
		}
	else if (i != j && j == k)
		{
		// 122
		double p2 = p[j];
		double w2 = w[j];
		return ((-1.0 + 2.0*p2)*v*w2)/(-1.0 + expLamV) + ((1.0 + p2*(-2.0 + lamV))*w2)/lam;
		}
	else if (i != j && j != k && i != k)
		{
		// 123
		double p2 = p[k];
		double w2 = w[k];
		return (p2*w2*(-2.0 + lamV*((expLamV+1.0)/(expLamV-1.0))))/lam;
		}
	return 0.0;
}





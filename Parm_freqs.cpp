#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "Alignment.h"
#include "MbRandom.h"
#include "Model.h"
#include "Parm_freqs.h"
#include "Parm_tree.h"
#include "Util.h"



Freqs::Freqs(MbRandom *rp, Alignment *ap, Model* mp) {

	ranPtr = rp;
	modelPtr = mp;
	numCodons = ap->getNumSenseCodons();
	f.resize(numCodons);
	a.resize(numCodons);
	for (int i=0; i<numCodons; i++)
		a[i] = 1.0;
	ranPtr->dirichletRv(a, f);
	print();
}

Freqs::Freqs(Freqs &c) {

	clone(c);
}

Freqs::~Freqs(void) {

}

Freqs& Freqs::operator=(const Freqs &f) {

	if (this != &f)
		clone(f);
	return *this;
}

double Freqs::change(void) {
	
	// update conditional likelihood and transition probability flags
	Topology* t = modelPtr->getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();
	
	// tuning parameters of the move
	int numChangedFreqs = 5;
	double alpha0 = 500.0;
	
	// pick some frequencies that will be changed
	std::vector<int> temp(numCodons);
	for (int i=0; i<numCodons; i++)
		temp[i] = i;
	std::vector<int> freqIds(numChangedFreqs);
	for (int i=0; i<numChangedFreqs; i++)
		{
		int whichFreq = (int)(ranPtr->uniformRv()*(numCodons-i));
		freqIds[i] = temp[whichFreq];
		temp[whichFreq] = temp[numCodons-i-1];
		temp[numCodons-i-1] = freqIds[i];
		}
		
	// fill in a vector of the current frequencies
	std::vector<double> curF(numChangedFreqs+1);
	std::vector<double> newF(numChangedFreqs+1);
	for (int i=0; i<numChangedFreqs+1; i++)
		{
		curF[i] = 0.0;
		newF[i] = 0.0;
		}
	for (int i=0; i<numCodons; i++)
		{
		// is the codon one of the codons to have its frequency changed?
		bool isCodonSpecial = false;
		int whichOffsetInFreqIds;
		for (int j=0; j<numChangedFreqs; j++)
			{
			if (freqIds[j] == i)
				{
				isCodonSpecial = true;
				whichOffsetInFreqIds = j;
				break;
				}
			}
		// fill in the frequency
		if (isCodonSpecial == true)
			curF[whichOffsetInFreqIds] = f[i];
		else 
			curF[numChangedFreqs] += f[i];
		}

	// fill in the dirichlet prior parameter
	std::vector<double> forwardA(numChangedFreqs+1);
	std::vector<double> reverseA(numChangedFreqs+1);
	for (int i=0; i<numChangedFreqs+1; i++)
		forwardA[i] = curF[i] * alpha0;
		
	// choose new frequencies
	ranPtr->dirichletRv(forwardA, newF);
	for (int i=0; i<numChangedFreqs+1; i++)
		reverseA[i] = newF[i] * alpha0;
	normalizeFreqs(newF, 0.000001);
	
	// fill in the new frequencies
	for (int i=0; i<numCodons; i++)
		{
		// is the codon one of the codons to have its frequency changed?
		bool isCodonSpecial = false;
		int whichOffsetInFreqIds;
		for (int j=0; j<numChangedFreqs; j++)
			{
			if (freqIds[j] == i)
				{
				isCodonSpecial = true;
				whichOffsetInFreqIds = j;
				break;
				}
			}
		// fill in the frequency
		if (isCodonSpecial == true)
			f[i] = newF[whichOffsetInFreqIds];
		else 
			f[i] *= newF[numChangedFreqs] / curF[numChangedFreqs];
		}
		
#	if 0
	double sum1 = 0.0, sum2 = 0.0;
	std::cerr << "Frequencies:" << std::endl;
	for (int i=0; i<freqIds.size(); i++)
		{
		std::cerr << std::setw(2) << i << " -- " << std::setw(2) << freqIds[i] << " " << std::fixed << std::setprecision(10) << curF[i] << " " << std::fixed << std::setprecision(10) << newF[i] << std::endl;
		sum1 += curF[i];
		sum2 += newF[i];
		}
	sum1 += curF[numChangedFreqs];
	sum2 += newF[numChangedFreqs];
	std::cerr << std::setw(2) << numChangedFreqs << " -- " << "  " << " " << std::fixed << std::setprecision(10) << curF[numChangedFreqs] << " " << std::fixed << std::setprecision(10) << newF[numChangedFreqs] << std::endl;
	std::cerr << std::setw(2) << "  " << "    " << "  " << " " << std::fixed << std::setprecision(10) << sum1 << " " << std::fixed << std::setprecision(10) << sum2 << std::endl;
#	endif
	
	// calculate the hastings ratio/Jacobian
	double x = ranPtr->lnDirichletPdf(reverseA, curF) - ranPtr->lnDirichletPdf(forwardA, newF) + (numCodons-numChangedFreqs)*log(newF[numChangedFreqs]/curF[numChangedFreqs]);
	return x;
}

double Freqs::lnProbability(void) {

	return ranPtr->lnDirichletPdf(a, f);
}

void Freqs::normalizeFreqs(std::vector<double>& v, double threshhold) {

	double sum1 = 0.0, sum2 = 0.0;
	for (int i=0; i<v.size(); i++)
		{
		if (v[i] >= threshhold)
			sum1 += v[i];
		else
			sum2 += threshhold;
		}
	double scaler = (1.0 - sum2) / sum1;
	for (int i=0; i<v.size(); i++)
		{
		if (v[i] < threshhold)
			v[i] = threshhold;
		else 
			v[i] *= scaler;
		}
}

void Freqs::print(void) {

	std::cout << "Codon frequencies: ";
	for (int i=0; i<numCodons; i++)
		std::cout << std::fixed << std::setprecision(10) << f[i] << " ";
	std::cout << std::endl;
}

void Freqs::clone(const Freqs &c) {

	ranPtr    = c.ranPtr;
	modelPtr  = c.modelPtr;
	numCodons = c.numCodons;
	f.resize(c.f.size());
	a.resize(c.a.size());
	for (int i=0; i<numCodons; i++)
		{
		a[i] = c.a[i];
		f[i] = c.f[i];
		}
}

CodonFrequencies::CodonFrequencies(MbRandom *rp, std::string pn, Model* mp, Alignment *ap) : Parm(rp, pn, mp) {

	freqs[0] = new Freqs(rp, ap, modelPtr);
	freqs[1] = new Freqs( *freqs[0] );
}

CodonFrequencies::~CodonFrequencies(void) {
	
	delete freqs[0];
	delete freqs[1];
}

double CodonFrequencies::lnPriorRatio(void) {

	return freqs[activeState]->lnProbability() - freqs[Util::flip(activeState)]->lnProbability();
}

double CodonFrequencies::change(void) {

	numAttemptedChanges++;
	return freqs[activeState]->change();
}

void CodonFrequencies::print(void) {

	return freqs[activeState]->print();
}

void CodonFrequencies::keep(void) {

	*freqs[Util::flip(activeState)] = *freqs[activeState];
}

void CodonFrequencies::restore(void) {

	*freqs[activeState] = *freqs[Util::flip(activeState)];
}

std::string CodonFrequencies::getParameterStr(void) {

	 std::vector<double> f = freqs[activeState]->getFreqs();
	 std::string pStr = "";
	 for (int i=0; i<f.size(); i++)
		{
		char tempCh[50];
		sprintf(tempCh, "%1.3lf\t", f[i]);
		std::string tempStr = tempCh;
		pStr += tempStr;
		}
	return pStr;
}

std::string CodonFrequencies::getParameterHeader(void) {

	 std::vector<double> f = freqs[activeState]->getFreqs();
	 std::string pStr = "";
	 for (int i=0; i<f.size(); i++)
		{
		char tempCh[50];
		sprintf(tempCh, "Pi(%d)\t", i+1);
		std::string tempStr = tempCh;
		pStr += tempStr;
		}
	return pStr;
}




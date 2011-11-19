#include <iomanip>
#include <iostream>
#include "MbRandom.h"
#include "Model.h"
#include "Parm_omegawts.h"
#include "Parm_tree.h"
#include "Util.h"



OmegaPrior::OmegaPrior(MbRandom *rp, int k, Model* mp) {

	ranPtr        = rp;
	modelPtr      = mp;
	numCategories = k;
	dndsWts.resize(numCategories);
	a.resize(numCategories);
	
	for (int i=0; i<numCategories; i++)
		a[i] = 1.0;
	ranPtr->dirichletRv(a, dndsWts);
}

OmegaPrior::OmegaPrior(OmegaPrior &d) {

	clone(d);
}

OmegaPrior& OmegaPrior::operator=(const OmegaPrior &d) {

	if (this != &d)
		clone(d);
	return *this;
}

double OmegaPrior::change(void) {

	// update conditional likelihood and transition probability flags
	Topology* t = modelPtr->getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();

	double alpha0 = 100.0;
	
	std::vector<double> forwardA(numCategories);
	std::vector<double> reverseA(numCategories);
	std::vector<double> newV(numCategories);

	for (int i=0; i<numCategories; i++)
		forwardA[i] = dndsWts[i] * alpha0;
	ranPtr->dirichletRv(forwardA, newV);
	
	for (int i=0; i<numCategories; i++)
		reverseA[i] = newV[i] * alpha0;
		
	double lnHastingsRatio = ranPtr->lnDirichletPdf(reverseA, dndsWts) - ranPtr->lnDirichletPdf(forwardA, newV);
		
	for (int i=0; i<numCategories; i++)
		dndsWts[i] = newV[i];
	
	return lnHastingsRatio;
}

double OmegaPrior::lnProbability(void) {

	return ranPtr->lnDirichletPdf(a, dndsWts);
}

void OmegaPrior::print(void) {

	std::cout << "Omega Weights: ";
	for (int i=0; i<numCategories; i++)
		std::cout << std::fixed << std::setprecision(10) << dndsWts[i] << " ";
	std::cout << std::endl;
}

void OmegaPrior::clone(const OmegaPrior &d) {

	ranPtr = d.ranPtr;
	modelPtr = d.modelPtr;
	numCategories = d.numCategories;
	dndsWts.resize(numCategories);
	a.resize(numCategories);
	for (int i=0; i<numCategories; i++)
		{
		a[i] = d.a[i];
		dndsWts[i] = d.dndsWts[i];
		}
}

OmegaWts::OmegaWts(MbRandom *rp, std::string pn, Model* mp, int k) : Parm(rp, pn, mp) {

	omegaPrior[0] = new OmegaPrior(ranPtr, k, modelPtr);
	omegaPrior[1] = new OmegaPrior( *omegaPrior[0] );
}

OmegaWts::~OmegaWts(void) {

	delete omegaPrior[0];
	delete omegaPrior[1];
}

double OmegaWts::lnPriorRatio(void) {

	return omegaPrior[activeState]->lnProbability() - omegaPrior[Util::flip(activeState)]->lnProbability();
}

double OmegaWts::change(void) {

	numAttemptedChanges++;
	return omegaPrior[activeState]->change();
}

void OmegaWts::print(void) {

	omegaPrior[activeState]->print();
}

void OmegaWts::keep(void) {

	*omegaPrior[Util::flip(activeState)] = *omegaPrior[activeState];
}

void OmegaWts::restore(void) {

	*omegaPrior[activeState] = *omegaPrior[Util::flip(activeState)];
}

std::string OmegaWts::getParameterStr(void) {

	 std::vector<double> w = omegaPrior[activeState]->getDndsWts();
	 std::string pStr = "";
	 for (int i=0; i<w.size(); i++)
		{
		char tempCh[50];
		sprintf(tempCh, "%1.3lf\t", w[i]);
		std::string tempStr = tempCh;
		pStr += tempStr;
		}
	return pStr;
}

std::string OmegaWts::getParameterHeader(void) {

	 std::vector<double> w = omegaPrior[activeState]->getDndsWts();
	 std::string pStr = "";
	 for (int i=0; i<w.size(); i++)
		{
		char tempCh[50];
		sprintf(tempCh, "OmegaWt(%d)\t", i+1);
		std::string tempStr = tempCh;
		pStr += tempStr;
		}
	return pStr;
}


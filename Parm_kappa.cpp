#include <iomanip>
#include <iostream>
#include "MbRandom.h"
#include "Model.h"
#include "Parm_kappa.h"
#include "Parm_tree.h"
#include "Util.h"



TiTv::TiTv(MbRandom *rp, Model* mp) {

	ranPtr = rp;
	modelPtr = mp;
	k = ranPtr->exponentialRv(1.0) / ranPtr->exponentialRv(1.0);
}

TiTv::TiTv(TiTv &t) {

	clone(t);
}

TiTv& TiTv::operator=(const TiTv &t) {

	if (this != &t)
		clone(t);
	return *this;
}

double TiTv::change(void) {

	// update conditional likelihood and transition probability flags
	Topology* t = modelPtr->getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();

	// change the transition/transversion rate ratio
	double tuning = log(4.0);
	double oldK = k;
	double newK = oldK * exp(tuning*(ranPtr->uniformRv()-0.5));
	k = newK;
	return log(newK) - log(oldK);
}

double TiTv::lnProbability(void) {

	return -2.0 * log(1.0 + k);
}

void TiTv::print(void) {

	std::cout << "Kappa = " << std::fixed << std::setprecision(5) << k << std::endl;
}

void TiTv::clone(const TiTv &t) {

	ranPtr   = t.ranPtr;
	modelPtr = t.modelPtr;
	k        = t.k;
}

Kappa::Kappa(MbRandom *rp, std::string pn, Model* mp) : Parm(rp, pn, mp) {

	titvs[0] = new TiTv(ranPtr, modelPtr);
	titvs[1] = new TiTv( *titvs[0] );
}

Kappa::~Kappa(void) {

	delete titvs[0];
	delete titvs[1];
}

double Kappa::lnPriorRatio(void) {

	return titvs[activeState]->lnProbability() - titvs[Util::flip(activeState)]->lnProbability();
}

double Kappa::change(void) {

	numAttemptedChanges++;
	return titvs[activeState]->change();
}

void Kappa::print(void) {

	titvs[activeState]->print();
}

void Kappa::keep(void) {

	*titvs[Util::flip(activeState)] = *titvs[activeState];
}

void Kappa::restore(void) {

	*titvs[activeState] = *titvs[Util::flip(activeState)];
}

std::string Kappa::getParameterStr(void) {

	double k = titvs[activeState]->getRate();
	char tempCh[50];
	sprintf(tempCh, "%1.3lf\t", k);
	std::string pStr = tempCh;
	return pStr;
}

std::string Kappa::getParameterHeader(void) {

	std::string pStr = "Kappa\t";
	return pStr;
}


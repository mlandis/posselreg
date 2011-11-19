#include <iomanip>
#include <iostream>
#include "MbRandom.h"
#include "Model.h"
#include "Parm_lambda.h"
#include "Parm_tree.h"
#include "Util.h"



SwRate::SwRate(MbRandom *rp, Model* mp, double swLambda) {

	ranPtr = rp;
	modelPtr = mp;
	exponentialPriorParm = swLambda;
	lam = ranPtr->exponentialRv(exponentialPriorParm);
}

SwRate::SwRate(SwRate &r) {

	clone(r);
}

SwRate& SwRate::operator=(const SwRate &r) {

	if (this != &r)
		clone(r);
	return *this;
}

double SwRate::change(void) {

	// update conditional likelihood and transition probability flags
	Topology* t = modelPtr->getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();

	// change the switching rate
	double tuning = log(4.0);
	double oldL = lam;
	double newL = oldL * exp(tuning*(ranPtr->uniformRv()-0.5));
	lam = newL;
	return log(newL) - log(oldL);
}

double SwRate::lnProbability(void) {

	return ranPtr->lnExponentialPdf(exponentialPriorParm, lam);
}

void SwRate::print(void) {

	std::cout << "Lambda = " << std::fixed << std::setprecision(5) << lam << std::endl;
}

void SwRate::clone(const SwRate &r) {

	ranPtr   = r.ranPtr;
	modelPtr = r.modelPtr;
	lam      = r.lam;
	exponentialPriorParm = r.exponentialPriorParm;
}

Lambda::Lambda(MbRandom *rp, std::string pn, Model* mp, double swLambda) : Parm(rp, pn, mp) {

	swrate[0] = new SwRate(ranPtr, modelPtr, swLambda);
	swrate[1] = new SwRate( *swrate[0] );
}

Lambda::~Lambda(void) {

	delete swrate[0];
	delete swrate[1];
}

double Lambda::lnPriorRatio(void) {

	return swrate[activeState]->lnProbability() - swrate[Util::flip(activeState)]->lnProbability();
}

double Lambda::change(void) {

	numAttemptedChanges++;
	return swrate[activeState]->change();
}

void Lambda::print(void) {

	swrate[activeState]->print();
}

void Lambda::keep(void) {

	*swrate[Util::flip(activeState)] = *swrate[activeState];
}

void Lambda::restore(void) {

	*swrate[activeState] = *swrate[Util::flip(activeState)];
}

std::string Lambda::getParameterStr(void) {

	double r = swrate[activeState]->getRate();
	char tempCh[50];
	sprintf(tempCh, "%1.3lf\t", r);
	std::string pStr = tempCh;
	return pStr;
}

std::string Lambda::getParameterHeader(void) {

	std::string pStr = "Lambda\t";
	return pStr;
}


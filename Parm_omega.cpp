#include <algorithm>
#include <iomanip>
#include <iostream>
#include "MbRandom.h"
#include "Model.h"
#include "Parm_omega.h"
#include "Parm_tree.h"
#include "Util.h"



Dnds::Dnds(MbRandom *rp, int k, Model* mp) {

	ranPtr        = rp;
	modelPtr      = mp;
	numCategories = k;
	dndsVals.resize(numCategories);
	double dS = ranPtr->exponentialRv(1.0);
	for (int i=0; i<numCategories; i++)
		dndsVals[i] = ranPtr->exponentialRv(1.0) / dS;
	sort(dndsVals.begin(), dndsVals.end());
}

Dnds::Dnds(Dnds &d) {

	clone(d);
}

Dnds& Dnds::operator=(const Dnds &d) {

	if (this != &d)
		clone(d);
	return *this;
}

double Dnds::change(void) {

	// update conditional likelihood and transition probability flags
	Topology* t = modelPtr->getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();

	// set the maximum distance of the move
	double x = 0.1;
	
	// propose new values for all of the omegas
	bool isValid = true;
	std::vector<double> newVals(numCategories);
	std::vector<double> offsets(numCategories);
	do
		{
		// pick the values for each of the categories
		for (int i=0; i<numCategories; i++)
			offsets[i] = x*(ranPtr->uniformRv()-0.5);
			
		// calculate the new points
		for (int i=0; i<numCategories; i++)
			newVals[i] = dndsVals[i] + offsets[i];
		
		// decide if the new point is valid
		isValid = true;
		for (int i=0; i<numCategories; i++)
			{
			if (newVals[i] < 0.0)
				{
				isValid = false;
				break;
				}
			for (int j=i+1; j<numCategories; j++)
				{
				if (newVals[i] > newVals[j])
					{
					isValid = false;
					break;
					}
				}
			if (isValid == false)
				break;
			}
	
		} while (isValid == false);
	
	// commit the change
	for (int i=0; i<numCategories; i++)
		dndsVals[i] = newVals[i];
	
	return 0.0;
}

double Dnds::lnProbability(void) {

	// OK up to a constant, I think!
	double sum = 1.0;
	for (int i=0; i<numCategories; i++)
		sum += dndsVals[i];
	double x = -(numCategories + 1) * log(sum);
	return x;
}

void Dnds::print(void) {

	std::cout << "Omegas: ";
	for (int i=0; i<numCategories; i++)
		std::cout << std::fixed << std::setprecision(10) << dndsVals[i] << " ";
	std::cout << std::endl;
}

void Dnds::clone(const Dnds &d) {

	ranPtr = d.ranPtr;
	modelPtr = d.modelPtr;
	numCategories = d.numCategories;
	dndsVals.resize(numCategories);
	for (int i=0; i<numCategories; i++)
		dndsVals[i] = d.dndsVals[i];
}

Omega::Omega(MbRandom *rp, std::string pn, Model* mp, int k) : Parm(rp, pn, mp) {

	dnds[0] = new Dnds(ranPtr, k, modelPtr);
	dnds[1] = new Dnds( *dnds[0] );
}

Omega::~Omega(void) {

	delete dnds[0];
	delete dnds[1];
}

double Omega::lnPriorRatio(void) {

	return dnds[activeState]->lnProbability() - dnds[Util::flip(activeState)]->lnProbability();
}

double Omega::change(void) {

	numAttemptedChanges++;
	return dnds[activeState]->change();
}

void Omega::print(void) {

	dnds[activeState]->print();
}

void Omega::keep(void) {

	*dnds[Util::flip(activeState)] = *dnds[activeState];
}

void Omega::restore(void) {

	*dnds[activeState] = *dnds[Util::flip(activeState)];
}

std::string Omega::getParameterStr(void) {

	 std::vector<double> w = dnds[activeState]->getDnds();
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

std::string Omega::getParameterHeader(void) {

	 std::vector<double> w = dnds[activeState]->getDnds();
	 std::string pStr = "";
	 for (int i=0; i<w.size(); i++)
		{
		char tempCh[50];
		sprintf(tempCh, "Omega(%d)\t", i+1);
		std::string tempStr = tempCh;
		pStr += tempStr;
		}
	return pStr;
}



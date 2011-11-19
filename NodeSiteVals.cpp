#include <iomanip>
#include <iostream>
#include "NodeSiteVals.h"




NodeSiteVals::NodeSiteVals(void) {

}

NodeSiteVals::~NodeSiteVals(void) {

}

void NodeSiteVals::addSample(int s, double p, double d) {

	sampleNum.push_back( s );
	probPos.push_back( p );
	aveDnds.push_back( d );
}

void NodeSiteVals::print(void) {

	for (int i=0; i<probPos.size(); i++)
		std::cout << std::fixed << std::setprecision(4) << probPos[i] << " ";
	std::cout << std::endl;
	for (int i=0; i<aveDnds.size(); i++)
		std::cout << std::fixed << std::setprecision(4) << aveDnds[i] << " ";
	std::cout << std::endl;
}

double NodeSiteVals::averageDnds(int discardSamplesUpTo) {

	int beginningSample = sampleNum.size() - 1;
	for (int i=0; i<sampleNum.size(); i++)
		{
		if ( sampleNum[i] > discardSamplesUpTo )
			{
			beginningSample = i;
			break;
			}
		}

	double sum = 0.0;
	int n = 0;
	for (int i=beginningSample; i<aveDnds.size(); i++)
		{
		sum += aveDnds[i];
		n++;
		}
		
	return (sum / n);
}

double NodeSiteVals::averagePosProb(int discardSamplesUpTo) {

	int beginningSample = sampleNum.size() - 1;
	for (int i=0; i<sampleNum.size(); i++)
		{
		if ( sampleNum[i] > discardSamplesUpTo )
			{
			beginningSample = i;
			break;
			}
		}

	double sum = 0.0;
	int n = 0;
	for (int i=beginningSample; i<aveDnds.size(); i++)
		{
		sum += probPos[i];
		n++;
		}
		
	return (sum / n);
}

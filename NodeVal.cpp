#include <iomanip>
#include <iostream>
#include "NodeVal.h"




NodeVals::NodeVals(void) {

	count = 0;
}

NodeVals::~NodeVals(void) {

}

void NodeVals::addSample(int s, double v) {

	count++;
	sampleNum.push_back( s );
	brlens.push_back( v );
}

void NodeVals::print(void) {

	std::cout << "Count = " << count << std::endl;
	for (int i=0; i<brlens.size(); i++)
		std::cout << std::fixed << std::setprecision(4) << brlens[i] << " ";
	std::cout << std::endl;
}

double NodeVals::averageBrlen(int discardSamplesUpTo) {

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
	for (int i=beginningSample; i<brlens.size(); i++)
		{
		sum += brlens[i];
		n++;
		}
		
	return (sum / n);
}

int NodeVals::numSamplesWithThisBipartition(int discardSamplesUpTo) {

	int num = 0;
	for (int i=0; i<sampleNum.size(); i++)
		if ( sampleNum[i] > discardSamplesUpTo )
			num++;
	return num;
}

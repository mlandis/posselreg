#include "NodeSiteKey.h"
#include "MbBitfield.h"



NodeSiteKey::NodeSiteKey(MbBitfield &bf, int s) {

	site = s;
	partDescription = new MbBitfield( bf.dim() );
	partDescription->inject(bf);
}

NodeSiteKey::~NodeSiteKey(void) {

}

bool NodeSiteKey::operator<(const NodeSiteKey &a) const {

	//std::cout << "comparing: " << *partDescription << " (" << site << ") - " << *(a.partDescription) << " (" << a.site << ")" << std::endl;

	if (a.site < site)
		{
		return true;
		}
	else if (a.site == site)
		{
		MbBitfield* myPrt = partDescription;
		MbBitfield* aPrt  = a.partDescription;
		for (int i=0; i<myPrt->dim(); i++)
			{
			if ( myPrt->isBitSet(i) == true && aPrt->isBitSet(i) == false )
				return true;
			else if ( myPrt->isBitSet(i) == false && aPrt->isBitSet(i) == true )
				return false;
			}
		}
	return false;
}

void NodeSiteKey::print(void) {

	std::cout << *partDescription << " " << site;
}


#include "NodeKey.h"
#include "MbBitfield.h"



NodeKey::NodeKey(MbBitfield &bf) {

	partDescription = new MbBitfield( bf.dim() );
	partDescription->inject(bf);
}

NodeKey::~NodeKey(void) {

}

bool NodeKey::operator<(const NodeKey &a) const {

	MbBitfield* myPrt = partDescription;
	MbBitfield* aPrt  = a.partDescription;
	for (int i=0; i<myPrt->dim(); i++)
		{
		if ( myPrt->isBitSet(i) == true && aPrt->isBitSet(i) == false )
			return true;
		else if ( myPrt->isBitSet(i) == false && aPrt->isBitSet(i) == true )
			return false;
		}
	return false;
	/*if ( *(a.partDescription) < *partDescription )
		{
		return true;
		}
	return false;*/
}

void NodeKey::print(void) {

	std::cout << *partDescription;
}


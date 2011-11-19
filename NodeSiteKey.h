#ifndef NODESITEKEY_H
#define NODESITEKEY_H

#include "MbBitfield.h"

class NodeSiteKey {

	public:
						   NodeSiteKey(MbBitfield &bf, int s);
						   ~NodeSiteKey(void);
					bool   operator<(const NodeSiteKey &a) const;
					void   print(void);
			 MbBitfield*   getPartition(void) { return partDescription; }

	private:
	         MbBitfield*   partDescription;
			         int   site;
};


#endif
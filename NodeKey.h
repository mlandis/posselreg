#ifndef NODEKEY_H
#define NODEKEY_H

#include "MbBitfield.h"

class NodeKey {

	public:
						   NodeKey(MbBitfield &bf);
						   ~NodeKey(void);
					bool   operator<(const NodeKey &a) const;
					void   print(void);
			 MbBitfield*   getPartition(void) { return partDescription; }

	private:
	         MbBitfield*   partDescription;
};


#endif
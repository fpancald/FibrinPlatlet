#ifndef WLC_FORCE_H_
#define WLC_FORCE_H_

#include "SystemStructures.h"

void WLC_Force(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams);
	
void GetStrainParameters(NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,  
	GeneralParams& generalParams,
	DomainParams& domainParams);




#endif

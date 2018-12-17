#ifndef PLT_VOL_EXC_FORCE_H_
#define PLT_VOL_EXC_FORCE_H_

#include <vector>
#include "SystemStructures.h"

void Plt_Vol_Exc_Force(
    NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams,
    PltInfoVecs& pltInfoVecs,
    AuxVecs& auxVecs);

#endif
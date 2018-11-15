#ifndef PLTVLMPUSHDEVICE_H_
#define PLTVLMPUSHDEVICE_H_

#include <vector>
#include "SystemStructures.h"

void PltVlmPushOnDevice(
  NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams,
  PltInfoVecs& pltInfoVecs,
  AuxVecs& auxVecs);

#endif
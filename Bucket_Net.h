#ifndef BUCKET_NET_H_
#define BUCKET_NET_H_

#include "SystemStructures.h"

void init_dim_general(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);

void init_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);


void build_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);


void extend_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);



#endif

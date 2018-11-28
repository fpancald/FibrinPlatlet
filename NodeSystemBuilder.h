
/*
 * NodeSystemBuilder.h
 *
 *  Created on: 25 авг. 2014 г.
 *      Author: yan
 */

#ifndef NODESYSTEMBUILDER_H_
#define NODESYSTEMBUILDER_H_

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "glm/glm.hpp"
#include <list>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

class NodeSystemDevice;

class NodeSystemBuilder {
public:
	NodeSystemBuilder(double _epsilon, double _dt);
	~NodeSystemBuilder();

	//the reason we need to have each buildNode contain previous and next node id's is
	//that we will use this structure to set preferred edge angle to incorporate bending.
	struct BuildNode {

		std::vector<unsigned> next;
		std::vector<unsigned> prev;
		unsigned id;

		BuildNode(unsigned id_) : id(id_) {}	 //the only thing set is the main id
	};
	//platelet build struct
	struct BuildPlt {

		unsigned id;

		BuildPlt(unsigned id_) : id(id_) {}	 //the only thing set is the main id
	};



public:

	double epsilon, dt;

	unsigned originNodeCount = 0;
	unsigned originLinkCount;
	unsigned numEdges;
	unsigned numNodes;
	unsigned torsionCount = 0;
	double defaultMass = 10;
	double defaultResistance = 3.769911184308;
	double defaultSpringStiffness = 1.0;
	double defaultTorsionSpringStiffness = 1.0;
	double defaultPersistanceLength = 1.0;
	double defaultTemperature = 300.0; // 300' kelvin ~ 27' celsius
	double defaultContourLengthMultiplier = 1.0;
	double defaultBoltzmannConstant = 1.3806488e-8;
	double defaultUnitsPerExtraNode = 1.0;
	double defaultLinkDiameter = 0.1; //used for fiber-fiber linking threshold distance.
	double compressionPercent = 1;
	unsigned defaultExtraNodesPerEdge = 0;

	bool linking = true;
	bool useConstantNumberOfExtraNodes = true;
	bool wormLikeEnabled = false;

	double defaultSubNodeMass = 1.0;
	bool useExtraNodes = false;

	//platelets
	unsigned originPltCount = 0;
	unsigned numPlts;
	double defaultPltMass = 100;
	double pltRForce = 1.4;
	double pltR = 0.1;
	unsigned plt_other_intrct = 100;
	unsigned plt_tndrl_intrct = 5;
	double pltForce = 1.0;
	double pltDensity = 0.0;
	bool pltfrcfld = false;
	bool plttndrl = true;
	bool pltonplt = false;

	//std::shared_ptr<NodeSystemDevice> host_ptr_devNodeSystem;
	std::vector<std::shared_ptr<BuildNode>> buildNodes;
	//platelets
	std::vector<std::shared_ptr<BuildPlt>> buildPlts;

	thrust::host_vector<unsigned> hostNodeIds;
	//platelets
	thrust::host_vector<unsigned> hostPltIds;

	thrust::host_vector<unsigned> hostWLCEdgeLeft;
	thrust::host_vector<unsigned> hostWLCEdgeRight;
	thrust::host_vector<double> hostWLCLenZero;

	//only used if subnodes are used.

	thrust::host_vector<unsigned> hostSpringDivisionCount;
	thrust::host_vector<unsigned> hostWLCEdgeSubLeft;
	thrust::host_vector<unsigned> hostWLCEdgeSubRight;
	thrust::host_vector<double> hostWLCLenSubZero;

	thrust::host_vector<unsigned> hostTorsionEdgeLeft;
	thrust::host_vector<unsigned> hostTorsionEdgeCenter;
	thrust::host_vector<unsigned> hostTorsionEdgeRight;
	thrust::host_vector<double> hostTorsionAngleZero;

	thrust::host_vector<double> hostPosX;
	thrust::host_vector<double> hostPosY;
	thrust::host_vector<double> hostPosZ;
	//platelets
	thrust::host_vector<double> hostPltPosX;
	thrust::host_vector<double> hostPltPosY;
	thrust::host_vector<double> hostPltPosZ;



	thrust::host_vector<bool> hostIsNodeFixed;
	//platelets
	thrust::host_vector<bool> hostIsPltFixed;

	std::vector<glm::dvec3> nodePositions;
	//platelets
	std::vector<glm::dvec3> pltPositions;

	unsigned addNode(__attribute__ ((unused)) double , glm::dvec3 );
	//platelets
	unsigned addPlt(__attribute__ ((unused)) double , glm::dvec3 );

	void putLinearSpring(__attribute__ ((unused)) unsigned,__attribute__ ((unused)) unsigned);
	void putWormLikeSpring(unsigned, unsigned);
	void putSpring(unsigned, unsigned);
	void putSubSpring(unsigned, unsigned);//updates edges if using subnodes
	void putTorsionSpring(unsigned, unsigned, unsigned);
	void addSubnodes(void);
	std::list<glm::dvec3> fillSpace(glm::dvec3, glm::dvec3, unsigned);
	void generateBuildNodesTriplets(void);
	void fixNode(unsigned);
	//platelets
	void fixPlt(unsigned);
	//void setSystemForParallelComputation();
	std::shared_ptr<NodeSystemDevice> create();


};

#endif /* NODESYSTEMBUILDER_H_ */

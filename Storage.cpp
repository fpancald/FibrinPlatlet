

#include "System.h"
#include "System_Builder.h"
#include "SystemStructures.h"
#include "Storage.h"


Storage::Storage(std::weak_ptr<System> a_system,
	std::weak_ptr<SystemBuilder> b_system , __attribute__ ((unused)) const std::string& a_fileName) {
	//std::cout << "FDM constructor" << std::endl;

	system = a_system;
	builder = b_system;

};



void Storage::save_params(void) {
	std::shared_ptr<System> sys = system.lock();
	if (sys) {

		//first create a new file using the current network strain
		
		std::string format = ".sta";
		
		std::string strain =  std::to_string(sys->generalParams.currentTime);
		std::string initial = "Params/Param_";
		std::ofstream ofs;
		std::string Filename = initial + strain + format;
		ofs.open(Filename.c_str());



		//unsigned maxNeighborCount = sys->generalParams.maxNeighborCount;
		unsigned maxNodeCount = sys->generalParams.maxNodeCount;
		unsigned originalNodeCount = sys->generalParams.originNodeCount;
		unsigned originalEdgeCount = sys->generalParams.originLinkCount;
		unsigned edgeCountDiscretize = sys->generalParams.originEdgeCount;
		//Now first place strain
		ofs << std::setprecision(5) <<std::fixed<< "time " << sys->generalParams.currentTime<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minX " << sys->domainParams.minX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxX " << sys->domainParams.maxX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minY " << sys->domainParams.minY<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxY " << sys->domainParams.maxY<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minZ " << sys->domainParams.minX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxZ " << sys->domainParams.maxX<<std::endl;
		
		
		ofs << std::setprecision(5) <<std::fixed<< "original_node_count " << originalNodeCount <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "node_count_discretize " << maxNodeCount <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "original_edge_count " << originalEdgeCount <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "edge_count_discretize " << edgeCountDiscretize <<std::endl;
		
		//place nodes
		for (unsigned i = 0; i < sys->nodeInfoVecs.nodeLocX.size(); i++) {
			double x = sys->nodeInfoVecs.nodeLocX[i];
			double y = sys->nodeInfoVecs.nodeLocY[i];
			double z = sys->nodeInfoVecs.nodeLocZ[i];
			ofs << std::setprecision(5) <<std::fixed<< "node " << x << " " << y << " " << z <<std::endl;
		
		}
		
		//place plts
		for (unsigned i = 0; i < sys->pltInfoVecs.pltLocX.size(); i++) {
			double x = sys->pltInfoVecs.pltLocX[i];
			double y = sys->pltInfoVecs.pltLocY[i];
			double z = sys->pltInfoVecs.pltLocZ[i];
			ofs << std::setprecision(5) <<std::fixed<< "plt " << x << " " << y << " " << z <<std::endl;
		
		}
		//place force node is experiencing
		for (unsigned i = 0; i < sys->nodeInfoVecs.nodeLocX.size(); i++) {
			ofs << std::setprecision(5) <<std::fixed<< "force_on_node " << sys->nodeInfoVecs.sumForcesOnNode[i]<<std::endl;
		
		}

		//place original edges
		for (unsigned edge = 0; edge < sys->generalParams.originEdgeCount; edge++) {
			unsigned idL = sys->nodeInfoVecs.deviceEdgeLeft[edge];
			unsigned idR = sys->nodeInfoVecs.deviceEdgeRight[edge];
			ofs <<"original_edge_discretized " <<idL <<" "<< idR <<std::endl;
			
		}
				 
		//place added edges
		for (unsigned edge = sys->generalParams.originEdgeCount; edge < sys->generalParams.currentEdgeCount; edge++) {
			unsigned idL = sys->nodeInfoVecs.deviceEdgeLeft[edge];
			unsigned idR = sys->nodeInfoVecs.deviceEdgeRight[edge];
			ofs <<"added_edge " <<idL <<" "<< idR <<std::endl;
			
		}

		//original edge strain
		for (unsigned i = 0; i < sys->generalParams.originEdgeCount; i++ ){
			double val = sys->nodeInfoVecs.discretizedEdgeStrain[i];

			ofs << std::setprecision(5)<< std::fixed<<"original_edge_strain " << val <<std::endl;
		}
				
		//original edge alignment
		for (unsigned i = 0; i < sys->generalParams.originEdgeCount; i++ ){
			double val = sys->nodeInfoVecs.discretizedEdgeAlignment[i];
			ofs << std::setprecision(5)<< std::fixed<<"original_edge_alignment " << val <<std::endl;
		}

		//added edge strain
		for (unsigned i = sys->generalParams.originEdgeCount; i < sys->generalParams.currentEdgeCount; i++ ){
			double val = sys->nodeInfoVecs.discretizedEdgeStrain[i];
			ofs << std::setprecision(5)<< std::fixed<<"added_edge_strain " << val <<std::endl;
		}
		
		//added links per node.
		for (unsigned i = 0; i < sys->generalParams.maxNodeCount; i++ ){
			unsigned val = sys->wlcInfoVecs.currentNodeEdgeCountVector[i] - 
				sys->wlcInfoVecs.numOriginalNeighborsNodeVector[i];
			ofs << std::setprecision(5)<< std::fixed<<"bind_sites_per_node " << val <<std::endl;
		}



	}
};


void Storage::print_VTK_File() {

	std::shared_ptr<System> sys = system.lock();
	if (sys) {
		iteration+=1;
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = "AnimationTest/FibrinNetwork_";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());


		unsigned maxNodeCount = sys->generalParams.maxNodeCount;
		__attribute__ ((unused)) unsigned maxNeighborCount = (sys->generalParams).maxNeighborCount;

		unsigned numEdges = sys->generalParams.currentEdgeCount;

		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;


		ofs << "POINTS " << maxNodeCount << " float" << std::endl;
		for (unsigned i = 0; i< maxNodeCount; i++) { 
			double xPos = sys->nodeInfoVecs.nodeLocX[i];
			double yPos = sys->nodeInfoVecs.nodeLocY[i];
			double zPos = sys->nodeInfoVecs.nodeLocZ[i];

			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}
		//now plot particles


		unsigned numCells = numEdges;
		unsigned numNumsInCells = 3 * numEdges;


		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;

		for (unsigned edge = 0; edge < numEdges; edge++) {
			unsigned idL = sys->nodeInfoVecs.deviceEdgeLeft[edge];
			unsigned idR = sys->nodeInfoVecs.deviceEdgeRight[edge];

			ofs<< 2 << " " << idL << " " << idR << std::endl;
		}
	/*	for (unsigned idA = 0; idA < maxNodeCount; idA++ ){

			unsigned beginIndex = idA * maxNeighborCount;
			unsigned endIndex = beginIndex + maxNeighborCount;
			for (unsigned i = beginIndex; i < endIndex; i++) {//currentSpringCount is the length of index and value vectors
				unsigned idB = sys->wlcInfoVecs.globalNeighbors[i];//look through possible neighbors. May contain ULONG_MAX

				//counts only half
				if ((idA < idB) && (idB < maxNodeCount) ) {
					ofs<< 2 << " " << idA << " " << idB << std::endl;
				}
			}
		}*/

		ofs << "CELL_TYPES " << numCells << std::endl;
		for (unsigned i = 0; i< numEdges; i++) {
			ofs << 3 << std::endl; //edge joining two points
		}



		//
		ofs << "CELL_DATA " << numCells << std::endl;
		ofs << "SCALARS Fiber_Strain double " << std::endl;
		ofs << "LOOKUP_TABLE default "  << std::endl;
		for (unsigned edge = 0; edge < numEdges; edge++) {
			unsigned idA = sys->nodeInfoVecs.deviceEdgeLeft[edge];
			unsigned idB = sys->nodeInfoVecs.deviceEdgeRight[edge];

			unsigned begin = idA * sys->generalParams.maxNeighborCount;
			unsigned end = begin + sys->generalParams.maxNeighborCount;
			double L0;
			for (unsigned i = begin; i < end; i++) {
				unsigned idTemp = sys->wlcInfoVecs.globalNeighbors[i];
				if (idTemp == idB){
					L0 = sys->wlcInfoVecs.lengthZero[i];
				}
			}
			double xL = sys->nodeInfoVecs.nodeLocX[idA];
			double yL = sys->nodeInfoVecs.nodeLocY[idA];
			double zL = sys->nodeInfoVecs.nodeLocZ[idA];
			double xR = sys->nodeInfoVecs.nodeLocX[idB];
			double yR = sys->nodeInfoVecs.nodeLocY[idB];
			double zR = sys->nodeInfoVecs.nodeLocZ[idB];

			double L1 = std::sqrt( (xL - xR)*(xL - xR)+(yL - yR)*(yL - yR)+(zL - zR)*(zL - zR));
			double strain = (L1 - L0) / L0;
			ofs << std::fixed << strain   << std::endl;

		}
	/*	for (unsigned idA = 0; idA < maxNodeCount; idA++ ){

			unsigned beginIndex = idA * maxNeighborCount;
			unsigned endIndex = beginIndex + maxNeighborCount;
			for (unsigned i = beginIndex; i < endIndex; i++) {//currentSpringCount is the length of index and value vectors
				unsigned idB = sys->wlcInfoVecs.globalNeighbors[i];//look through possible neighbors. May contain ULONG_MAX

				if ((idA < idB) && (idB < maxNodeCount) ) {
					__attribute__ ((unused)) unsigned index = idA * maxNeighborCount + idB;
					double L0 = sys->wlcInfoVecs.lengthZero[i];
					double xL = sys->nodeInfoVecs.nodeLocX[idA];
					double yL = sys->nodeInfoVecs.nodeLocY[idA];
					double zL = sys->nodeInfoVecs.nodeLocZ[idA];
					double xR = sys->nodeInfoVecs.nodeLocX[idB];
					double yR = sys->nodeInfoVecs.nodeLocY[idB];
					double zR = sys->nodeInfoVecs.nodeLocZ[idB];



					double L1 = std::sqrt( (xL - xR)*(xL - xR)+(yL - yR)*(yL - yR)+(zL - zR)*(zL - zR));
					double strain = (L1 - L0) / L0;
					ofs << std::fixed << strain   << std::endl;
				}
			}
		}*/

		ofs.close();

	}

	//now print platelets 
	if ((sys)) {
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = "AnimationTest/Platelet";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());
		
	
		unsigned maxPltCount = sys->generalParams.maxPltCount;
		
		double xPos;
		double yPos;
		double zPos;
		
		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
		
		 
		ofs << "POINTS " << maxPltCount  << " float" << std::endl;
		for (unsigned i = 0; i< maxPltCount; i++) {
			xPos = sys->pltInfoVecs.pltLocX[i];
			yPos = sys->pltInfoVecs.pltLocY[i];
			zPos = sys->pltInfoVecs.pltLocZ[i];

			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}

		//std::cout<<'here1'<<std::flush;
		
		unsigned numCells = 1;

		unsigned numNumsInCells = 1 + maxPltCount;

		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;
		
		//place edges as cells of type 2. 
		ofs<< maxPltCount;
		for (unsigned point = 0; point < maxPltCount; point++ ){
			ofs<< " " << point;
		}
		ofs<<" "<< std::endl;

		ofs << "CELL_TYPES " << numCells << std::endl;  
		//set edges and last set scattered points
				
		ofs << 2 << std::endl;//scatter points for capsid
		
	}
	
	//now print platelets with attatchments
	if ((sys)) {
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = "AnimationTest/PlateletConn";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());
		
	
		unsigned maxPltCount = sys->generalParams.maxPltCount;

		unsigned num_connections = sys->pltInfoVecs.numConnections;
		
		double xPos;
		double yPos;
		double zPos;
		
		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
		
		 
		ofs << "POINTS " << maxPltCount + num_connections << " float" << std::endl;
		for (unsigned i = 0; i< maxPltCount; i++) {
			xPos = sys->pltInfoVecs.pltLocX[i];
			yPos = sys->pltInfoVecs.pltLocY[i];
			zPos = sys->pltInfoVecs.pltLocZ[i];

			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}

		//set location for nodes that plt is connected to
		//ie  
		for (unsigned i = 0; i < num_connections; i++ ) {
			unsigned node_id = sys->pltInfoVecs.nodeImagingConnection[i];

			xPos = sys->nodeInfoVecs.nodeLocX[node_id];
			yPos = sys->nodeInfoVecs.nodeLocY[node_id];
			zPos = sys->nodeInfoVecs.nodeLocZ[node_id];
			
			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		
		}


		//std::cout<<'here1'<<std::flush;
		
		unsigned numCells = 1;
		numCells += num_connections;//add conections cells for edges

		unsigned numNumsInCells = 1 + maxPltCount;
		numNumsInCells += 3 * num_connections;//3 numbers per edge

		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;
		
		//place edges as cells of type 2. 
		ofs<< maxPltCount;
		for (unsigned point = 0; point < maxPltCount; point++ ){
			ofs<< " " << point;
		}
		ofs<<" "<< std::endl;

		//std::cout<<'here2'<<std::flush;
		for (unsigned edge = 0; edge < num_connections; edge++ ){

			unsigned node_Id = maxPltCount + edge;
			unsigned plt_id = sys->pltInfoVecs.pltImagingConnection[edge];
				
			ofs <<2<< " "<< node_Id << " "<< plt_id <<std::endl;
		}
		ofs << "CELL_TYPES " << numCells << std::endl;  
		//set edges and last set scattered points
				
		ofs << 2 << std::endl;//scatter points for capsid
		
	//	std::cout<<'here3'<<std::flush;
		for (unsigned edge = 0; edge< num_connections; edge++ ){
			ofs<< 3 <<std::endl;
		}
		ofs.close();
	}

};

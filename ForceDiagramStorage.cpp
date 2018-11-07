

#include "NodeSystemDevice.h"
#include "NodeSystemBuilder.h"
#include "SystemStructures.h"
#include "ForceDiagramStorage.h"


ForceDiagramStorage::ForceDiagramStorage(std::weak_ptr<NodeSystemDevice> a_system,
	std::weak_ptr<NodeSystemBuilder> b_system , const std::string& a_fileName) {
	//std::cout << "FDM constructor" << std::endl;

	system = a_system;
	builder = b_system;
/*	bn = a_fileName; //this will be used later to open files
	std::ofstream statesOutput(a_fileName + ".sta");
	std::ofstream statesOutputStrain(a_fileName + "_Strain.sta");

	std::shared_ptr<NodeSystemDevice> sysA = system.lock();
	std::shared_ptr<NodeSystemBuilder> sysB = builder.lock();

	if ((sysA) && (sysB) ){
		unsigned maxNodeCount = sysA->generalParams.maxNodeCount;
		__attribute__ ((unused)) unsigned maxNeighborCount = sysA->generalParams.maxNeighborCount;

		statesOutput << "node_count " << maxNodeCount << '\n';
		statesOutput << "origin_node_count " << sysA->generalParams.originNodeCount << '\n';
		statesOutput << "origin_link_count " << sysA->generalParams.originLinkCount << '\n';
		statesOutput << "sub_node_count " << sysA->generalParams.subNodeCount << std::endl;//system->getSubNodesSize() << '\n';
		statesOutput << "link_count " << sysA->generalParams.originEdgeCount << '\n';

		for (unsigned edge = 0; edge < sysB->hostWLCEdgeLeft.size(); edge++) {
			unsigned idLeft = sysB->hostWLCEdgeLeft[edge];
			unsigned idRight = sysB->hostWLCEdgeRight[edge];
			statesOutput << '\n' << idLeft << ' ' << idRight;
		}

	}


	statesOutput.close();*/
};

void ForceDiagramStorage::updateStrain() {

};

void ForceDiagramStorage::updateTotalStrain(void) {
	/*std::shared_ptr<NodeSystemDevice> sys = system.lock();
	if (sys) {


		std::string format = ".sta";
		std::string strain =  std::to_string(currentStrain);
		std::string initial = "StrainTest/Strain_";
		std::ofstream ofs;
		std::string Filename = initial + format;
		ofs.open(Filename.c_str());



		unsigned maxNeighborCount = sys->generalParams.maxNeighborCount;
		unsigned maxNodeCount = sys->generalParams.maxNodeCount;
		unsigned originalNodeCount = sys->generalParams.originNodeCount;
		unsigned originalEdgeCount = sys->generalParams.originLinkCount;
		unsigned edgeCountDiscretize = sys->generalParams.originEdgeCount;
		//Now first place strain
		ofs << std::setprecision(5) <<std::fixed<< "network_strain " << currentStrain<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minX " << sys->domainParams.minX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxX " << sys->domainParams.maxX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minY " << sys->domainParams.minY<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxY " << sys->domainParams.maxY<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minZ " << sys->domainParams.minX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxZ " << sys->domainParams.maxX<<std::endl;




	}*/
}


void ForceDiagramStorage::print_VTK_File() {

	std::shared_ptr<NodeSystemDevice> sys = system.lock();
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
		unsigned maxNeighborCount = (sys->generalParams).maxNeighborCount;

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
		/*for (unsigned idA = 0; idA < maxNodeCount; idA++ ){

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
		ofs << "SCALARS magnitude double " << std::endl;
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
		/*for (unsigned idA = 0; idA < maxNodeCount; idA++ ){

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
	if (sys) {
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = "AnimationTest/Platelet_";
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


		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;


		ofs << "POINTS " << maxPltCount << " float" << std::endl;
		for (unsigned i = 0; i< maxPltCount; i++) {
			double xPos = sys->pltInfoVecs.pltLocX[i];
			double yPos = sys->pltInfoVecs.pltLocY[i];
			double zPos = sys->pltInfoVecs.pltLocZ[i];

			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}
		//now plot particles


		unsigned numCells = 1;
		unsigned numNumsInCells = maxPltCount+1;


		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;
		//place edges as cells of type 2. 
		ofs<< maxPltCount << " ";
		for (unsigned point = 0; point < maxPltCount; point++ ){
			ofs<< " " << point;
		}
		ofs<<" "<< std::endl;

		ofs << "CELL_TYPES " << numCells << std::endl;
		ofs << 2 << std::endl;//scatter points for capsid
		


		ofs.close();

	}
};

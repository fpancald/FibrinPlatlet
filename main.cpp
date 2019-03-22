

#include <iomanip>
#include <string>
#include <memory>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <inttypes.h>
#include <cstddef>

#include "pugixml/include/pugixml.hpp"

#include "System.h"

#include "System_Builder.h"
#include "Storage.h"



std::shared_ptr<System> createSystem(const char* schemeFile, std::shared_ptr<SystemBuilder> builder)	{
	pugi::xml_document doc;
	pugi::xml_parse_result parseResult = doc.load_file(schemeFile);

	if (!parseResult) {
		std::cout << "parse error in createNodeSystem: " << parseResult.description() << std::endl;
		return nullptr;
	}
	pugi::xml_node root = doc.child("data");
	pugi::xml_node nodes = root.child("nodes");
	pugi::xml_node plts = root.child("plts");
	pugi::xml_node links = root.child("links");
	pugi::xml_node props = root.child("settings");

	//first, we'll input settings
	if (!(root && nodes && links)) {
		std::cout << "couldn't find nessesary data\n";
		//return false;
	}

	//default settings in NSB.h
	if (auto p = props.child("resistance_fibrin"))
		builder->viscousDamp_Fibrin = (p.text().as_double());

	if (auto p = props.child("resistance_plt"))
		builder->viscousDamp_Plt = (p.text().as_double());

	if (auto p = props.child("spring-stiffness"))
		builder->defaultSpringStiffness = (p.text().as_double());

	if (auto p = props.child("torsion-stiffness"))
		builder->defaultTorsionSpringStiffness = (p.text().as_double());

	if (auto p = props.child("persistance-length"))
		builder->defaultPersistanceLength = (p.text().as_double());

	if (auto p = props.child("absolute-temperature"))
		builder->defaultTemperature = (p.text().as_double());

	if (auto p = props.child("contour-length-multiplier"))
		builder->defaultContourLengthMultiplier = (p.text().as_double());

	if (auto p = props.child("units-per-extra-node"))
		builder->defaultUnitsPerExtraNode = (p.text().as_double());

	if (auto p = props.child("extra-nodes-per-edge"))
		builder->defaultExtraNodesPerEdge = (p.text().as_uint());

	if (auto p = props.child("use-extra-nodes"))
		builder->useExtraNodes = (p.text().as_bool());

	if (auto p = props.child("constant-extra-nodes"))
		builder->useConstantNumberOfExtraNodes = (p.text().as_bool());

	if (auto p = props.child("link-diameter"))
		builder->defaultLinkDiameter = (p.text().as_double());

	if (auto p = props.child("worm-like"))
		builder->wormLikeEnabled = (p.text().as_bool());

	if (auto p = props.child("use-linking"))
		builder->linking = (p.text().as_bool());
////////////////////////////////////////////////////
//platelets parameters
	if (auto p = props.child("plt_mass"))
		builder->defaultPltMass = (p.text().as_double());

	if (auto p = props.child("plt_force"))
		builder->pltForce = (p.text().as_double());

	if (auto p = props.child("plt_other_intrct"))
		builder->plt_other_intrct = (p.text().as_uint());

	if (auto p = props.child("plt_tndrl_intrct"))
		builder->plt_tndrl_intrct = (p.text().as_uint());

	if (auto p = props.child("plt_r"))
		builder->pltR = (p.text().as_double());

	if (auto p = props.child("plt_r_force")) {
		builder->pltRForce = (p.text().as_double());
	}
	if (auto p = props.child("plt_r_adhesion")) {
		double RAdhesion=(p.text().as_double());
		if (RAdhesion>0.0 && RAdhesion<1.0){
			builder->pltRAdhesion = RAdhesion;
		}
		else{
			std::cout << "parse error: platelet adhesion radius fraction is not valid\n";
			return 0;
		}
	}


	if (auto p = props.child("plt_density")) {
		builder->pltDensity = (p.text().as_double());
		std::cout << "setting density: " << builder->pltDensity << std::endl;
	}

	if (auto p = props.child("use-pltforcefield")){
		builder->pltfrcfld = (p.text().as_bool());
		std::cout<<"frcFld: "<< builder->pltfrcfld<<std::endl;
	}
	if (auto p = props.child("use-plttndrl")){
		builder->plttndrl = (p.text().as_bool());	
		std::cout<<"plttndrl: "<< builder->plttndrl<<std::endl;
	}
	if (auto p = props.child("use-pltrelease")){
		builder->pltrelease = (p.text().as_bool());	
		std::cout<<"pltrelease: "<< builder->pltrelease<<std::endl;
	}
	if (auto p = props.child("use-plthandhand")){
		builder->plthandhand = (p.text().as_bool());	
		std::cout<<"plthandhand: "<< builder->plthandhand<<std::endl;
	}
	
	if (auto p = props.child("use-pltonplt")){
		builder->pltonplt = (p.text().as_bool());
		std::cout<<"plt_interact: "<< builder->pltonplt<<std::endl;
	}

	std::cout << "builder ptr address: " << builder << std::endl;
	std::vector<unsigned> originNodes;
//buid nodes
	double mass;
	double x, y, z; //variables to be used reading in data.
	double defaultMass = nodes.attribute("default-mass").as_double(-1.0);
	builder->defaultMass = defaultMass;

	//debuggin printl
	//std::cout<<"default mass: " << defaultMass << std::endl;
	for (auto node = nodes.child("node"); node; node = node.next_sibling("node")) {
		mass = node.attribute("mass").as_double(defaultMass);
		if (mass < 0.0) {
			std::cout << "parse error: node mass is undefined\n";
			return 0;
		}
		const char* text = node.text().as_string();

		//std::cout<<"mass: " << mass << std::endl;
		if (3 != sscanf(text, "%lf %lf %lf", &x, &y, &z)) {
			std::cout << "parse node error\n";
			return 0;
		}
		__attribute__ ((unused)) int unused = builder->addNode(mass, glm::dvec3(x, y, z));
	}

	unsigned from, to;
	for (auto link = links.child("link"); link; link = link.next_sibling("link")) {
		if (2 != sscanf(link.text().as_string(""), "%u %u" , &from, &to)) {
			std::cout << "parse link error\n";
			return 0;
		}
		builder->putSpring(from, to); //adds edges into saved vectors

	}

	std::cout << "post springs" << std::endl;


	//replace false entries with true if node is fixed.
	pugi::xml_node fixedRoot = root.child("fixed");
	if (fixedRoot) {
		for (auto node = fixedRoot.child("node"); node; node = node.next_sibling("node"))
			builder->fixNode(node.text().as_uint());
	}

	std::cout << "post fixed" << std::endl;

//build platelets

	__attribute__ ((unused)) double pltmass;
	double pltx, plty, pltz; //variables to be used reading in data.
	
	//only use platelet input if density is zero
	for (auto plt = plts.child("plt"); plt; plt = plt.next_sibling("plt")) {

		const char* text = plt.text().as_string();

		//std::cout<<"mass: " << mass << std::endl;
		if (3 != sscanf(text, "%lf %lf %lf", &pltx, &plty, &pltz)) {
			std::cout << "parse plt error\n";
			return 0;
		}
		__attribute__ ((unused)) int unused = builder->addPlt(builder->defaultPltMass, glm::dvec3(pltx, plty, pltz));
	}
	



	//replace false entries with true if node is fixed.
	pugi::xml_node fixedPltRoot = root.child("pltfixed");
	if (fixedPltRoot) {
		for (auto plt = fixedPltRoot.child("plt"); plt; plt = plt.next_sibling("plt"))
			builder->fixPlt(plt.text().as_uint());
	}

	std::cout << "post fixed" << std::endl;

	//last, set and add non resistance and non external constraints.
	auto model = builder->create();

	std::cout << "model built" << "\n";

	return model;
}


std::string generateOutputFileName(std::string inputFileName)
{
	time_t now;
	const int MAX_DATE = 64;
	char theDate[MAX_DATE];
	theDate[0] = '\0';

	now = time(nullptr);

	if (now != -1) {
		strftime(theDate, MAX_DATE, "_%Y.%m.%d_%H-%M-%S", gmtime(&now));
		return inputFileName + theDate;
	}
	return "";
}



void run(int argc, char** argv) {

	time_t t0,t1;
	t0 = time(0);

	double epsilon = 0.01;
	double timeStep = 0.001;

	for (int i = -1; i < argc-1; i++) {

		std::string arg = argv[i];
		unsigned pos = arg.find('=');

		std::string key = arg.substr(0, pos);
		std::string val = arg.substr(pos + 1);

		std::cout<<"argc: "<< argc <<std::endl;
		std::cout<<"arg: "<< arg <<std::endl;
		std::cout<<"pos: "<< pos <<std::endl;
		std::cout<<"key: "<< key <<std::endl;
		std::cout<<"val: "<< val <<std::endl;

		if (key == "-dt") {
			timeStep = std::atof(val.c_str());
			std::cout<<"setting timestep: "<< timeStep << std::endl;
			continue;
		}
		if (key == "-eps") {
			epsilon = std::atof(val.c_str());
			std::cout<<"setting epsilon: "<< epsilon << std::endl;
			continue;
		}
	}


		auto builder = std::make_shared<SystemBuilder>(epsilon, timeStep);

		auto system = createSystem(argv[argc-1], builder);

		auto outputFileName = generateOutputFileName(argv[argc-1]);

		//once the system is set, we'll store the initial values via the ptr system.
		//Storage storage( system, outputFileName);
		auto storage = std::make_shared<Storage>(system, builder, outputFileName);


		std::cout << "assigning fdiagram in main" << std::endl;
		system->assignStorage(storage);
		std::cout << "post fdiagram in main" << std::endl;

		std::cout << "solving system in main" << std::endl;
		system->solveSystem();


	t1 = time(0);  //current time at the end of solving the system.
	int total,hours,min,sec;
	total = difftime(t1,t0);
	hours = total / 3600;
	min = (total % 3600) / 60;
	sec = (total % 3600) % 60;
	std::cout << "Total time hh: " << hours << " mm:" << min << " ss:" << sec <<"\n";

}


int main(int argc, char** argv)
{
	std::cout << argc;
	std::cout << std::endl;

	run(argc - 2, argv + 2);

	return 0;
}

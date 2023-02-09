#include "simulator.h"
#include <iostream>
/*#include "consts.h"                    // defined constants
#include "macros.h"                    // macros used throughout SWMM
#include "enums.h"                     // enumerated variables
#include "error.h"                     // error message codes
#include "datetime.h"                  // date/time functions
#include "objects.h"                   // definitions of SWMM's data objects
#include "funcs.h"                     // declaration of all global functions
#include "text.h"                      // listing of all text strings 
#define  EXTERN                        // defined as 'extern' in headers.h
#include "globals.h"                   // declaration of all global variables

ar* argv[]);

#include "swmm5.h"                     // declaration of exportable functions
extern "C" { int Nnodes[MAX_NODE_TYPES];}
extern "C" {TNode* Node;}
*/
void initializeSWMM(void);
void closeSWMM(void);
void updateSources(void);
void close2D(void);
void initialize2D(std::string);
void initializeSources2D(void);


void setSourceLocation(void);
float flows[1] = { 0
};
 Simulator sim;
 float flow_input;
 //extern "C" int Nnodes[MAX_NODE_TYPES];
 int main(int argc, char *argv[]) {

	 std::cout << "start";

	 for (int i = 0; i < 1; i++){

		 //flow_input = atof(argv[2]);
		 flow_input = flows[i];
		 std::cout << "The flow is: " << std::endl;
		 std::cout << "going to initialize2D" <<std::endl;
		 //initialize2D(argv);
		
		std::string cfg;
		cfg =argv[1];// "./SampleData/sample.cfg"; 
		initialize2D(cfg.c_str());
		 std::cout << "This will run for " << sim.tf << " seconds" << std::endl;
		 while (sim.t < sim.tf && sim.steady_state){

			 //updateSources();
			 sim.RunSimulation();
			 //std::cout << "Updating sources...";

		 }

		
		 close2D();
 }
		/////////////////////////////////////////////////////////////////////////////////////////
    return 0;
}

void initialize2D(std::string cfg){
	std::cout << "start Initialize" << std::endl ;
	sim.OpenSimulation(cfg.c_str(), flow_input);
}

void close2D(){
	sim.CloseSimulation();
}

void initializeSources2D(){
	sim.InitializeSources(2);
}
void setSourceLocation(){
	int j;
	for (j = 0; j <= 1; j++){
		float x, y;
		switch (j){
		case 0:
			x = -80.951499;
			y = 34.037;
			break; 
		case 1:
			x = -80.951499;
			y = 34.037;
			break;

		}

		sim.SetSourceLocation(j, x, y);
	}
}

void updateSources(){
	int j;

	for (j = 0; j <= 1; j++){
		float flow = flow_input;
		switch (j){
		case 0:
			//This is Overcreek 
			flow /= 35.3147f;
			flow *= 1.0f;
			break;
		case 1:

			//This is nowhere...
			flow /= 35.3147f;
			flow *= .0f;
			break;

		}
		//std::cout << "now setting source" << std::endl;
		sim.setSourceValue(j, flow);
		//printf("%.2f ", Node[j].overflow);
	}
	sim.updateSources();
	
}


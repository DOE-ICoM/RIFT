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

		 //float t;

		 //sim.opensimulation(argv[1]);
		 //while (sim.t<sim.tf){
		 //	t = sim.runsimulation();
		 //}
		 //sim.closesimulation();


		 //Uncomment here for the 1D/2D
		 /*
		 /////////////////////////////////////////////////////////////////////////////////
		 initializeSWMM();
		 initialize2D();
		 long newHour, oldHour = 0;
		 long theDay, theHour;
		 DateTime elapsedTime = 0.0;
		 double swmmTime = 0.0;

		 while(sim.t<sim.tf){
		 if (sim.t <= swmmTime){
		 sim.RunSimulation();
		 }
		 if (swmmTime <= sim.t){
		 swmm_step(&elapsedTime);
		 swmmTime = elapsedTime*24*60*60;
		 updateSources();
		 }
		 }

		 */




		 //        do
		 //        {

		 //int crap;
		 //crap = Nnodes[JUNCTION];
		 //            swmm_step(&elapsedTime);

		 //            newHour = (long)(elapsedTime * 24.0);
		 //            if ( newHour > oldHour )
		 //            {
		 //                theDay = (long)elapsedTime;
		 //                theHour = (long)((elapsedTime - floor(elapsedTime)) * 24.0);
		 //               // writecon("\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		 //                sprintf(Msg, "%-5d hour: %-2d", theDay, theHour);
		 //               // writecon(Msg);
		 //                oldHour = newHour;
		 //            }


		 //

		 //        } while ( elapsedTime > 0.0 && !ErrorCode );

		 //			closeSWMM();
		 close2D();
 }
		/////////////////////////////////////////////////////////////////////////////////////////


		/*runTime = difftime(time(0), start);*/
        /*printf("\n\n... EPA-SWMM completed in %.2f seconds.", runTime);*/



    return 0;
}

//void initialize2D(char* argv[]){
void initialize2D(std::string cfg){
	std::cout << "start Initialize" << std::endl ;
	//sim.OpenSimulation(argv[1],flow_input);
	sim.OpenSimulation(cfg.c_str(), flow_input);
	/*Used for the LDRD Project
	initializeSources2D();
	setSourceLocation();
	updateSources();
	*/
}
/*
void initializeSWMM(){
	char *inputFile;
    char *reportFile;
    char *binaryFile;
    char blank[] = "";
    time_t start;
    double runTime;
	start = time(0);
	  inputFile = "d:/projects/nisac/flood_development/1d2d/swmm_nowhere/example1.inp";
        reportFile = "c:/temp/report.rpt";
        

		///SWMM_RUN///////////////////////////////////////////////////////////////////////
		
	
    int j;
    // --- open the files & read input data
    ErrorCode = 0;
    swmm_open(inputFile, reportFile , binaryFile);
	
    // --- run the simulation if input data OK
   
    swmm_start(TRUE);

	for (j=0; j<=Nnodes[JUNCTION];j++){
		printf("\nNode %s X %.2f Y %.2f",Node[j].ID,Node[j].Xcoord,Node[j].Ycoord);
	}
       
}
*/
void close2D(){
	sim.CloseSimulation();
}
/*
void closeSWMM(){
	   // --- clean up
        swmm_end();
    

    // --- report results
    if ( Fout.mode == SCRATCH_FILE ) swmm_report();                            //(5.0.016 - LR)

    // --- close the system
    swmm_close();
}
*/
/*
void updateSources(){
	int j;
	for (j=0; j<=Nnodes[JUNCTION];j++){

		sim.setSourceValue(j,Node[j].overflow);
		//printf("%.2f ", Node[j].overflow);
	}
}
void setSourceLocation(){
	int j;
	for (j=0; j<=Nnodes[JUNCTION];j++){
		sim.SetSourceLocation(j,Node[j].Xcoord,Node[j].Ycoord);
	}
}

void initializeSources2D(){
	sim.InitializeSources(Nnodes[JUNCTION]);
}
*/
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


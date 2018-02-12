
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"

void ScriptsGenerator() {
   
    gROOT->SetStyle("Plain");

    const unsigned int total_number_of_events = 500000;
    const unsigned int number_of_scripts = 10;

    const double aw_fw_ratio = 0.123344784871;

    double cpu_per_event_simul[4];
    double cpu_per_event_reco[4];
    int number_of_events[4];

    number_of_events[0] = total_number_of_events;
    number_of_events[1] = total_number_of_events;
    number_of_events[2] = int(((double)total_number_of_events*aw_fw_ratio)/(aw_fw_ratio+1.));
    number_of_events[3] = int(((double)total_number_of_events)/(aw_fw_ratio+1.));

    cpu_per_event_simul[0] = 0.1; 
    cpu_per_event_simul[1] = 0.115; 
    cpu_per_event_simul[2] = 0.13; 
    cpu_per_event_simul[3] = 0.19; 

    cpu_per_event_reco[0] = 0.1305;
    cpu_per_event_reco[1] = 0.1505; 
    cpu_per_event_reco[2] = 0.21; 
    cpu_per_event_reco[3] = 0.35;//not known yet 

    char name[128];
    char simul_names[4][64] = {"source_bulk","source_surface","tracker_aw","tracker_fw"};

    int seeds[number_of_scripts][4];

    unsigned int nn = 0;
    ifstream myfile ("seeds.txt");
    if(myfile.is_open()) {
	while ( myfile.good() && nn < number_of_scripts) {
	  myfile >> seeds[nn][0] >> seeds[nn][1] >> seeds[nn][2] >> seeds[nn][3];  
	  std::cout << seeds[nn][0] <<" "<< seeds[nn][1] <<" "<< seeds[nn][2] <<" "<< seeds[nn][3] << std::endl;
	  nn++;
	}
	myfile.close();
    }
    else std::cout << "ERROR: seed file not found!" << std::endl;

    for (int i_script=0;i_script<number_of_scripts;i_script++) {
	
	sprintf(name,"simul_config_%s_%d",simul_names[0],i_script);
	ofstream myfile (name);
	if (myfile.is_open()) {
	  myfile << "#@key_label  \"name\"\n";
	  myfile << "#@meta_label \"type\"\n\n";
	  myfile << "[name=\"flsimulate\" type=\"flsimulate::section\"]\n";
	  myfile << "numberOfEvents : integer = 500000\n\n";
	  myfile << "[name=\"flsimulate.variantService\" type=\"flsimulate::section\"]\n";
	  myfile << "settings : string[4] = \"primary_events:generator=Bi214_Po214\" \n";
	  myfile << "\"vertexes:generator=source_pads_bulk\" \n";
	  myfile << "\"geometry:layout=Basic\" \n";
	  myfile << "\"simulation:output_profile=all_details\" \n\n";
	  myfile << "[name=\"flsimulate.simulation\" type=\"flsimulate::section\"]\n";
	  myfile << "rngEventGeneratorSeed         : integer = "<<seeds[i_script][0]<<"\n";
	  myfile << "rngVertexGeneratorSeed        : integer = "<<seeds[i_script][1]<<"\n";
	  myfile << "rngGeant4GeneratorSeed        : integer = "<<seeds[i_script][2]<<"\n";
	  myfile << "rngHitProcessingGeneratorSeed : integer = "<<seeds[i_script][3]<<"\n\n";
	  myfile << " #output_profile : string[1] = \"all_details\"";
	  myfile.close();
	}
	else cout << "Unable to open file" << endl;

      }
}

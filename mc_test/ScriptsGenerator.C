
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

    const unsigned int total_number_of_events = 1000000;
    const unsigned int number_of_scripts = 10;
    const bool do_simul = true, do_reco = true;

    const double aw_fw_ratio = 0.123344784871;

    double cpu_per_event_simul[4];
    double cpu_per_event_reco[4];
    int number_of_events[4];

    number_of_events[0] = total_number_of_events;
    number_of_events[1] = total_number_of_events;
    number_of_events[2] = int(((double)total_number_of_events*aw_fw_ratio)/(aw_fw_ratio+1.));
    number_of_events[3] = int(((double)total_number_of_events)/(aw_fw_ratio+1.)) + 1;
    for(int ii=0;ii<4;ii++) number_of_events[ii] /= number_of_scripts;

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
    char vertex_generator_name[4][64] = {"source_pads_bulk","source_pads_surface","field_wire_surface","anode_wire_surface"};

    if(do_simul) {

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
      
      sprintf(name,"./script_launcher_BiPo_simul_%dev.csh",total_number_of_events);
      ofstream f_launcher (name);
      if (f_launcher.is_open()) {
	
	f_launcher << "#!/bin/csh\n\n";
	f_launcher << "newgroup --temp nemo\n\n";
	
	for(int i_simul=0;i_simul<4;i_simul++) {
	  for (int i_script=0;i_script<number_of_scripts;i_script++) {
	    
	    sprintf(name,"./configs/temp/simul_config_BiPo_%s_%d",simul_names[i_simul],i_script);
	    ofstream f_config (name);
	    if (f_config.is_open()) {
	      f_config << "#@key_label  \"name\"\n";
	      f_config << "#@meta_label \"type\"\n\n";
	      f_config << "[name=\"flsimulate\" type=\"flsimulate::section\"]\n";
	      f_config << "numberOfEvents : integer = "<<number_of_events[i_simul]<<"\n\n";
	      f_config << "[name=\"flsimulate.variantService\" type=\"flsimulate::section\"]\n";
	      f_config << "settings : string[4] = \"primary_events:generator=Bi214_Po214\" \n";
	      f_config << "\"vertexes:generator="<<vertex_generator_name[i_simul]<<"\" \n";
	      f_config << "\"geometry:layout=Basic\" \n";
	      f_config << "\"simulation:output_profile=all_details\" \n\n";
	      f_config << "[name=\"flsimulate.simulation\" type=\"flsimulate::section\"]\n";
	      f_config << "rngEventGeneratorSeed         : integer = "<<seeds[i_script][0]<<"\n";
	      f_config << "rngVertexGeneratorSeed        : integer = "<<seeds[i_script][3]<<"\n";
	      f_config << "rngGeant4GeneratorSeed        : integer = "<<seeds[i_script][1]<<"\n";
	      f_config << "rngHitProcessingGeneratorSeed : integer = "<<seeds[i_script][2]<<"\n\n";
	      f_config << "output_profile : string[1] = \"all_details\"\n";
	      f_config.close();
	    }
	    else cout << "Unable to open file" << endl;
	    
	    sprintf(name,"./scripts/temp/script_simul_BiPo_%s_%d",simul_names[i_simul],i_script);
	    ofstream f_script (name);
	    if (f_script.is_open()) {
	      f_script << "#!/bin/csh\n\n";
	      f_script << "cd $HOME\n\n";
	      f_script << "source ConfigNemoLyon.csh\n\n";
	      f_script << "cd analysis/mc_test/\n\n";
	      f_script << "flsimulate -c $SIMUL_CONFIG_PATH/simul_config_BiPo_"<< simul_names[i_simul] <<"_"<< i_script <<" -o $NEMO_PATH/analysis/mc_test/brio_files/temp/BiPo_"<< simul_names[i_simul] <<"_"<< total_number_of_events <<"ev_"<< i_script <<".brio\n";
	      f_script.close();
	    }
	    else cout << "Unable to open file" << endl;
	    
	    int cpu = (int)(cpu_per_event_simul[i_simul] * number_of_events[i_simul] * 2.);
	    if(cpu > 100000) { cout << "WARNING: cpu time exceeds 100000" << endl; cpu = 100000; }
	    f_launcher << "qsub -P P_nemo -l xrootd=1,sps=1,ct="<< cpu <<",fsize=30G,s_rss=10G -j y "<< name <<"\n"; 
	  }
	}
	f_launcher << "exit\n";
	f_launcher.close();
      }
      else cout << "Unable to open file" << endl;
    }

    if(do_reco) {
      
      sprintf(name,"./script_launcher_BiPo_reco_%dev",total_number_of_events);
      ofstream f_launcher (name);
      if (f_launcher.is_open()) {
	
	f_launcher << "#!/bin/csh\n\n";
	f_launcher << "newgroup --temp nemo\n\n";
	
	for(int i_simul=0;i_simul<4;i_simul++) {
	  for (int i_script=0;i_script<number_of_scripts;i_script++) {
	    
	    sprintf(name,"./configs/temp/reco_config_BiPo_%s_%d",simul_names[i_simul],i_script);
	    ofstream f_config (name);
	    if (f_config.is_open()) {

	      f_config << "#@key_label  \"name\"\n";
	      f_config << "#@meta_label \"type\"\n\n";
	      f_config << "[name=\"flreconstruct\" type=\"flreconstruct::section\"]\n";
	      f_config << "#@config Basic setup\n";
	      f_config << "# #@description Reconstruction version (automatic: extracted from input metadata)\n\n";
	      f_config << "# experimentalSetupUrn : string = \"urn:snemo:demonstrator:setup:1.0\"\n";
	      f_config << "#@description Number of events to reconstruct (default: 0 = no limit)\n";
	      f_config << "numberOfEvents : integer = "<<number_of_events[i_simul]<<"\n\n";
	      f_config << "[name=\"flreconstruct.plugins\" type=\"flreconstruct::section\"]\n";
	      f_config << "plugins : string[5] = \"Falaise_CAT\" \\ \n";
	      f_config << "\"TrackFit\" \\ \n";
	      f_config << "\"Falaise_TrackFit\" \\ \n";
	      f_config << "\"Falaise_ChargedParticleTracking\" \\ \n";
	      f_config << "\"BiPo\" \\ \n\n";
	      f_config << "[name=\"pipeline\" type=\"dpp::chain_module\"] \n";
	      f_config << "modules : string[6] = \"CalibrateTracker\" \\ \n";
	      f_config << "\"CalibrateCalorimeters\" \\ \n";
	      f_config << "\"CATTrackerClusterizer\" \\ \n";
	      f_config << "\"TrackFitting\" \\ \n";
	      f_config << "\"ChargedParticleTracker\" \\ \n";
	      f_config << "\"BiPo\"\n\n";
	      f_config << "[name=\"CalibrateTracker\" type=\"snemo::processing::mock_tracker_s2c_module\"]\n";
	      f_config << "logging.priority : string = \"warning\"\n\n";
	      f_config << "[name=\"CalibrateCalorimeters\" type=\"snemo::processing::mock_calorimeter_s2c_module\"]\n";
	      f_config << "logging.priority : string = \"warning\"\n\n";
	      f_config << "[name=\"CATTrackerClusterizer\" type=\"snemo::reconstruction::cat_tracker_clustering_module\"]\n";
	      f_config << "TPC.processing_prompt_hits    : boolean = true\n";
	      f_config << "TPC.processing_delayed_hits   : boolean = true\n\n";
	      f_config << "#@description CAT logging level\n";
	      f_config << "CAT.level : string = \"mute\"\n\n";
	      f_config << "[name=\"TrackFitting\" type=\"snemo::reconstruction::trackfit_tracker_fitting_module\"]\n";
	      f_config << "fitting_models		: string[2] = \"helix\" \"line\"\n";
	      f_config << "line.only_guess		: string[4] = \"BB\" \"TB\" \"BT\" \"TT\"\n";
	      f_config << "helix.only_guess	: string[8] = \"BBB\" \"BBT\" \"BTT\" \"TBB\" \"TBT\" \"TTB\" \"BTB\" \"TTT\"\n\n";
	      f_config << "#@description Fit the clusters of delayed Geiger hits\n";
	      f_config << "line.guess.fit_delayed_clusters : boolean = 1\n\n";
	      f_config << "#@description Track fit adds start time as an additionnal parameter to the fit (needs a calibration driver)\n";
	      f_config << "line.fit.fit_start_time    : boolean = 1\n\n";
	      f_config << "#@description Track fit recomputes the drift distance from drift time (needs a calibration driver)\n";
	      f_config << "line.fit.using_drift_time  : boolean = 0\n\n";
	      f_config << "[name=\"ChargedParticleTracker\" type=\"snemo::reconstruction::charged_particle_tracking_module\"]\n";
	      f_config << "drivers : string[4] = \"VED\" \"CCD\" \"CAD\" \"AFD\"\n";
	      f_config << "AFD.minimal_delayed_time                : real as time = 10 us\n";
	      f_config << "AFD.minimal_cluster_xy_search_distance  : real as length = 40 cm\n";
	      f_config << "AFD.minimal_cluster_z_search_distance   : real as length = 30 cm\n\n";
	      f_config << "[name=\"BiPo\" type=\"BiPo\"]\n";
	      f_config << "logging.priority : string = \"warning\"\n\n";
	      f_config << "#@description ROOT filename\n";
	      f_config << "filename : string as path = \"$NEMO_PATH/analysis/mc_test/root_files/temp/ BiPo_"<< simul_names[i_simul] <<"_"<< total_number_of_events <<"ev_"<< i_script <<".root\n";
	      f_config.close();
	    }
	    else cout << "Unable to open file" << endl;
	    
	    sprintf(name,"./scripts/temp/script_reco_BiPo_%s_%d",simul_names[i_simul],i_script);
	    ofstream f_script (name);
	    if (f_script.is_open()) {
	      f_script << "#!/bin/csh\n\n";
	      f_script << "cd $HOME\n\n";
	      f_script << "source ConfigNemoLyon.csh\n\n";
	      f_script << "cd analysis/mc_test/\n\n";
	      f_script << "flreconstruct -p $RECO_CONFIG_PATH/simul_config_BiPo_"<< simul_names[i_simul] <<"_"<< i_script <<" -i $NEMO_PATH/analysis/mc_test/brio_files/temp/BiPo_"<< simul_names[i_simul] <<"_"<< total_number_of_events <<"ev_"<< i_script <<".brio -o $NEMO_PATH/analysis/mc_test/brio_files/temp/BiPo_"<< simul_names[i_simul] <<"_"<< total_number_of_events <<"ev_"<< i_script <<".reco.brio\n";
	      f_script.close();
	    }
	    else cout << "Unable to open file" << endl;
	    
	    int cpu = (int)(cpu_per_event_reco[i_simul] * number_of_events[i_simul] * 2.);
	    if(cpu > 100000) { cout << "WARNING: cpu time exceeds 100000" << endl; cpu = 100000; }
	    f_launcher << "qsub -P P_nemo -l xrootd=1,sps=1,ct="<< cpu <<",fsize=30G,s_rss=10G -j y "<< name <<"\n"; 
	  }
	}
	f_launcher << "exit\n";
	f_launcher.close();
      }
      else cout << "Unable to open file" << endl;
    }




    

}

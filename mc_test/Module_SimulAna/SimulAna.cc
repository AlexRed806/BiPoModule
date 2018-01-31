//
//  SimulAna.cpp
//  
//
//  Created by Alessandro Minotti on 23/01/2018.
//

#include "SimulAna.h"

//#include "</Users/minotti/SuperNEMO/sw/Falaise/install/include/falaise/snemo/datamodels/particle_track_data.h"
#include <datamodels/particle_track_data.h>
#include <datamodels/tracker_trajectory_data.h>
#include <datamodels/tracker_trajectory_solution.h>
#include <datamodels/base_trajectory_pattern.h>
#include <datamodels/particle_track.h>

#include <bayeux/mctools/simulated_data.h>
#include <bayeux/mctools/base_step_hit.h>
#include <bayeux/geomtools/i_shape_1d.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

DPP_MODULE_REGISTRATION_IMPLEMENT(SimulAna,"SimulAna");

SimulAna::SimulAna() : dpp::base_module() {}

SimulAna::~SimulAna() {
    
    this->reset();
}

void SimulAna::initialize(const datatools::properties & setup_,
                          datatools::service_manager & service_mgr_,
                          dpp::module_handle_dict_type & module_dict_) {
    
    dpp::base_module::_common_initialize(setup_);

    //std::cout << "Module initialisation succesfull" << std::endl;

    // ROOT Tree        
    std::string root_filename = "SimulAna.root";
    if (setup_.has_key("filename")) {
        root_filename = setup_.fetch_string("filename");
        // To allow environment variable
        datatools::fetch_path_with_env(root_filename);
    }

    _root_file_ = new TFile(root_filename.c_str(), "RECREATE");
    _root_file_->cd();
    
    _root_tree_simulated_electrons_ = new TTree("tree_se", "tree_se");
    _root_tree_simulated_electrons_->SetDirectory(_root_file_);
    _root_tree_simulated_electrons_->Branch("sd.primary_electron_kinetic_energy", &_root_variables_simulated_electrons_.kinetic_energy);
    _root_tree_simulated_electrons_->Branch("sd.primary_electron_total_energy", &_root_variables_simulated_electrons_.total_energy);

    _root_tree_reconstructed_electrons_ = new TTree("tree_re", "tree_re");
    _root_tree_reconstructed_electrons_->SetDirectory(_root_file_);
    _root_tree_reconstructed_electrons_->Branch("sd.primary_electron_kinetic_energy", &_root_variables_simulated_electrons_.kinetic_energy);
    _root_tree_reconstructed_electrons_->Branch("ptd.reconstructed_electron_n_geiger_hits", &_root_variables_reconstructed_electrons_.n_geiger_hits);
    _root_tree_reconstructed_electrons_->Branch("ptd.reconstructed_electron_track_length", &_root_variables_reconstructed_electrons_.track_length);


    /*
    _root_histograms_.n_geiger_hits_electrons
      = new TH1I("n_geiger_hits_electrons","n_geiger_hits_electrons",100,0,50);
    _root_histograms_.track_length_electrons
      = new TH1F("track_length_electrons","track_length_electrons",100,0,3000);
    _root_histograms_.energy_deposit_electrons
      = new TH1F("energy_deposit_electrons","energy_deposit_electrons",100,0,4);
    */

    // Counters    
    _number_of_simulated_electrons_ = 0;
    _number_of_simulated_alphas_ = 0;
    _number_of_simulated_gammas_ = 0;

    _number_of_electrons_ = 0;
    _number_of_foil_electrons_ = 0;
    _number_of_wall_electrons_ = 0;
    _number_of_negative_charge_electrons_ = 0;

    _number_of_alphas_ = 0;
    _number_of_delayed_alphas_ = 0;
    _number_of_foil_alphas_ = 0;
    _number_of_nowall_alphas_ = 0;

    this->_set_initialized(true);    
}


dpp::base_module::process_status SimulAna::process(datatools::things & record_)
{
  //std::cout << "Module process called" << std::endl;

  // Event-wise counters //
  _number_of_event_electrons_ = 0;
  _number_of_event_alphas_ = 0;
  _number_of_event_gammas_ = 0;
    
    if (! record_.has("SD")) {
        //std::cout << "Ain't got SD" << std::endl;
        return PROCESS_CONTINUE;
    }
    else {
        
        const mctools::simulated_data & my_sd
        = record_.get<mctools::simulated_data>("SD");
        
        //std::cout << "Event time: " << my_sd.get_time() << std::endl;
        
        const mctools::simulated_data::primary_event_type my_primary_event
        = my_sd.get_primary_event();
        //std::cout << "Event time: " << my_primary_event.get_time() << std::endl;
        
        
        unsigned int N_primary = my_primary_event.get_number_of_particles();
        
        //std::cout << "Number of primary particles: " << N_primary << std::endl;
        
        genbb::primary_particle my_primary_particle[N_primary];
        for(unsigned int i_particle=0;i_particle<N_primary;i_particle++) {
            
            my_primary_particle[i_particle] = my_primary_event.get_particle(i_particle);
            //std::cout << "Particle type: " << my_primary_particle[i_particle].get_type() << std::endl;
            
            if(my_primary_particle[i_particle].is_electron()) {

	      _number_of_event_electrons_ ++;                
	      _number_of_simulated_electrons_ ++;

	      if (get_logging_priority() == datatools::logger::PRIO_TRACE) {
		DT_LOG_TRACE(get_logging_priority(), "Primary particle is electron");
                    
		//std::cout << "Electron Ek: " << my_primary_particle[i_particle].get_kinetic_energy() <<      " MeV" << std::endl;
	      }
              
	      _root_variables_simulated_electrons_.kinetic_energy = my_primary_particle[i_particle].get_kinetic_energy();
	      _root_variables_reconstructed_electrons_.kinetic_energy = my_primary_particle[i_particle].get_kinetic_energy();
	      _root_variables_simulated_electrons_.total_energy = my_primary_particle[i_particle].get_total_energy();
	      _root_tree_simulated_electrons_->Fill();
              
            }
	    else if(my_primary_particle[i_particle].is_alpha()) {

	      _number_of_event_alphas_ ++;
	      _number_of_simulated_alphas_ ++;

	    }
	    else if(my_primary_particle[i_particle].is_gamma()) {

	      _number_of_event_gammas_ ++;
	      _number_of_simulated_gammas_ ++;

	    }
        }
	
    }
    

    
    if (! record_.has("PTD")) {
        return PROCESS_CONTINUE;
    }
    const snemo::datamodel::particle_track_data & my_ptd = record_.get<snemo::datamodel::particle_track_data>("PTD");
        
    if (! my_ptd.has_particles()) {
        return PROCESS_CONTINUE;
    }

    //////////  Condition to select 1e1a events  //////////
    if (_number_of_event_electrons_ != 1 || _number_of_event_alphas_ != 1) {
      return PROCESS_CONTINUE;
    }
        
    //////////  LOOP OVER PARTICLES  //////////
    const snemo::datamodel::particle_track_data::particle_collection_type & my_particles = my_ptd.get_particles();
        
    for (const auto & ip : my_particles) {
        const snemo::datamodel::particle_track & my_pt = ip.get();
        
        is_helix = false;
	is_straight = false;
	is_prompt = false;
	is_delayed = false;
        does_touch_foil = false;
        does_touch_wall = false;
        has_negtive_charge = false;

        if(my_pt.get_charge()==snemo::datamodel::particle_track::charge_type::NEGATIVE)
            has_negtive_charge = true;
        
        //////////  STUDY TRAJECTORY  //////////
        snemo::datamodel::tracker_trajectory my_tj = my_pt.get_trajectory();
        
        if(my_tj.get_pattern().get_pattern_id() == "helix")
	  is_helix = true;
	else if(my_tj.get_pattern().get_pattern_id() == "line")
	  is_straight = true;
	
	_root_variables_reconstructed_electrons_.track_length = my_tj.get_pattern().get_shape().get_length();

	if(my_tj.has_cluster()) {
	  snemo::datamodel::tracker_cluster my_cl = my_tj.get_cluster();
	  _root_variables_reconstructed_electrons_.n_geiger_hits = my_cl.get_number_of_hits();
	  if(my_cl.is_delayed()) is_delayed = true;
	  else is_prompt = true;
	}
        /// end of study trajectory ///
            
        //////////  LOOP OVER VERTICES  //////////
        const snemo::datamodel::particle_track::vertex_collection_type & my_vertices = my_pt.get_vertices();
                
        for (const auto & iv : my_vertices) {
                
            if(snemo::datamodel::particle_track::vertex_is (iv.get(), snemo::datamodel::particle_track::vertex_type::VERTEX_ON_SOURCE_FOIL))
                does_touch_foil = true;
                    
            if(snemo::datamodel::particle_track::vertex_is (iv.get(), snemo::datamodel::particle_track::vertex_type::VERTEX_ON_MAIN_CALORIMETER))
                does_touch_wall = true;
        
        } /// end loop over vertices ///

        /*
        if (snemo::datamodel::particle_track::particle_has_negative_charge(my_pt)) {
                _root_variables_.electric_charge = "negative";
        }
            
        if (snemo::datamodel::particle_track::particle_has_negative_charge(my_pt)) {
            if (get_logging_priority() == datatools::logger::PRIO_TRACE) {
                DT_LOG_TRACE(get_logging_priority(), "Particle has negative charge");
                my_pt.tree_dump();
            }
        }*/
        
        if(is_prompt) { _number_of_electrons_++;
            if(does_touch_foil) {_number_of_foil_electrons_++;
                if(does_touch_wall) { _number_of_wall_electrons_++;
		  if(has_negtive_charge) {

		    _number_of_negative_charge_electrons_++;

		    //_root_variables_reconstructed_electrons_.n_geiger_hits = my_primary_particle[i_particle].get_total_energy();
		    //_root_variables_reconstructed_electrons_.track_length = my_primary_particle[i_particle].get_total_energy();
		    //_root_variables_reconstructed_electrons_.kinetic_energy = my_primary_particle[i_particle].get_kinetic_energy();
		    _root_tree_reconstructed_electrons_->Fill();

		    //_root_histograms_.n_geiger_hits_electrons->Fill(1)
		    //_root_histograms_.track_length_electrons->Fill(1)
		    //_root_histograms_.energy_depsit_electrons->Fill(1)
		  }
		}
            }
        }
        if(is_delayed) { _number_of_alphas_++;
	  if(is_delayed) {_number_of_delayed_alphas_++;
            if(does_touch_foil) {_number_of_foil_alphas_++;
	      if(!does_touch_wall) { _number_of_nowall_alphas_++;
	      }
            }
	  }
	}

        
    } /// end loop over particles ///
    


    
    return PROCESS_OK;
}

void SimulAna::reset() {
    
    if (_root_file_) {
        // write the output, finished streaming
        _root_file_->cd();
        _root_tree_simulated_electrons_->Write();
        _root_tree_reconstructed_electrons_->Write();
        _root_file_->Close();
        // clean up
        delete _root_file_;
	_root_tree_simulated_electrons_ = 0;
	_root_tree_reconstructed_electrons_ = 0;
        _root_file_ = 0;
    }
    
    SimulAna::print_results();
    
    this->_set_initialized(false);
}


void SimulAna::print_results() {

   std::cout << "Number of electrons: " << _number_of_electrons_ <<
     " (" << (double)_number_of_electrons_/(double)_number_of_simulated_electrons_<<
      " +/- " << pow((double)_number_of_electrons_,0.5)/(double)_number_of_simulated_electrons_<<
      " of total simulated electrons" << std::endl;
    std::cout << "Number of electrons touching the foil: " << _number_of_foil_electrons_ <<
      " (" << (double)_number_of_foil_electrons_/(double)_number_of_simulated_electrons_<<
      " +/- " << pow((double)_number_of_foil_electrons_,0.5)/(double)_number_of_simulated_electrons_<<
      " of total simulated electrons" << std::endl;
    std::cout << "Number of electrons touching the wall: " << _number_of_wall_electrons_ <<
      " (" << (double)_number_of_wall_electrons_/(double)_number_of_simulated_electrons_<<
      " +/- " << pow((double)_number_of_wall_electrons_,0.5)/(double)_number_of_simulated_electrons_<<
      " of total simulated electrons" << std::endl;
    std::cout << "Number of electrons with negatie charge: " << _number_of_negative_charge_electrons_ <<
      " (" << (double)_number_of_negative_charge_electrons_/(double)_number_of_simulated_electrons_<<
      " +/- " << pow((double)_number_of_negative_charge_electrons_,0.5)/(double)_number_of_simulated_electrons_<<
      " of total simulated electrons" << std::endl;

   std::cout << "Number of alphas: " << _number_of_alphas_ <<
      " (" << (double)_number_of_alphas_/(double)_number_of_simulated_alphas_<<
      " +/- " << pow((double)_number_of_alphas_,0.5)/(double)_number_of_simulated_alphas_<<
      " of total simulated alphas" << std::endl;
    std::cout << "Number of alphas that are delayed: " << _number_of_delayed_alphas_<<
      " (" << (double)_number_of_delayed_alphas_/(double)_number_of_simulated_alphas_<<
      " +/- " << pow((double)_number_of_delayed_alphas_,0.5)/(double)_number_of_simulated_alphas_<<
      " of total simulated alphas" << std::endl;
    std::cout << "Number of alphas touching the foil: " << _number_of_foil_alphas_ <<
      " (" << (double)_number_of_foil_alphas_/(double)_number_of_simulated_alphas_<<
      " +/- " << pow((double)_number_of_foil_alphas_,0.5)/(double)_number_of_simulated_alphas_<<
      " of total simulated alphas" << std::endl;
    std::cout << "Number of alphas not touching the wall: " << _number_of_nowall_alphas_ <<
      " (" << (double)_number_of_nowall_alphas_/(double)_number_of_simulated_alphas_<<
      " +/- " << pow((double)_number_of_nowall_alphas_,0.5)/(double)_number_of_simulated_alphas_<<
      " of total simulated alphas" << std::endl;

}


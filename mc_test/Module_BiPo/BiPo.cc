//
//  BiPo.cpp
//  
//  Created by Alessandro Minotti on 05/02/2018.
//

#include "BiPo.h"

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

DPP_MODULE_REGISTRATION_IMPLEMENT(BiPo,"BiPo");

BiPo::BiPo() : dpp::base_module() {}

BiPo::~BiPo() {
    
    this->reset();
}

void BiPo::initialize(const datatools::properties & setup_,
                          datatools::service_manager & service_mgr_,
                          dpp::module_handle_dict_type & module_dict_) {
    
    dpp::base_module::_common_initialize(setup_);

    // ROOT Tree        
    std::string root_filename = "BiPo.root";
    if (setup_.has_key("filename")) {
        root_filename = setup_.fetch_string("filename");
        // To allow environment variable
        datatools::fetch_path_with_env(root_filename);
    }

    _root_file_ = new TFile(root_filename.c_str(), "RECREATE");
    _root_file_->cd();
    
    _root_tree_simulated_electrons_ = new TTree("tree_mc_elec", "tree_mc_elec");
    _root_tree_simulated_electrons_->SetDirectory(_root_file_);
    _root_tree_simulated_electrons_->Branch("sd.primary_electrons_kinetic_energy", &_root_variables_simulated_particles_.kinetic_energy);
    _root_tree_simulated_electrons_->Branch("sd.primary_electrons_total_energy", &_root_variables_simulated_particles_.total_energy);
    
    _root_tree_simulated_alphas_ = new TTree("tree_mc_alpha", "tree_mc_alpha");
    _root_tree_simulated_alphas_->SetDirectory(_root_file_);
    _root_tree_simulated_alphas_->Branch("sd.primary_alphas_kinetic_energy", &_root_variables_simulated_particles_.kinetic_energy);
    _root_tree_simulated_alphas_->Branch("sd.primary_alphas_total_energy", &_root_variables_simulated_particles_.total_energy);
    _root_tree_simulated_alphas_->Branch("sd.primary_alphas_time", &_root_variables_simulated_particles_.emission_time);
    
    _root_tree_reconstructed_electrons_source_sel_ = new TTree("tree_rec_elec_source", "tree_rec_elec_source");
    _root_tree_reconstructed_electrons_source_sel_->SetDirectory(_root_file_);
    //_root_tree_reconstructed_electrons_source_sel_->Branch("sd.reconstructed_electrons_kinetic_energy", &_root_variables_simulated_particles_.kinetic_energy);
    _root_tree_reconstructed_electrons_source_sel_->Branch("ptd.reconstructed_electrons_n_geiger_hits", &_root_variables_reconstructed_particles_.n_geiger_hits);
    _root_tree_reconstructed_electrons_source_sel_->Branch("ptd.reconstructed_electrons_track_length", &_root_variables_reconstructed_particles_.track_length);
    
    _root_tree_reconstructed_electrons_tracker_sel_ = new TTree("tree_rec_elec_tracker", "tree_rec_elec_tracker");
    _root_tree_reconstructed_electrons_tracker_sel_->SetDirectory(_root_file_);
    //_root_tree_reconstructed_electrons_tracker_sel_->Branch("sd.reconstructed_electrons_kinetic_energy", &_root_variables_simulated_particles_.kinetic_energy);
    _root_tree_reconstructed_electrons_tracker_sel_->Branch("ptd.reconstructed_electrons_n_geiger_hits", &_root_variables_reconstructed_particles_.n_geiger_hits);
    _root_tree_reconstructed_electrons_tracker_sel_->Branch("ptd.reconstructed_electrons_track_length", &_root_variables_reconstructed_particles_.track_length);
    
    _root_tree_reconstructed_alphas_source_sel_ = new TTree("tree_rec_alpha_source", "tree_rec_alpha_source");
    _root_tree_reconstructed_alphas_source_sel_->SetDirectory(_root_file_);
    _root_tree_reconstructed_alphas_source_sel_->Branch("ptd.reconstructed_alphas_n_geiger_hits", &_root_variables_reconstructed_particles_.n_geiger_hits);
    _root_tree_reconstructed_alphas_source_sel_->Branch("ptd.reconstructed_alphas_track_length", &_root_variables_reconstructed_particles_.track_length);
    
    _root_tree_reconstructed_alphas_tracker_sel_ = new TTree("tree_rec_alpha_tracker", "tree_rec_alpha_tracker");
    _root_tree_reconstructed_alphas_tracker_sel_->SetDirectory(_root_file_);
    _root_tree_reconstructed_alphas_tracker_sel_->Branch("ptd.reconstructed_alphas_n_geiger_hits", &_root_variables_reconstructed_particles_.n_geiger_hits);
    _root_tree_reconstructed_alphas_tracker_sel_->Branch("ptd.reconstructed_alphas_track_length", &_root_variables_reconstructed_particles_.track_length);
    
    _root_tree_fitted_alphas_ = new TTree("tree_fit_alpha", "tree_fit_alpha");
    _root_tree_fitted_alphas_->SetDirectory(_root_file_);
    _root_tree_fitted_alphas_->Branch("ptd.reconstructed_alphas_n_geiger_hits", &_root_variables_reconstructed_particles_.n_geiger_hits);
    _root_tree_fitted_alphas_->Branch("ptd.reconstructed_alphas_track_length", &_root_variables_reconstructed_particles_.track_length);
    _root_tree_fitted_alphas_->Branch("ptd.reconstructed_alphas_delta_t",
       &_root_variables_reconstructed_particles_.delta_t_true_fit);
    
    _root_tree_reconstructed_1e1a_topology_source_sel_ = new TTree("tree_rec_1e1a_source", "tree_rec_1e1a_source");
    _root_tree_reconstructed_1e1a_topology_source_sel_->SetDirectory(_root_file_);
    _root_tree_reconstructed_1e1a_topology_source_sel_->Branch("ptd.reconstructed_alphas_track_length", &_root_variables_topologies_.alpha_track_length);
    _root_tree_reconstructed_1e1a_topology_source_sel_->Branch("ptd.reconstructed_1e1a_delta_t", &_root_variables_topologies_.delta_t_prompt_delayed);
    
    _root_tree_reconstructed_1e1a_topology_tracker_sel_ = new TTree("tree_rec_1e1a_tracker", "tree_rec_1e1a_tracker");
    _root_tree_reconstructed_1e1a_topology_tracker_sel_->SetDirectory(_root_file_);
    _root_tree_reconstructed_1e1a_topology_tracker_sel_->Branch("ptd.reconstructed_alphas_track_length", &_root_variables_topologies_.alpha_track_length);
    _root_tree_reconstructed_1e1a_topology_tracker_sel_->Branch("ptd.reconstructed_1e1a_delta_t", &_root_variables_topologies_.delta_t_prompt_delayed);
    
    
    // Counters    
    _number_of_simulated_electrons_ = 0;
    _number_of_simulated_alphas_ = 0;
    _number_of_simulated_gammas_ = 0;

    _number_of_simulated_1e1a_ = 0;

    _number_of_electrons_ = 0;
    _number_of_prompt_electrons_ = 0;
    _number_of_helix_electrons_ = 0;
    _number_of_foil_electrons_ = 0;
    _number_of_wall_electrons_ = 0;
    _number_of_anywall_electrons_ = 0;
    _number_of_negative_charge_electrons_ = 0;

    _number_of_alphas_ = 0;
    _number_of_delayed_alphas_ = 0;
    _number_of_foil_alphas_ = 0;
    _number_of_nowall_alphas_ = 0;
    _number_of_nofoil_alphas_ = 0;
    
    n_wierdos = 0;

    this->_set_initialized(true);    
}


dpp::base_module::process_status BiPo::process(datatools::things & record_) {
    
  // Event-wise counters //
  _number_of_event_electrons_ = 0;
  _number_of_event_alphas_ = 0;
  _number_of_event_gammas_ = 0;
    
    std::cout << "---------------------------------" << std::endl;
    std::cout << "---------   NEW EVENT   ---------" << std::endl;
    std::cout << "---------------------------------" << std::endl;


    ///////////////////////////////////////////
    ////////////  MONTE CARLO DATA  ///////////
    ///////////////////////////////////////////

    if (!record_.has("SD"))
        return PROCESS_CONTINUE;
    
    const mctools::simulated_data & my_sd = record_.get<mctools::simulated_data>("SD");
    
    const mctools::simulated_data::primary_event_type my_primary_event = my_sd.get_primary_event();
    
    /// Loop over primary particles ///
    unsigned int n_primary = my_primary_event.get_number_of_particles();
    genbb::primary_particle my_primary_particle[n_primary];

    for(unsigned int i_particle=0;i_particle<n_primary;i_particle++) {
        
        my_primary_particle[i_particle] = my_primary_event.get_particle(i_particle);
        
        _root_variables_simulated_particles_.kinetic_energy = my_primary_particle[i_particle].get_kinetic_energy();
        _root_variables_simulated_particles_.total_energy = my_primary_particle[i_particle].get_total_energy();
        _root_variables_simulated_particles_.emission_time = my_primary_particle[i_particle].get_time();

        if(my_primary_particle[i_particle].is_electron()) {

	      _number_of_event_electrons_ ++;                
	      _number_of_simulated_electrons_ ++;
            
          electron_true_ekin = _root_variables_simulated_particles_.kinetic_energy;
          electron_true_t = _root_variables_simulated_particles_.emission_time;
	      
          _root_tree_simulated_electrons_->Fill();
            
          //_root_variables_reconstructed_particles_.kinetic_energy = my_primary_particle[i_particle].get_kinetic_energy();
        }
        
	    else if(my_primary_particle[i_particle].is_alpha()) {

	      _number_of_event_alphas_ ++;
	      _number_of_simulated_alphas_ ++;
            
          alpha_true_t = _root_variables_simulated_particles_.emission_time;
            
          _root_tree_simulated_alphas_->Fill();
	    }
        
	    else if(my_primary_particle[i_particle].is_gamma()) {

	      _number_of_event_gammas_ ++;
	      _number_of_simulated_gammas_ ++;
        }
    } // end of loop over particles

    
    ///////////////////////////////////////////
    //////////  RECONSTRUCTION DATA  //////////
    ///////////////////////////////////////////
    
    if (!record_.has("PTD"))
        return PROCESS_CONTINUE;
    
    const snemo::datamodel::particle_track_data & my_ptd = record_.get<snemo::datamodel::particle_track_data>("PTD");
        
    if (!my_ptd.has_particles() || my_ptd.has_non_associated_calorimeters())
        return PROCESS_CONTINUE;
    
    //////////  Condition to select 1e1a events  //////////
    if (_number_of_event_electrons_ != 1 || _number_of_event_alphas_ != 1) {
        return PROCESS_CONTINUE;
    }
    _number_of_simulated_1e1a_ ++; //increment counter of searched topology
    
    
    //////////  LOOP OVER PARTICLES  //////////
    const snemo::datamodel::particle_track_data::particle_collection_type & my_particles = my_ptd.get_particles();
    
    _got_alpha_source_sel_ = false;
    _got_electron_source_sel_ = false;
    _got_alpha_tracker_sel_ = false;
    _got_electron_tracker_sel_ = false;
    
    for (const auto & ip : my_particles) {
        const snemo::datamodel::particle_track & my_pt = ip.get();
        
        is_helix = false;
        is_straight = false;
        is_prompt = false;
        is_delayed = false;
        does_touch_foil = false;
        does_touch_wall = false;
        does_touch_calo = false;
        has_negative_charge = false;
        has_no_charge = false;
        has_positive_charge = false;

        //////////  STUDY CHARGE  //////////
        if(my_pt.get_charge()==snemo::datamodel::particle_track::charge_type::NEGATIVE)
            has_negative_charge = true;
        
        else if(my_pt.get_charge()==snemo::datamodel::particle_track::charge_type::UNDEFINED)
            has_no_charge = true;
        
        else if(my_pt.get_charge()==snemo::datamodel::particle_track::charge_type::POSITIVE)
            has_positive_charge = true;

        else {n_wierdos++;}
            
        //////////  STUDY TRAJECTORY  //////////
        snemo::datamodel::tracker_trajectory my_tj = my_pt.get_trajectory();
        
        if(my_tj.get_pattern().get_pattern_id() == "helix") {
            is_helix = true;
            
            _root_variables_reconstructed_particles_.fitted_time = my_tj.get_auxiliaries().fetch_real_scalar("t0");
            electron_fitted_t = _root_variables_reconstructed_particles_.fitted_time;
        }
        else if(my_tj.get_pattern().get_pattern_id() == "line") {
            is_straight = true;
	
        //if(my_tj.has_auxiliaries())
            _root_variables_reconstructed_particles_.fitted_time = my_tj.get_auxiliaries().fetch_real_scalar("t0");
            alpha_fitted_t = _root_variables_reconstructed_particles_.fitted_time;
            _root_variables_topologies_.alpha_track_length = my_tj.get_pattern().get_shape().get_length();
        }
        _root_variables_reconstructed_particles_.track_length = my_tj.get_pattern().get_shape().get_length();

        
        if(my_tj.has_cluster()) {
          snemo::datamodel::tracker_cluster my_cl = my_tj.get_cluster();
          _root_variables_reconstructed_particles_.n_geiger_hits = my_cl.get_number_of_hits();
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
            
            if(snemo::datamodel::particle_track::vertex_is (iv.get(), snemo::datamodel::particle_track::vertex_type::VERTEX_ON_MAIN_CALORIMETER)
               ||snemo::datamodel::particle_track::vertex_is (iv.get(), snemo::datamodel::particle_track::vertex_type::VERTEX_ON_X_CALORIMETER)
               ||snemo::datamodel::particle_track::vertex_is (iv.get(), snemo::datamodel::particle_track::vertex_type::VERTEX_ON_GAMMA_VETO) )
                does_touch_calo = true;
        
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
        
        if(is_delayed && is_straight && has_no_charge) {
            
            if(alpha_fitted_t != 0 && alpha_true_t !=0)
                _root_variables_reconstructed_particles_.delta_t_true_fit = alpha_true_t - alpha_fitted_t;
            else _root_variables_reconstructed_particles_.delta_t_true_fit = pow(-1,0.5);
            _root_tree_fitted_alphas_->Fill();
        }
        
        //////////  CUT FLOW: SOURCE SEL  //////////
        if(is_helix) { _number_of_electrons_++;
            if(does_touch_foil) { _number_of_foil_electrons_++;
                if(does_touch_wall) { _number_of_wall_electrons_++;
                    if(has_negative_charge) {
                        _got_electron_source_sel_ = true;
                        _number_of_negative_charge_electrons_++;
                        //t_electron_source_sel = 0;
                        _root_tree_reconstructed_electrons_source_sel_->Fill();
                    }
                }
            }
        }
        if(is_delayed) { _number_of_alphas_++;
            if(is_delayed) {_number_of_delayed_alphas_++;
                if(does_touch_foil) {_number_of_foil_alphas_++;
                    if(!does_touch_wall) {
                        _got_alpha_source_sel_ = true;
                        _number_of_nowall_alphas_++;
                        //t_alpha_source_sel = 0;
                        _root_tree_reconstructed_alphas_source_sel_->Fill();
                    }
                }
            }
        }
        if(_got_electron_source_sel_ && _got_alpha_source_sel_) {
            
            _root_variables_topologies_.alpha_track_length = alpha_track_length;
            _root_variables_topologies_.delta_t_prompt_delayed = alpha_fitted_t - electron_fitted_t;
            _root_tree_reconstructed_1e1a_topology_source_sel_->Fill();
        }
        
        //////////  CUT FLOW: TRACKER SEL  //////////
        if(is_prompt) { _number_of_prompt_electrons_++;
            if(does_touch_wall) {
                _got_electron_tracker_sel_ = true;
                if(does_touch_calo) {
                    _number_of_anywall_electrons_++;
                    //t_electron_tracker_sel = 0;
                    _root_tree_reconstructed_electrons_tracker_sel_->Fill();
                }
            }
        }
        if(is_delayed) {
            if(is_delayed) {
                if(!does_touch_wall) {
                    if(!does_touch_foil) {
                        _got_alpha_tracker_sel_ = true;
                        _number_of_nofoil_alphas_++;
                        //t_alpha_tracker_sel = 0;
                        _root_tree_reconstructed_alphas_tracker_sel_->Fill();
                    }
                }
            }
        }
        if(_got_electron_tracker_sel_ && _got_alpha_tracker_sel_) {
            
            _root_variables_topologies_.alpha_track_length = alpha_track_length;
            _root_variables_topologies_.delta_t_prompt_delayed = alpha_fitted_t - electron_fitted_t;
            _root_tree_reconstructed_1e1a_topology_tracker_sel_->Fill();
        }


    } /// end loop over particles ///
    
    
    if(_got_electron_source_sel_ && _got_alpha_source_sel_) {
        
        
        
        
    }
    
    if(_got_electron_tracker_sel_ && _got_alpha_tracker_sel_) {
        
        
        
        
    }

    return PROCESS_OK;
}

void BiPo::reset() {
    
    if (_root_file_) {
        // write the output, finished streaming
        _root_file_->cd();
        _root_tree_simulated_electrons_->Write();
        _root_tree_simulated_alphas_->Write();
        _root_tree_reconstructed_electrons_source_sel_->Write();
        _root_tree_reconstructed_electrons_tracker_sel_->Write();
        _root_tree_reconstructed_alphas_source_sel_->Write();
        _root_tree_reconstructed_alphas_tracker_sel_->Write();
        _root_tree_fitted_alphas_->Write();
        _root_tree_reconstructed_1e1a_topology_source_sel_->Write();
        _root_tree_reconstructed_1e1a_topology_tracker_sel_->Write();
        _root_file_->Close();
        // clean up
        delete _root_file_;
        _root_tree_simulated_electrons_ = 0;
        _root_tree_reconstructed_electrons_source_sel_ = 0;
        _root_tree_reconstructed_electrons_tracker_sel_ = 0;
        _root_tree_reconstructed_alphas_source_sel_ = 0;
        _root_tree_reconstructed_alphas_tracker_sel_ = 0;
        _root_file_ = 0;
    }
    
    BiPo::print_results();
    
    this->_set_initialized(false);
}


void BiPo::print_results() {

    std::cout << "------  RESULTS OF THE CUT FLOW FOR SOURCE SELECTION  ------" << std::endl;
    
   std::cout << "Number of electrons: " << _number_of_electrons_ <<
     " - (" << 100*(double)_number_of_electrons_/(double)_number_of_simulated_1e1a_<<
      " +/- " << 100*pow((double)_number_of_electrons_,0.5)/(double)_number_of_simulated_1e1a_<<
      ")% of total simulated electrons" << std::endl;
    std::cout << "Number of electrons touching the foil: " << _number_of_foil_electrons_ <<
      " - (" << 100*(double)_number_of_foil_electrons_/(double)_number_of_simulated_1e1a_<<
      " +/- " << 100*pow((double)_number_of_foil_electrons_,0.5)/(double)_number_of_simulated_1e1a_<<
      ")% of total simulated electrons" << std::endl;
    std::cout << "Number of electrons touching the wall: " << _number_of_wall_electrons_ <<
      " - (" << 100*(double)_number_of_wall_electrons_/(double)_number_of_simulated_1e1a_<<
      " +/- " << 100*pow((double)_number_of_wall_electrons_,0.5)/(double)_number_of_simulated_1e1a_<<
      ")% of total simulated electrons" << std::endl;
    std::cout << "Number of electrons with negative charge: " << _number_of_negative_charge_electrons_ <<
      " - (" << 100*(double)_number_of_negative_charge_electrons_/(double)_number_of_simulated_1e1a_<<
      " +/- " << 100*pow((double)_number_of_negative_charge_electrons_,0.5)/(double)_number_of_simulated_1e1a_<<
      ")% of total simulated electrons" << std::endl;

   std::cout << "Number of alphas: " << _number_of_alphas_ <<
      " - (" << 100*(double)_number_of_alphas_/(double)_number_of_simulated_1e1a_<<
      " +/- " << 100*pow((double)_number_of_alphas_,0.5)/(double)_number_of_simulated_1e1a_<<
      ")% of total simulated alphas" << std::endl;
    std::cout << "Number of alphas that are delayed: " << _number_of_delayed_alphas_<<
      " - (" << 100*(double)_number_of_delayed_alphas_/(double)_number_of_simulated_1e1a_<<
      " +/- " << 100*pow((double)_number_of_delayed_alphas_,0.5)/(double)_number_of_simulated_1e1a_<<
      ")% of total simulated alphas" << std::endl;
    std::cout << "Number of alphas touching the foil: " << _number_of_foil_alphas_ <<
      " - (" << 100*(double)_number_of_foil_alphas_/(double)_number_of_simulated_1e1a_<<
      " +/- " << 100*pow((double)_number_of_foil_alphas_,0.5)/(double)_number_of_simulated_1e1a_<<
      ")% of total simulated alphas" << std::endl;
    std::cout << "Number of alphas not touching the wall: " << _number_of_nowall_alphas_ <<
      " - (" << 100*(double)_number_of_nowall_alphas_/(double)_number_of_simulated_1e1a_<<
      " +/- " << 100*pow((double)_number_of_nowall_alphas_,0.5)/(double)_number_of_simulated_1e1a_<<
      ")% of total simulated alphas" << std::endl;
    
    std::cout << "------  RESULTS OF THE CUT FLOW FOR TRACKER SELECTION  ------" << std::endl;

    std::cout << "Number of electrons that are prompt: " << _number_of_prompt_electrons_ <<
    " - (" << 100*(double)_number_of_prompt_electrons_/(double)_number_of_simulated_1e1a_<<
    " +/- " << 100*pow((double)_number_of_prompt_electrons_,0.5)/(double)_number_of_simulated_1e1a_<<
    ")% of total simulated electrons" << std::endl;
    std::cout << "Number of electrons touching ANY wall: " << _number_of_anywall_electrons_ <<
    " - (" << 100*(double)_number_of_anywall_electrons_/(double)_number_of_simulated_1e1a_<<
    " +/- " << 100*pow((double)_number_of_anywall_electrons_,0.5)/(double)_number_of_simulated_1e1a_<<
    ")% of total simulated electrons" << std::endl;
    
    std::cout << "Number of alphas that are delayed: " << _number_of_delayed_alphas_<<
    " - (" << 100*(double)_number_of_delayed_alphas_/(double)_number_of_simulated_1e1a_<<
    " +/- " << 100*pow((double)_number_of_delayed_alphas_,0.5)/(double)_number_of_simulated_1e1a_<<
    ")% of total simulated alphas" << std::endl;
    std::cout << "Number of alphas not touching the wall nor the foil: " << _number_of_nofoil_alphas_ <<
    " - (" << 100*(double)_number_of_nofoil_alphas_/(double)_number_of_simulated_1e1a_<<
    " +/- " << 100*pow((double)_number_of_nofoil_alphas_,0.5)/(double)_number_of_simulated_1e1a_<<
    ")% of total simulated alphas" << std::endl;
    
    std::cout << "Number of wierdos: " << n_wierdos << std::endl;

}


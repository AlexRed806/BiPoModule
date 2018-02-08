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
    
    // root trees with mc true information
    
    _root_tree_simulated_electrons_ = new TTree("tree_mc_elec", "tree_mc_elec");
    _root_tree_simulated_electrons_->SetDirectory(_root_file_);
    _root_tree_simulated_electrons_->Branch("sd.primary_electrons_kinetic_energy", &_root_variables_simulated_particles_.kinetic_energy);
    _root_tree_simulated_electrons_->Branch("sd.primary_electrons_total_energy", &_root_variables_simulated_particles_.total_energy);
    _root_tree_simulated_electrons_->Branch("sd.primary_electrons_time", &_root_variables_simulated_particles_.emission_time);

    _root_tree_simulated_alphas_ = new TTree("tree_mc_alpha", "tree_mc_alpha");
    _root_tree_simulated_alphas_->SetDirectory(_root_file_);
    _root_tree_simulated_alphas_->Branch("sd.primary_alphas_kinetic_energy", &_root_variables_simulated_particles_.kinetic_energy);
    _root_tree_simulated_alphas_->Branch("sd.primary_alphas_total_energy", &_root_variables_simulated_particles_.total_energy);
    _root_tree_simulated_alphas_->Branch("sd.primary_alphas_time", &_root_variables_simulated_particles_.emission_time);
    
    // root trees with ptd reconstructed particles information
    
    _root_tree_reconstructed_electrons_source_sel_ = new TTree("tree_rec_elec_source", "tree_rec_elec_source");
    _root_tree_reconstructed_electrons_source_sel_->SetDirectory(_root_file_);
    _root_tree_reconstructed_electrons_source_sel_->Branch("ptd.reconstructed_electrons_n_geiger_hits", &_root_variables_reconstructed_particles_.n_geiger_hits);
    _root_tree_reconstructed_electrons_source_sel_->Branch("ptd.reconstructed_electrons_track_length", &_root_variables_reconstructed_particles_.track_length);
    
    _root_tree_reconstructed_electrons_tracker_sel_ = new TTree("tree_rec_elec_tracker", "tree_rec_elec_tracker");
    _root_tree_reconstructed_electrons_tracker_sel_->SetDirectory(_root_file_);
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
    _root_tree_reconstructed_1e1a_topology_source_sel_->Branch("ptd.reconstructed_alphas_n_geiger_hits", &_root_variables_reconstructed_particles_.n_geiger_hits);
    _root_tree_reconstructed_1e1a_topology_source_sel_->Branch("ptd.reconstructed_1e1a_delta_t", &_root_variables_topologies_.delta_t_prompt_delayed);
    
    _root_tree_reconstructed_1e1a_topology_tracker_sel_ = new TTree("tree_rec_1e1a_tracker", "tree_rec_1e1a_tracker");
    _root_tree_reconstructed_1e1a_topology_tracker_sel_->SetDirectory(_root_file_);
    _root_tree_reconstructed_1e1a_topology_tracker_sel_->Branch("ptd.reconstructed_alphas_track_length", &_root_variables_topologies_.alpha_track_length);
    _root_tree_reconstructed_1e1a_topology_tracker_sel_->Branch("ptd.reconstructed_alphas_n_geiger_hits", &_root_variables_topologies_.alpha_n_geiger_hits);
    _root_tree_reconstructed_1e1a_topology_tracker_sel_->Branch("ptd.reconstructed_1e1a_delta_t", &_root_variables_topologies_.delta_t_prompt_delayed);
    
    
    // Counters
    
    _number_of_prompts_ = 0;
    _number_of_electrons_ = 0;
    _number_of_foil_electrons_ = 0;
    _number_of_wall_electrons_ = 0;
    _number_of_anywall_electrons_ = 0;
    _number_of_negative_charge_electrons_ = 0;

    _number_of_delayeds_ = 0;
    _number_of_alphas_ = 0;
    _number_of_foil_alphas_ = 0;
    _number_of_foil_nocalo_alphas_ = 0;
    _number_of_nofoil_nocalo_alphas_ = 0;
    _number_of_nofoil_alphas_ = 0;
    
    _number_of_1e1a_source_sel_ = 0;
    _number_of_1e1a_tracker_sel_ = 0;
    
    n_wierdos = 0;

    this->_set_initialized(true);    
}


dpp::base_module::process_status BiPo::process(datatools::things & record_) {
    
    //std::cout << "---------------------------------" << std::endl;
    //std::cout << "---------   NEW EVENT   ---------" << std::endl;
    //std::cout << "---------------------------------" << std::endl;

    /// creating my event to store all the info
    
    event _event_;
    
    
    ///////////////////////////////////////////
    ////////////  MONTE CARLO DATA  ///////////
    ///////////////////////////////////////////

    if (!record_.has("SD"))
        return PROCESS_CONTINUE;
    
    /// initialize the SD bank
    
    const mctools::simulated_data & my_sd = record_.get<mctools::simulated_data>("SD");
    
    const mctools::simulated_data::primary_event_type my_primary_event = my_sd.get_primary_event();
    
    
    /// Loop over primary particles ///
    unsigned int n_primary = my_primary_event.get_number_of_particles();
    genbb::primary_particle my_primary_particle[n_primary];
    
    for(unsigned int i_particle=0;i_particle<n_primary;i_particle++) {
        
        mc_particle _mc_particle_;
        
        my_primary_particle[i_particle] = my_primary_event.get_particle(i_particle);
        
        _mc_particle_.kinetic_energy = my_primary_particle[i_particle].get_kinetic_energy();
        _mc_particle_.emission_time = my_primary_particle[i_particle].get_time()/1000.;
        
        _root_variables_simulated_particles_.kinetic_energy = my_primary_particle[i_particle].get_kinetic_energy();
        _root_variables_simulated_particles_.total_energy = my_primary_particle[i_particle].get_total_energy();
        _root_variables_simulated_particles_.emission_time = my_primary_particle[i_particle].get_time()/1000.;

        if(my_primary_particle[i_particle].is_electron()) {
            
            _mc_particle_.type = "electron";
	      
            _root_tree_simulated_electrons_->Fill();
            
            _event_.number_of_electrons ++;
            _simulation_.number_of_electrons ++;
        }
        
	    else if(my_primary_particle[i_particle].is_alpha()) {
            
            _mc_particle_.type = "alpha";
            
            _root_tree_simulated_alphas_->Fill();
            
            _event_.number_of_alphas ++;
            _simulation_.number_of_alphas ++;
	    }
        
	    else if(my_primary_particle[i_particle].is_gamma()) {
            
            _mc_particle_.type = "gamma";

            _event_.number_of_gammas ++;
            _simulation_.number_of_gammas ++;
        }
        
        _event_.event_mc_particles.push_back(_mc_particle_);
        
    } // end of loop over particles

    
    ///////////////////////////////////////////
    //////////  RECONSTRUCTION DATA  //////////
    ///////////////////////////////////////////
    
    /// search and point the PTD cluster
    if (!record_.has("PTD"))
        return PROCESS_CONTINUE;
    
    const snemo::datamodel::particle_track_data & my_ptd = record_.get<snemo::datamodel::particle_track_data>("PTD");
    
    
    //////////  Condition to select 1e1a events  //////////
    if (_event_.number_of_electrons != 1 || _event_.number_of_alphas != 1)
        return PROCESS_CONTINUE;

    _simulation_.number_of_1e1a ++; //increment counter of searched topology
    
    
    /// check if there are reconstructed particles ///
    if (!my_ptd.has_particles())
        return PROCESS_CONTINUE;
    
    
    //////////  LOOP OVER PARTICLES  //////////
    const snemo::datamodel::particle_track_data::particle_collection_type & my_particles = my_ptd.get_particles();
    
    for (const auto & ip : my_particles) {
        const snemo::datamodel::particle_track & my_pt = ip.get();
        
        particle _particle_;

        //////////  STUDY CHARGE  //////////
        if(my_pt.get_charge()==snemo::datamodel::particle_track::charge_type::NEGATIVE)
            _particle_.charge = -1;
        
        else if(my_pt.get_charge()==snemo::datamodel::particle_track::charge_type::UNDEFINED)
            _particle_.charge = 0;

        else if(my_pt.get_charge()==snemo::datamodel::particle_track::charge_type::POSITIVE)
            _particle_.charge = 1;

        else {n_wierdos++;}
            
        //////////  STUDY TRAJECTORY  //////////
        if(!my_pt.has_trajectory()) std::cout << "Particle without trajectory?!?" << std::endl;
        
        snemo::datamodel::tracker_trajectory my_tj = my_pt.get_trajectory();
        
        
        /// first, check trajectory
        
        _particle_.track_length = my_tj.get_pattern().get_shape().get_length();
        
        /// if trajectory is helix, the reconstruced time is retrieved from first calorimeter hit (if present)
        
        if(my_tj.get_pattern().get_pattern_id() == "helix") {
        
            _particle_.trajectory_pattern = "curved";
            
            double temp;

            if(my_pt.has_associated_calorimeter_hits())
                temp = my_pt.get_associated_calorimeter_hits()[0].get().get_time() - 2.18484482E-9; //correction for average path
                
            else temp = pow(-1,0.5);
            
            _particle_.reconstructed_time = temp;
        }
        
        /// if trajectory is straight, the reconstruced time is retrieved from fit on trajectory

        else if(my_tj.get_pattern().get_pattern_id() == "line") {
            
            _particle_.trajectory_pattern = "straight";
            
            double temp = my_tj.get_auxiliaries().fetch_real_scalar("t0");
            
            if(temp > 2500.) temp /= 1000.; // THERE IS A BIG PROBLEM WITH UNITS HERE!!!
            
            // if alpha finder has done the fit, the fitted time is not good
            if(my_tj.get_cluster().get_number_of_hits() < 3) temp = pow(-1,0.5);
            
            _particle_.reconstructed_time = temp;
            _root_variables_reconstructed_particles_.fitted_time = temp;

        }
        else {n_wierdos++;}
        //std::cout << std::endl;

        _root_variables_reconstructed_particles_.track_length = my_tj.get_pattern().get_shape().get_length();
        
        
        /// then, loop over vertices ///
        const snemo::datamodel::particle_track::vertex_collection_type & my_vertices = my_pt.get_vertices();
        
        if(my_vertices.size() > 2) std::cout << "More than two vertices?!?" << std::endl;
        
        for (const auto & iv : my_vertices) {
            
            if(snemo::datamodel::particle_track::vertex_is (iv.get(), snemo::datamodel::particle_track::vertex_type::VERTEX_ON_SOURCE_FOIL))
                _particle_.does_hit_source_foil = true;
            
            else if(snemo::datamodel::particle_track::vertex_is (iv.get(), snemo::datamodel::particle_track::vertex_type::VERTEX_ON_MAIN_CALORIMETER))
                _particle_.does_hit_main_calo = true;
            
            else if(snemo::datamodel::particle_track::vertex_is (iv.get(), snemo::datamodel::particle_track::vertex_type::VERTEX_ON_X_CALORIMETER))
                _particle_.does_hit_x_calo = true;
            
            else if(snemo::datamodel::particle_track::vertex_is (iv.get(), snemo::datamodel::particle_track::vertex_type::VERTEX_ON_GAMMA_VETO))
                _particle_.does_hit_gamma_veto = true;
            
        } /// end loop over vertices ///

        
        /// and finally, check associated cluster ///

        if(my_tj.has_cluster()) {
            
          snemo::datamodel::tracker_cluster my_cl = my_tj.get_cluster();
            
            _particle_.number_of_geiger_hits = my_cl.get_number_of_hits();
            
            if(my_cl.is_delayed()) _particle_.is_delayed = true;
            else _particle_.is_delayed = false;

          _root_variables_reconstructed_particles_.n_geiger_hits = my_cl.get_number_of_hits();
        }
        else std::cout << "Trajectory without cluster?!?" << std::endl;
        /// end of study trajectory ///
            


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
        
        if(_particle_.is_delayed && _particle_.trajectory_pattern=="straight" && _particle_.charge==0) {
            
            for(std::vector<mc_particle>::iterator it = _event_.event_mc_particles.begin(); it != _event_.event_mc_particles.end(); ++it) {
                
                if( (*it).type == "alpha" && (*it).emission_time > 0.)
                    _root_variables_reconstructed_particles_.delta_t_true_fit = _root_variables_reconstructed_particles_.fitted_time - (*it).emission_time;
                
            }
            _root_tree_fitted_alphas_->Fill();
        }
        
        
        //////////  CUT FLOW  //////////
        
        if(!_particle_.is_delayed) { _number_of_prompts_++;
            if(_particle_.trajectory_pattern=="curved") { _number_of_electrons_++;
                if(_particle_.does_hit_source_foil) { _number_of_foil_electrons_++;
                    if(_particle_.does_hit_main_calo && !(_particle_.does_hit_x_calo || _particle_.does_hit_gamma_veto)) { _number_of_wall_electrons_++;
                        if(_particle_.charge==-1) {
                            _particle_.is_electron_source_sel = true;
                            _number_of_negative_charge_electrons_++;
                            _root_tree_reconstructed_electrons_source_sel_->Fill();
                        }
                    }
                }
            }
            if(_particle_.does_hit_main_calo || _particle_.does_hit_x_calo || _particle_.does_hit_gamma_veto) {
                _particle_.is_electron_tracker_sel = true;
                _number_of_anywall_electrons_++;
            }
        }
        
        if(_particle_.is_delayed) { _number_of_delayeds_++;
            if(_particle_.trajectory_pattern=="straight") { _number_of_alphas_++;
                if(_particle_.does_hit_source_foil) { _number_of_foil_alphas_++;
                    if(!(_particle_.does_hit_main_calo || _particle_.does_hit_x_calo || _particle_.does_hit_gamma_veto)) {
                        _particle_.is_alpha_source_sel = true;
                        _number_of_foil_nocalo_alphas_++;
                        _root_tree_reconstructed_alphas_source_sel_->Fill();
                    }
                }
                else { _number_of_nofoil_alphas_++;
                    if(!(_particle_.does_hit_main_calo || _particle_.does_hit_x_calo || _particle_.does_hit_gamma_veto)) {
                        _particle_.is_alpha_tracker_sel = true;
                        _number_of_nofoil_nocalo_alphas_++;
                        _root_tree_reconstructed_alphas_tracker_sel_->Fill();
                    }
                }
            }
        }
        ////////////////////////////////

        /// store particle in topology vector
        _event_.event_particles.push_back(_particle_);
        
    } /// end loop over particles ///


    /// now look for topologies
    
    if( _event_.event_particles.size() ==2 && (
                                               (_event_.event_particles[0].is_electron_source_sel && _event_.event_particles[1].is_alpha_source_sel)
                                               || (_event_.event_particles[1].is_electron_source_sel && _event_.event_particles[0].is_alpha_source_sel) ) ) {
        
        _number_of_1e1a_source_sel_++;
    
        double temp, temp1;
        
        for(std::vector<particle>::iterator it = _event_.event_particles.begin(); it != _event_.event_particles.end(); ++it) {
        
            if( (*it).is_alpha_source_sel ) {
        
                _root_variables_topologies_.alpha_track_length = (*it).track_length;
                _root_variables_topologies_.alpha_n_geiger_hits = (*it).number_of_geiger_hits;
                _event_.delayed_time = (*it).reconstructed_time;
            }
            else _event_.prompt_time = (*it).reconstructed_time;
        
        }
        
        _root_variables_topologies_.delta_t_prompt_delayed = _event_.delayed_time - _event_.prompt_time;

        _root_tree_reconstructed_1e1a_topology_source_sel_->Fill();
    }
    
    if( _event_.event_particles.size() ==2 && (
                                               (_event_.event_particles[0].is_electron_tracker_sel && _event_.event_particles[1].is_alpha_tracker_sel)
                                               || (_event_.event_particles[1].is_electron_tracker_sel && _event_.event_particles[0].is_alpha_tracker_sel) ) ) {
        _number_of_1e1a_tracker_sel_++;
        
        double temp, temp1;
        
        for(std::vector<particle>::iterator it = _event_.event_particles.begin(); it != _event_.event_particles.end(); ++it) {
            
            if( (*it).is_alpha_tracker_sel ) {
                
                _root_variables_topologies_.alpha_track_length = (*it).track_length;
                _root_variables_topologies_.alpha_n_geiger_hits = (*it).number_of_geiger_hits;
                _event_.delayed_time = (*it).reconstructed_time;
            }
            else _event_.prompt_time = (*it).reconstructed_time;
            
        }
        
        _root_variables_topologies_.delta_t_prompt_delayed = _event_.delayed_time - _event_.prompt_time;
        
        _root_tree_reconstructed_1e1a_topology_tracker_sel_->Fill();
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
    
   std::cout << "Number of prompts: " << _number_of_prompts_ <<
     " - (" << 100*(double)_number_of_prompts_/(double)_simulation_.number_of_1e1a<<
      " +/- " << 100*pow((double)_number_of_prompts_,0.5)/(double)_simulation_.number_of_1e1a<<
      ")% of total simulated events" << std::endl;
    std::cout << "Number of electrons: " << _number_of_electrons_ <<
    " - (" << 100*(double)_number_of_electrons_/(double)_simulation_.number_of_1e1a<<
    " +/- " << 100*pow((double)_number_of_electrons_,0.5)/(double)_simulation_.number_of_1e1a<<
    ")% of total simulated events" << std::endl;
    std::cout << "Number of electrons touching the foil: " << _number_of_foil_electrons_ <<
      " - (" << 100*(double)_number_of_foil_electrons_/(double)_simulation_.number_of_1e1a<<
      " +/- " << 100*pow((double)_number_of_foil_electrons_,0.5)/(double)_simulation_.number_of_1e1a<<
      ")% of total simulated events" << std::endl;
    std::cout << "Number of electrons touching the wall: " << _number_of_wall_electrons_ <<
      " - (" << 100*(double)_number_of_wall_electrons_/(double)_simulation_.number_of_1e1a<<
      " +/- " << 100*pow((double)_number_of_wall_electrons_,0.5)/(double)_simulation_.number_of_1e1a<<
      ")% of total simulated events" << std::endl;
    std::cout << "Number of electrons with negative charge: " << _number_of_negative_charge_electrons_ <<
      " - (" << 100*(double)_number_of_negative_charge_electrons_/(double)_simulation_.number_of_1e1a<<
      " +/- " << 100*pow((double)_number_of_negative_charge_electrons_,0.5)/(double)_simulation_.number_of_1e1a<<
      ")% of total simulated events" << std::endl;

    std::cout << "Number of delayeds: " << _number_of_delayeds_ <<
    " - (" << 100*(double)_number_of_delayeds_/(double)_simulation_.number_of_1e1a<<
    " +/- " << 100*pow((double)_number_of_delayeds_,0.5)/(double)_simulation_.number_of_1e1a<<
    ")% of total simulated events" << std::endl;
   std::cout << "Number of alphas: " << _number_of_alphas_ <<
      " - (" << 100*(double)_number_of_alphas_/(double)_simulation_.number_of_1e1a<<
      " +/- " << 100*pow((double)_number_of_alphas_,0.5)/(double)_simulation_.number_of_1e1a<<
      ")% of total simulated events" << std::endl;
    std::cout << "Number of alphas touching the foil: " << _number_of_foil_alphas_ <<
      " - (" << 100*(double)_number_of_foil_alphas_/(double)_simulation_.number_of_1e1a<<
      " +/- " << 100*pow((double)_number_of_foil_alphas_,0.5)/(double)_simulation_.number_of_1e1a<<
      ")% of total simulated events" << std::endl;
    std::cout << "Number of alphas not touching any calo: " << _number_of_foil_nocalo_alphas_ <<
      " - (" << 100*(double)_number_of_foil_nocalo_alphas_/(double)_simulation_.number_of_1e1a<<
      " +/- " << 100*pow((double)_number_of_foil_nocalo_alphas_,0.5)/(double)_simulation_.number_of_1e1a<<
      ")% of total simulated events" << std::endl;
    
    std::cout << "Number of 1e1a topologies: " << _number_of_1e1a_source_sel_ <<
    " - (" << 100*(double)_number_of_1e1a_source_sel_/(double)_simulation_.number_of_1e1a<<
    " +/- " << 100*pow((double)_number_of_1e1a_source_sel_,0.5)/(double)_simulation_.number_of_1e1a<<
    ")% of total simulated events" << std::endl;
    
    std::cout << "------  RESULTS OF THE CUT FLOW FOR TRACKER SELECTION  ------" << std::endl;

    std::cout << "Number of prompts: " << _number_of_prompts_ <<
    " - (" << 100*(double)_number_of_prompts_/(double)_simulation_.number_of_1e1a<<
    " +/- " << 100*pow((double)_number_of_prompts_,0.5)/(double)_simulation_.number_of_1e1a<<
    ")% of total simulated events" << std::endl;
    std::cout << "Number of electrons touching ANY calo: " << _number_of_anywall_electrons_ <<
    " - (" << 100*(double)_number_of_anywall_electrons_/(double)_simulation_.number_of_1e1a<<
    " +/- " << 100*pow((double)_number_of_anywall_electrons_,0.5)/(double)_simulation_.number_of_1e1a<<
    ")% of total simulated events" << std::endl;
    
    std::cout << "Number of delayeds: " << _number_of_delayeds_ <<
    " - (" << 100*(double)_number_of_delayeds_/(double)_simulation_.number_of_1e1a<<
    " +/- " << 100*pow((double)_number_of_delayeds_,0.5)/(double)_simulation_.number_of_1e1a<<
    ")% of total simulated events" << std::endl;
    std::cout << "Number of alphas: " << _number_of_alphas_ <<
    " - (" << 100*(double)_number_of_alphas_/(double)_simulation_.number_of_1e1a<<
    " +/- " << 100*pow((double)_number_of_alphas_,0.5)/(double)_simulation_.number_of_1e1a<<
    ")% of total simulated events" << std::endl;
    std::cout << "Number of alphas not touching the foil: " << _number_of_nofoil_alphas_ <<
    " - (" << 100*(double)_number_of_nofoil_alphas_/(double)_simulation_.number_of_1e1a<<
    " +/- " << 100*pow((double)_number_of_nofoil_alphas_,0.5)/(double)_simulation_.number_of_1e1a<<
    ")% of total simulated events" << std::endl;
    std::cout << "Number of alphas not touching any calo: " << _number_of_nofoil_nocalo_alphas_ <<
    " - (" << 100*(double)_number_of_nofoil_nocalo_alphas_/(double)_simulation_.number_of_1e1a<<
    " +/- " << 100*pow((double)_number_of_nofoil_nocalo_alphas_,0.5)/(double)_simulation_.number_of_1e1a<<
    ")% of total simulated events" << std::endl;
    
    std::cout << "Number of 1e1a topologies: " << _number_of_1e1a_tracker_sel_ <<
    " - (" << 100*(double)_number_of_1e1a_tracker_sel_/(double)_simulation_.number_of_1e1a<<
    " +/- " << 100*pow((double)_number_of_1e1a_tracker_sel_,0.5)/(double)_simulation_.number_of_1e1a<<
    ")% of total simulated alphas" << std::endl;
    
    
    std::cout << "Number of wierdos: " << n_wierdos << std::endl;

}


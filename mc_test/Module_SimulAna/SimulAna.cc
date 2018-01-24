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

#include <TFile.h>
#include <TTree.h>

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
    
    std::string root_filename = "SimulAna.root";
    if (setup_.has_key("filename")) {
        root_filename = setup_.fetch_string("filename");
        // To allow environment variable
        datatools::fetch_path_with_env(root_filename);
    }
    
    _root_file_ = new TFile(root_filename.c_str(), "RECREATE");
    _root_file_->cd();
    
    _root_tree_ = new TTree("tree", "tree");
    _root_tree_->SetDirectory(_root_file_);
    _root_tree_->Branch("sd.primary_electron_kinetic_energy", &_root_variables_.kinetic_energy);
    _root_tree_->Branch("sd.primary_electron_total_energy", &_root_variables_.total_energy);

    
    _number_of_electrons_ = 0;
    _number_of_foil_electrons_ = 0;
    _number_of_wall_electrons_ = 0;
    _number_of_negative_charge_electrons_ = 0;
    
    this->_set_initialized(true);
    
}


dpp::base_module::process_status SimulAna::process(datatools::things & record_)
{
    //std::cout << "Module process called" << std::endl;
    
    /*if (! record_.has("PTD")) {
        std::cout << "Ain't got PTD" << std::endl;
        return PROCESS_CONTINUE;
    }*/
    
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
                
                if (get_logging_priority() == datatools::logger::PRIO_TRACE) {
                    DT_LOG_TRACE(get_logging_priority(), "Primary particle is electron");
                    
                    //std::cout << "Electron Ek: " << my_primary_particle[i_particle].get_kinetic_energy() <<      " MeV" << std::endl;
                }
                
                _root_variables_.kinetic_energy = my_primary_particle[i_particle].get_kinetic_energy();
                
                _root_variables_.total_energy = my_primary_particle[i_particle].get_total_energy();
                
                _root_tree_->Fill();
                
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
        
    //////////  LOOP OVER PARTICLES  //////////
    const snemo::datamodel::particle_track_data::particle_collection_type & my_particles = my_ptd.get_particles();
        
    for (const auto & ip : my_particles) {
        const snemo::datamodel::particle_track & my_pt = ip.get();
        
        is_helix = false;
        does_touch_foil = false;
        does_touch_wall = false;
        has_negtive_charge = false;
        
        if(my_pt.get_charge()==snemo::datamodel::particle_track::charge_type::NEGATIVE)
            has_negtive_charge = true;
        
        //////////  STUDY TRAJECTORY  //////////
        snemo::datamodel::tracker_trajectory my_trajectory = my_pt.get_trajectory();
        
        if(my_trajectory.get_pattern().get_pattern_id() == "helix")
            is_helix = true;
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
        
        if(is_helix) { _number_of_electrons_++;
            if(does_touch_foil) {_number_of_foil_electrons_++;
                if(does_touch_wall) { _number_of_wall_electrons_++;
                    if(has_negtive_charge) _number_of_negative_charge_electrons_++;
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
        _root_tree_->Write();
        _root_file_->Close();
        // clean up
        delete _root_file_;
        _root_tree_ = 0;
        _root_file_ = 0;
    }
    
    std::cout << "Number of electrons: " << _number_of_electrons_ << std::endl;
    std::cout << "Number of electrons touching the foil: " << _number_of_foil_electrons_ << std::endl;
    std::cout << "Number of electrons touching the wall: " << _number_of_wall_electrons_ << std::endl;
    std::cout << "Number of electrons with negatie charge: " << _number_of_negative_charge_electrons_ << std::endl;
    
    this->_set_initialized(false);
}

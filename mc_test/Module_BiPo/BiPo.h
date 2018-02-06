//
//  BiPo.h
//  
//
//  Created by Alessandro Minotti on 23/01/2018.
//

#ifndef BiPo_h
#define BiPo_h

#include <stdio.h>
#include <string.h>

#include "bayeux/dpp/base_module.h"

// Forward declaration
class TFile;
class TTree;
class TH1I;
class TH1F;

// Structure to embed ROOT variable
struct root_variables_particles {
  //SD variables
  double kinetic_energy;
  double total_energy;
  double emission_time;
  //PTD variables
  int n_geiger_hits;
  double track_length;
  double fitted_time;
  double delta_t_true_fit;
};
struct root_variables_topologies {
  //int n_particles;
  //bool has_prompt;
  //bool has_delayed;
  double alpha_track_length;
  double delta_t_prompt_delayed;
};

class BiPo : public dpp::base_module {
    
public:
    
    BiPo();
    virtual ~BiPo();
    
    virtual void initialize(const datatools::properties& myConfig,
                            datatools::service_manager& flServices,
                            dpp::module_handle_dict_type& moduleDict);
    
    virtual dpp::base_module::process_status process(datatools::things& workItem);

    virtual void reset();
    void print_results();

private:

    TFile * _root_file_;
    TTree * _root_tree_simulated_electrons_;
    TTree * _root_tree_simulated_alphas_;
    TTree * _root_tree_reconstructed_electrons_source_sel_;
    TTree * _root_tree_reconstructed_electrons_tracker_sel_;
    TTree * _root_tree_reconstructed_alphas_source_sel_;
    TTree * _root_tree_reconstructed_alphas_tracker_sel_;
    TTree * _root_tree_fitted_alphas_;
    TTree * _root_tree_reconstructed_1e1a_topology_source_sel_;
    TTree * _root_tree_reconstructed_1e1a_topology_tracker_sel_;

    root_variables_particles _root_variables_simulated_particles_;
    root_variables_particles _root_variables_reconstructed_particles_;
    
    root_variables_topologies _root_variables_topologies_;
    
    bool is_helix;
    bool is_straight;
    bool is_prompt;
    bool is_delayed;
    bool does_touch_foil;
    bool does_touch_wall;
    bool does_touch_calo;
    bool has_negative_charge;
    bool has_no_charge;
    bool has_positive_charge;
    
    bool _got_alpha_source_sel_;
    bool _got_electron_source_sel_;
    bool _got_alpha_tracker_sel_;
    bool _got_electron_tracker_sel_;

    unsigned int _number_of_simulated_electrons_;
    unsigned int _number_of_simulated_alphas_;
    unsigned int _number_of_simulated_gammas_;
    unsigned int _number_of_simulated_1e1a_;

    unsigned int _number_of_event_electrons_;
    unsigned int _number_of_event_alphas_;
    unsigned int _number_of_event_gammas_;

    unsigned int _number_of_electrons_;
    unsigned int _number_of_prompt_electrons_;
    unsigned int _number_of_helix_electrons_;
    unsigned int _number_of_foil_electrons_;
    unsigned int _number_of_wall_electrons_;
    unsigned int _number_of_anywall_electrons_;
    unsigned int _number_of_negative_charge_electrons_;

    unsigned int _number_of_alphas_;
    unsigned int _number_of_delayed_alphas_;
    unsigned int _number_of_foil_alphas_;
    unsigned int _number_of_nowall_alphas_;
    unsigned int _number_of_nofoil_alphas_;

    unsigned int n_wierdos;
    
    //Temporary variables
    double alpha_true_t, alpha_fitted_t, electron_true_t, electron_fitted_t;
    double electron_true_ekin;

    DPP_MODULE_REGISTRATION_INTERFACE(BiPo)

    
};

#endif /* BiPo_h */

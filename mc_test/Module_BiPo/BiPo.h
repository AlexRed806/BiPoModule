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
struct root_variables {
  double kinetic_energy;
  double total_energy;
  int n_geiger_hits;
  double track_length;
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
    TTree * _root_tree_simulated_electrons_, * _root_tree_reconstructed_electrons_ , * _root_tree_reconstructed_alphas_;
    
    root_variables _root_variables_simulated_particles_, _root_variables_reconstructed_particles_;
    //root_histograms _root_histograms_;
    
    bool is_helix;
    bool is_straight;
    bool is_prompt;
    bool is_delayed;
    bool does_touch_foil;
    bool does_touch_wall;
    bool has_negative_charge;
    bool has_no_charge;
    bool has_positive_charge;

    unsigned int _number_of_simulated_electrons_;
    unsigned int _number_of_simulated_alphas_;
    unsigned int _number_of_simulated_gammas_;
    unsigned int _number_of_simulated_1e1a_;

    unsigned int _number_of_event_electrons_;
    unsigned int _number_of_event_alphas_;
    unsigned int _number_of_event_gammas_;

    unsigned int _number_of_electrons_;
    unsigned int _number_of_helix_electrons_;
    unsigned int _number_of_foil_electrons_;
    unsigned int _number_of_wall_electrons_;
    unsigned int _number_of_negative_charge_electrons_;

    unsigned int _number_of_alphas_;
    unsigned int _number_of_delayed_alphas_;
    unsigned int _number_of_foil_alphas_;
    unsigned int _number_of_nowall_alphas_;

    unsigned int n_wierdos;

    DPP_MODULE_REGISTRATION_INTERFACE(BiPo)

    
};

#endif /* BiPo_h */

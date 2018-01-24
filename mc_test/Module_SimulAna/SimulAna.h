//
//  SimulAna.h
//  
//
//  Created by Alessandro Minotti on 23/01/2018.
//

#ifndef SimulAna_h
#define SimulAna_h

#include <stdio.h>
#include <string.h>

#include "bayeux/dpp/base_module.h"

// Forward declaration
class TFile;
class TTree;

// Structure to embed ROOT variable
struct root_variables {
    double kinetic_energy;
    double total_energy;
};

class SimulAna : public dpp::base_module {
    
public:
    
    SimulAna();
    virtual ~SimulAna();
    
    virtual void initialize(const datatools::properties& myConfig,
                            datatools::service_manager& flServices,
                            dpp::module_handle_dict_type& moduleDict);
    
    virtual dpp::base_module::process_status process(datatools::things& workItem);

    virtual void reset();

private:

    TFile * _root_file_;
    TTree * _root_tree_;
    
    root_variables _root_variables_;
    
    bool is_helix;
    bool does_touch_foil;
    bool does_touch_wall;
    bool has_negtive_charge;
    
    unsigned int _number_of_electrons_;
    unsigned int _number_of_foil_electrons_;
    unsigned int _number_of_wall_electrons_;
    unsigned int _number_of_negative_charge_electrons_;
    
    DPP_MODULE_REGISTRATION_INTERFACE(SimulAna)

    
};

#endif /* SimulAna_h */

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

// forward declaration
class TFile;
class TTree;
class TH1I;
class TH1F;

// structures to embed ROOT variable
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
  double alpha_n_geiger_hits;
  double delta_t_prompt_delayed;
};

// structures to embed particles and events information
struct particle {
    std::string trajectory_pattern; //helix or line
    double reconstructed_time; //true emission time for SD, reconstructed for PTD (from calo for electrons, from fit for alphas)
    bool is_delayed; //1 if delayed, 0 if prompt
    int charge; //-1, 0, +1
    double track_length;
    unsigned int number_of_geiger_hits;
    bool does_hit_main_calo; //to be removed
    bool does_hit_x_calo; //to be removed
    bool does_hit_gamma_veto; //to be removed
    bool does_hit_source_foil; //to be removed
    std::vector<double> main_calo_vtx;
    std::vector<double> x_calo_vtx;
    std::vector<double> gamma_veto_vtx;
    std::vector<double> source_foil_vtx;
    bool is_electron_source_sel;
    bool is_alpha_source_sel;
    bool is_electron_tracker_sel;
    bool is_alpha_tracker_sel;
    
    particle(): trajectory_pattern("none"), reconstructed_time(pow(-1,0.5)), is_delayed(0), charge(pow(-1,0.5)), track_length(pow(-1,0.5)), number_of_geiger_hits(pow(-1,0.5)), does_hit_main_calo(0), does_hit_x_calo(0), does_hit_gamma_veto(0), does_hit_source_foil(0), is_electron_source_sel(0), is_alpha_source_sel(0), is_electron_tracker_sel(0), is_alpha_tracker_sel(0) {}
};
struct mc_particle {
    std::string type;
    double emission_time;
    double kinetic_energy;
};
struct event {
    unsigned int number_of_electrons;
    unsigned int number_of_alphas;
    unsigned int number_of_gammas;
    std::vector<mc_particle> event_mc_particles;
    std::vector<particle> event_particles;
    double prompt_time, delayed_time;
    
    event(): number_of_electrons(0), number_of_alphas(0), number_of_gammas(0), prompt_time(pow(-1,0.5)), delayed_time(pow(-1,0.5)) {}
};

struct simulation {
    unsigned int number_of_electrons = 0;
    unsigned int number_of_alphas = 0;
    unsigned int number_of_gammas = 0;
    unsigned int number_of_1e1a = 0;
    
    simulation(): number_of_electrons(0), number_of_alphas(0), number_of_gammas(0), number_of_1e1a(0) {}
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
    
    simulation _simulation_;
    
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

    unsigned int _number_of_prompts_;
    unsigned int _number_of_electrons_;
    unsigned int _number_of_foil_electrons_;
    unsigned int _number_of_wall_electrons_;
    unsigned int _number_of_anywall_electrons_;
    unsigned int _number_of_negative_charge_electrons_;

    unsigned int _number_of_delayeds_;
    unsigned int _number_of_alphas_;
    unsigned int _number_of_foil_alphas_;
    unsigned int _number_of_foil_nocalo_alphas_;
    unsigned int _number_of_nofoil_nocalo_alphas_;
    unsigned int _number_of_nofoil_alphas_;
    
    unsigned int _number_of_1e1a_source_sel_;
    unsigned int _number_of_1e1a_tracker_sel_;

    unsigned int n_wierdos;
    
    // temporary variables
    double alpha_true_t, alpha_fitted_t, electron_true_t, electron_fitted_t;
    double electron_true_ekin;

    DPP_MODULE_REGISTRATION_INTERFACE(BiPo)

    
};

#endif /* BiPo_h */

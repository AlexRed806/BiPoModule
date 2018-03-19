
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



//void TreeMerger(unsigned int n_evt, unsigned int n_script) {
void TreeMerger() {
   
    gROOT->SetStyle("Plain");

    const int n_simul = 4;
    const int n_trees = 9;
    const unsigned int total_number_of_events = 1000000;
    const unsigned int number_of_scripts = 20;

    const double aw_fw_ratio = 0.123344784871;

    char name[128];
    char simul_names[n_simul][64] = {"source_bulk","source_surface","tracker_aw","tracker_fw"};
    char tree_name[n_trees][64] = {"tree_mc_elec","tree_mc_alpha","tree_rec_elec_source","tree_rec_elec_tracker","tree_rec_alpha_source","tree_rec_alpha_tracker","tree_fit_alpha","tree_rec_1e1a_source","tree_rec_1e1a_tracker"};

    for(int i_simul=0;i_simul<n_simul;i_simul++) {
      
      TTree *input_tree[n_trees][number_of_scripts];
      TList *list[n_trees] = new TList;
      
      for (int i_script=0;i_script<number_of_scripts;i_script++) {
	
          sprintf(name,"./temp/BiPo_%s_%dev_%d.root",simul_names[i_simul],total_number_of_events,i_script);
	
          TFile *f = new TFile(name);
          
          for(int i_tree=0;i_tree<number_of_scripts;i_tree++) {
              
              input_tree[i_tree][i_script] = (TTree*)f->Get(tree_name[i_tree]);
              list->Add(input_tree[i_tree][i_script]);
          }
	
      }

      sprintf(name,"./BiPo_%s_%dev_merged.root",simul_names[i_simul],total_number_of_events);

      TFile *output_file = new TFile(name,"RECREATE");
      output_file->cd();
        TTree *output_tree[n_trees];
      for(int i_tree=0;i_tree<number_of_scripts;i_tree++) {

          output_tree[i_tree] = TTree::MergeTrees(list_source);
          output_tree[i_tree]->SetName(tree_name[i_tree]);
          output_tree[i_tree]->Write();
      }
    }
     


    TFile *f1 = new TFile("BiPo_tracker_aw_1000000ev_merged.root");
    TFile *f2 = new TFile("BiPo_tracker_fw_1000000ev_merged.root");
    
    TTree *t1[n_trees], *t2[n_trees], *t_out[n_trees];
    TList *l[n_trees];

    for(int i_tree=0;i_tree<number_of_scripts;i_tree++) {

        t1[i_tree] = (TTree*)f1->Get(tree_name[i_tree]);
        t2[i_tree] = (TTree*)f2->Get(tree_name[i_tree]);
        l[i_tree] = new TList;
        l[i_tree]->Add(t1[i_tree]);
        l[i_tree]->Add(t2[i_tree]);
    }

    TFile *fout = new TFile("BiPo_tracker_all_1000000ev_merged.root","RECREATE");
    fout->cd();
    
    for(int i_tree=0;i_tree<number_of_scripts;i_tree++) {

        t_out[i_tree] = TTree::MergeTrees(l[i_tree]);
        t_out[i_tree]->SetName(tree_name[i_tree]);
        t_out[i_tree]->Write();
    }
    fout->Close();


}

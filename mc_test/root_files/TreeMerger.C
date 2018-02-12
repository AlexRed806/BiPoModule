
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

void TreeMerger() {
   
    gROOT->SetStyle("Plain");

    const unsigned int total_number_of_events = 200000;
    const unsigned int number_of_scripts = 20;

    const double aw_fw_ratio = 0.123344784871;

    int number_of_events[4];

    number_of_events[0] = total_number_of_events;
    number_of_events[1] = total_number_of_events;
    number_of_events[2] = int(((double)total_number_of_events*aw_fw_ratio)/(aw_fw_ratio+1.));
    number_of_events[3] = int(((double)total_number_of_events)/(aw_fw_ratio+1.)) + 1;
    for(int ii=0;ii<4;ii++) number_of_events[ii] /= number_of_scripts;

    char name[128];
    char simul_names[4][64] = {"source_bulk","source_surface","tracker_aw","tracker_fw"};

    for(int i_simul=0;i_simul<4;i_simul++) {
      
      TTree *input_tree_source[number_of_scripts], *input_tree_tracker[number_of_scripts];
      TList *list_source = new TList;
      TList *list_tracker = new TList;
      
      for (int i_script=0;i_script<number_of_scripts;i_script++) {
	
	sprintf(name,"./temp/BiPo_%s_%dev_%d.root",simul_names[i_simul],total_number_of_events,i_script);
	
	TFile *f = new TFile(name);
	input_tree_source[i_script] = (TTree*)f->Get("tree_rec_1e1a_source");
	input_tree_tracker[i_script] = (TTree*)f->Get("tree_rec_1e1a_tracker");
	//TFile *root_file = new TFile(root_filename.c_str(), "READ");
	
	//f_config << "filename : string as path = \"$NEMO_PATH/analysis/mc_test/root_files/temp/BiPo_"<< simul_names[i_simul] <<"_"<< total_number_of_events <<"ev_"<< i_script <<".root";
	
        list_source->Add(input_tree_source[i_script]);
        list_tracker->Add(input_tree_tracker[i_script]);
	
      }

      sprintf(name,"./BiPo_%s_%dev_merged.root",simul_names[i_simul],total_number_of_events);

      TFile *output_file = new TFile(name,"RECREATE");
      output_file->cd();
      TTree *output_tree_source = TTree::MergeTrees(list_source);
      TTree *output_tree_tracker = TTree::MergeTrees(list_tracker);
      output_tree_source->SetName("tree_rec_1e1a_source");
      output_tree_tracker->SetName("tree_rec_1e1a_tracker");
      output_tree_source->Write();
      output_tree_tracker->Write();
    }
     


    TFile *f1 = new TFile("BiPo_tracker_aw_200000ev_merged.root");
    TFile *f2 = new TFile("BiPo_tracker_fw_200000ev_merged.root");
    TTree *t11 = (TTree*)f1->Get("tree_rec_1e1a_source");
    TTree *t12 = (TTree*)f1->Get("tree_rec_1e1a_tracker");
    TTree *t22 = (TTree*)f2->Get("tree_rec_1e1a_tracker");
    TTree *t21 = (TTree*)f2->Get("tree_rec_1e1a_source");
    TList *l1 = new TList;
    TList *l2 = new TList;
    l1->Add(t11);
    l1->Add(t21);
    l2->Add(t22);
    l2->Add(t12);
    TFile *f100 = new TFile("BiPo_tracker_all_200000ev_merged.root","RECREATE");
    root [16] f100->cd();
    TTree *output_tree_source = TTree::MergeTrees(l1);
    TTree *output_tree_tracker = TTree::MergeTrees(l2);
    output_tree_source->SetName("tree_rec_1e1a_source");
    output_tree_tracker->SetName("tree_rec_1e1a_tracker");
    output_tree_tracker->Write();
    output_tree_source->Write();
    f100->Close();


}


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

void Exposures() {

    using namespace RooFit;
    using namespace RooStats;
    
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(1111111);
    //gStyle->SetOptStat(0);
    
    gStyle->SetLegendTextSize(.03);
    gStyle->SetLabelFont(42,"XY");
    gStyle->SetLegendFont(132);
    gStyle->SetLegendBorderSize(1);
    gStyle->SetLegendFillColor(19);
    gStyle->SetStatFont(42);
    gStyle->SetTitleFont(132,"XY");
    gStyle->SetTitleOffset(1.25,"X");

    const int n_data = 8;
    const int n_type = 3;
    
    double n_days[n_data];
    double precision[n_type][n_data];
    
    char pdf_explicit_names[n_type][64] = {"fitted PDF of source bulk ^{214}Bi-Po","fitted PDF of source suface ^{214}Bi-Po","fitted PDF of tracker ^{214}Bi-Po"};
    
    int nn = 0;
    //ifstream myfile ("fitted_ab_source_sel.txt");
    ifstream myfile ("fitted_ab_tracker_sel.txt");
    if(myfile.is_open()) {
        while ( myfile.good() && nn < n_data) {
            myfile >> n_days[nn] >> precision[0][nn] >> precision[1][nn] >> precision[2][nn];
            //precision[0][nn]*=100;
            //precision[1][nn]*=100;
            //precision[2][nn]*=100;
            std::cout << n_days[nn] <<" "<< precision[0][nn] <<" "<< precision[1][nn] <<" "<< precision[2][nn] << std::endl;
            nn++;
        }
        myfile.close();
    }
    else std::cout << "ERROR: seed file not found!" << std::endl;
    
    
    TGraph* gr[n_type];
    
    TCanvas *c = new TCanvas("track_length","track_length",1961,344,700,502);
    TH1F *h = c->DrawFrame(0.8,0,150,150);
    //c->SetBorderMode(0);
    c->SetLogx();
    c->SetGridx();
    c->SetGridy();
    TLegend *l = new TLegend(0.14,0.65,0.55,0.89);

    for(int i_type=2;i_type<n_type;i_type++) {
        
        gr[i_type] = new TGraph(n_data,n_days,precision[i_type]);
        h->SetXTitle("exposure (days)");
        h->SetXTitle("exposure (days)");
        h->SetYTitle("precision (%)");
        gr[i_type]->SetLineWidth(4);
        gr[i_type]->SetLineColor(i_type+2);
        gr[i_type]->SetMarkerColor(i_type+2);
        gr[i_type]->SetMarkerStyle(7);
        //if(i_type==0) gr[i_type]->Draw("P");
        //else gr[i_type]->Draw("P");
        gr[i_type]->SetTitle(pdf_explicit_names[i_type]);
        gr[i_type]->Draw("CP");
        //gr[i_type]->SetLineWidth(1);
        l->AddEntry(gr[i_type],gr[i_type]->GetTitle(),"l");
        gr[i_type]->SetTitle("");
    }
    l->Draw();
    c->Modified();
    //c->Update();


    
}

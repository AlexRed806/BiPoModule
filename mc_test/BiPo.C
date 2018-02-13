
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

void BiPo(unsigned int n_ev, unsigned int n_ps, unsigned int n_day) {

    using namespace RooFit;
    using namespace RooStats;
    
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(1111111);
    //gStyle->SetOptStat(0);

    const int days_of_exposure = n_day;
    
    const unsigned int n_files = 3;
    const int n_events = n_ev;

    const double activities[n_files] = {0.0154,0.00018,0.0455};
    const double efficiencies[n_files] = {0.018342, 0.0866779, 0.0080547};
    //const double efficiencies[n_files] = {0.109888, 0.00555601, 0.0897256};
    const double given_exposure = 86400.*(double)days_of_exposure;
    const int n_pseudo = n_ps;
    
    const double fit_limit_coeffs[2] = {0.,5.};
    double fit_limits[n_files][2];
    
    double mc_exposures[n_files];
    for(int jj=0;jj<3;jj++) mc_exposures[jj] = (double)n_events / activities[jj];
    
    char name[56];
    char simul_names[3][54] = {"source_bulk","source_surface","tracker_all"};

    TH1D *h_track_length[n_files];
    for(int i_file=0; i_file<n_files; i_file++) {
        sprintf(name,"track_length_%s",simul_names[i_file]);
        h_track_length[i_file] = new TH1D(name,name,100,0,500);
    }
    TH1D *h_cumulative_track_length = new TH1D("cumulative_track_length","cumulative_track_length",100,0,500);

    TH1D *h_pseudo_track_length = new TH1D("pseudo_track_length","pseudo_track_length",100,0,500);

    for(int i_file=0; i_file<n_files; i_file++) {
        
        fit_limits[i_file][0] = activities[i_file]*fit_limit_coeffs[0];
        fit_limits[i_file][1] = activities[i_file]*fit_limit_coeffs[1];
    
        stringstream simul_name;
        simul_name << "./root_files/BiPo_"<<n_events<<"ev/BiPo_"<<simul_names[i_file]<<"_"<<n_events<<"ev_merged.root";
        cout << "Opening file " << simul_name.str() << endl;
 
        TFile *my_file = TFile::Open(simul_name.str().data());
        //TTree *tree_1e1a = (TTree*)my_file->Get("tree_rec_1e1a_source");
        //TTree *tree_1e1a = (TTree*)my_file->Get("tree_rec_1e1a_tracker");
        TChain *tree_1e1a = new TChain("tree_rec_1e1a_source");
        //TChain *tree_1e1a = new TChain("tree_rec_1e1a_tracker");
        tree_1e1a->Add(simul_name.str().data());
        
        
        int n_geiger;
        double track_length, delta_t;
        
        tree_1e1a->SetBranchAddress("ptd.reconstructed_alphas_n_geiger_hits",&n_geiger);
        tree_1e1a->SetBranchAddress("ptd.reconstructed_alphas_track_length",&track_length);
        tree_1e1a->SetBranchAddress("ptd.reconstructed_1e1a_delta_t",&delta_t);

        Long64_t nentries = tree_1e1a->GetEntries();
        for (Long64_t i=0;i<nentries;i++) {
            tree_1e1a->GetEntry(i);
            
            double exposure_correction = 1./(mc_exposures[i_file]);
            
            h_track_length[i_file]->Fill(track_length,exposure_correction);
            h_cumulative_track_length->Fill(track_length,exposure_correction);

        }

    }

    
    TCanvas *c = new TCanvas("track_length","track_length",1961,344,700,502);
    
    h_cumulative_track_length->SetTitle("");
    h_cumulative_track_length->GetYaxis()->SetTitle("");
    h_cumulative_track_length->GetXaxis()->SetTitle("#alpha track leength (mm)");

    h_cumulative_track_length->SetLineWidth(3);
    h_cumulative_track_length->SetLineColor(1);
        
    h_cumulative_track_length->Draw();
    
    
    TCanvas *c2 = new TCanvas("track_length_","track_length_",1961,344,700,502);
    for(int i_file=0; i_file<n_files; i_file++) {
            
        h_track_length[i_file]->SetLineWidth(3);
        h_track_length[i_file]->SetLineColor(i_file+2);

        if(i_file==0) h_track_length[i_file]->Draw();
        else h_track_length[i_file]->Draw("sames");
    }


    bool verbose(false);
    
    TH1D *h_fitted_activities[3];
    for(int i_file=0; i_file<n_files; i_file++) {

        sprintf(name,"fitted_actvity_%s",simul_names[i_file]);

        h_fitted_activities[i_file] = new TH1D(name,name,(int)(n_pseudo/10),fit_limits[i_file][0],fit_limits[i_file][1]);
    }
    
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    
    RooWorkspace w("w");
    
    RooRealVar track_length_obs("track_length_obs","track_length_obs",0,500);
    RooRealVar *fitted_activity[3];
    
    RooDataHist datahist("datahist","datahist",track_length_obs, h_cumulative_track_length);
    RooHistPdf pdf("pdf","pdf",track_length_obs,datahist,3);
    
    RooDataHist *datahists[3];
    RooHistPdf *pdfs[3];
    
    for(int i_file=0; i_file<n_files; i_file++) {
        
        stringstream dh_name, pdf_name, var_name;
        dh_name << "data_hist_" << simul_names[i_file];
        pdf_name << "pdf_" << simul_names[i_file];
        var_name << "fitted_activity_" << simul_names[i_file];
        
        datahists[i_file] = new RooDataHist(dh_name.str().data(),dh_name.str().data(),
                                            track_length_obs, h_track_length[i_file]);
        pdfs[i_file] = new RooHistPdf(pdf_name.str().data(),pdf_name.str().data(),track_length_obs,
                                      *datahists[i_file],3);
        //fitted_activity[i_file] = new RooRealVar(var_name.str().data(),var_name.str().data(),activities[i_file]*given_exposure*0.,activities[i_file]*given_exposure*10.);
        fitted_activity[i_file] = new RooRealVar(var_name.str().data(),var_name.str().data(),(given_exposure*efficiencies[i_file])*fit_limits[i_file][0],(given_exposure*efficiencies[i_file])*fit_limits[i_file][1]);
        //fitted_activity[i_file] = new RooRealVar(var_name.str().data(),var_name.str().data(),activities[i_file]*(given_exposure*efficiencies[i_file]),"s-1");
        fitted_activity[i_file]->setVal(activities[i_file]*(given_exposure*efficiencies[i_file]));

    }
    
    RooAddPdf model("model","model",RooArgSet(*pdfs[0],*pdfs[1],*pdfs[2]),RooArgList(*fitted_activity[0],*fitted_activity[1],*fitted_activity[2]));

    for(int i_pseudo=0;i_pseudo<n_pseudo;i_pseudo++) {
        
        RooDataSet * data = pdf.generate( track_length_obs, h_cumulative_track_length->Integral()*given_exposure);
        
        data->SetName("data");
        model.fitTo(*data,PrintLevel(-1));

        for(int i_file=0; i_file<n_files; i_file++) {
            
            h_fitted_activities[i_file]->Fill( fitted_activity[i_file]->getValV()/(given_exposure*efficiencies[i_file]) );
        }
        
        if(verbose) {
        
            TCanvas *c1 = new TCanvas("track_length_pseudo","track_length_pseudo",1961,344,700,502);
        
        
            RooPlot * pl = track_length_obs.frame(Title("pseudo track length distribution"));
            data->plotOn(pl);
            model.plotOn(pl);
            pl->SetTitle("");
            pl->Draw();

            cout << activities[0]*(given_exposure*efficiencies[0]) <<" "<< fitted_activity[0]->getValV() << endl;
            cout << activities[1]*(given_exposure*efficiencies[1]) <<" "<< fitted_activity[1]->getValV() << endl;
            cout << activities[2]*(given_exposure*efficiencies[2]) <<" "<< fitted_activity[2]->getValV() << endl;
        }
    }
    
    //IMPORTANT: fitted_activity[0]->getValV() is equivalent to h_track_length[0]->Integral()*given_exposure

        
    TCanvas *c_fa[3];
    for(int i_file=0;i_file<n_files;i_file++) {
        sprintf(name,"fitted_activity_%s",simul_names[i_file]);
        c_fa[i_file] = new TCanvas(name,name,1961,344,700,502);

        h_fitted_activities[i_file]->SetLineWidth(2);
        h_fitted_activities[i_file]->SetLineColor(i_file+2);
        h_fitted_activities[i_file]->Draw();

        TF1 *myfit = new TF1("myfit","gaus");

        h_fitted_activities[i_file]->Fit(myfit,"","",
                                         h_fitted_activities[i_file]->GetBinCenter(std::max((int)(n_pseudo/200),1)),
                                         h_fitted_activities[i_file]->GetBinCenter(h_fitted_activities[i_file]->GetNbinsX()+1));
        cout << "FITTED HISTOGRAM OF MOCK-DATA ACTIVITIES WITH GAUSSIAN FUNCTION" << endl;
        cout << "GAUSSIAN MEAN: " << myfit->GetParameter(1) << endl;
        cout << "GAUSSIAN SIGMA: " << myfit->GetParameter(2) << endl;
        cout << "ACCURACY OF ACTIVITY MEASUREMENT AFTER " << days_of_exposure << " DAYS: " << (myfit->GetParameter(2)/myfit->GetParameter(1))*100. << "%" <<endl;
    }

    
    
}

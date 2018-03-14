
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"

void BiPo(const unsigned int n_events, const unsigned int n_pseudo, const unsigned int days_of_exposure, const string type_of_sel) {

    using namespace RooFit;
    using namespace RooStats;
    
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1111);
    //gStyle->SetOptStat(1111111);
    gStyle->SetOptStat(0);
    
    gStyle->SetLegendTextSize(.03);
    gStyle->SetLabelFont(42,"XY");
    gStyle->SetLegendFont(132);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetStatFont(42);
    gStyle->SetTitleFont(132,"XY");

    const bool verbose = true;
    const unsigned int n_files = 3;
    const unsigned int first_file = 0;
    const unsigned int step = 1;
    const unsigned int poly_order = 3;
    const double n_bin_coeff = 10.;
    
    const double activities[n_files] = {0.0154,0.00018,0.0455};
    const double given_exposure = 86400.*(double)days_of_exposure;
    const double fit_limit_coeffs[2] = {-10.,50.};
    
    double efficiencies[n_files];
    double fit_limits[n_files][2];
    
    double mc_exposures[n_files];
    for(int jj=first_file;jj<n_files;jj+=step) mc_exposures[jj] = (double)n_events / activities[jj];
    
    char name[64], explicit_name[128];
    char simul_names[n_files][32] = {"source_bulk","source_surface","tracker_all"};
    char simul_explicit_names[n_files][64] = {"^{214}Bi-Po events in source bulk","^{214}Bi-Po events on source surface","^{214}Bi-Po events on tracker wires surface"};
    char pdf_explicit_names[n_files][64] = {"fitted PDF of source bulk ^{214}Bi-Po","fitted PDF of source suface ^{214}Bi-Po","fitted PDF of tracker ^{214}Bi-Po"};

    /// Root variables histograms definition ///
    
    TH1D *h_a_track_length[n_files];
    TH1D *h_a_n_geiger[n_files];
    TH1D *h_e_calo_energy[n_files];
    TH1D *h_delta_t[n_files];

    for(int i_file=first_file; i_file<n_files; i_file+=step) {
        
        sprintf(name,"a_track_length_%s",simul_names[i_file]);
        h_a_track_length[i_file] = new TH1D(name,simul_explicit_names[i_file],100,0,500);
        sprintf(name,"a_n_geiger_%s",simul_names[i_file]);
        h_a_n_geiger[i_file] = new TH1D(name,simul_explicit_names[i_file],10,0,10);
        sprintf(name,"e_calo_energy_%s",simul_names[i_file]);
        h_e_calo_energy[i_file] = new TH1D(name,simul_explicit_names[i_file],100,0,5);
        sprintf(name,"delta_t_%s",simul_names[i_file]);
        h_delta_t[i_file] = new TH1D(name,simul_explicit_names[i_file],100,0,2500000);

    }
    
    TH1D *h_cumulative_a_track_length = new TH1D("cumulative_a_track_length","^{214}Bi-Po events - combined vertexes",100,0,500);
    TH1D *h_cumulative_a_n_geiger = new TH1D("cumulative_a_n_geiger","^{214}Bi-Po events - combined vertexes",10,0,10);
    TH1D *h_cumulative_e_calo_energy = new TH1D("cumulative_e_calo_energy","^{214}Bi-Po events - combined vertexes",100,0,5);
    TH1D *h_cumulative_delta_t = new TH1D("cumulative_delta_t","^{214}Bi-Po events - combined vertexes",100,0,2500000);

    ///

    TH1D *h_pseudo_a_track_length = new TH1D("pseudo_a_track_length","pseudo_a_track_length",100,0,500);
    
    TH2D *h_vtx_top_view = new TH2D("h_vtx_top_view","h_vtx_top_view",100,-500,500,100,-3000,3000);
    TH2D *h_vtx_side_view = new TH2D("h_vtx_side_view","h_vtx_side_view",100,-500,500,100,-2000,2000);

    TH2D *hh_track_length_vs_source_dist = new TH2D("hh_track_length_vs_source_dist","hh_track_length_vs_source_dist",72,-10,350,51,0,510);

    for(int i_file=first_file; i_file<n_files; i_file+=step) {
        
        fit_limits[i_file][0] = activities[i_file]*fit_limit_coeffs[0];
        fit_limits[i_file][1] = activities[i_file]*fit_limit_coeffs[1];
    
        stringstream simul_name;
        simul_name << "./root_files/BiPo_"<<n_events<<"ev/BiPo_"<<simul_names[i_file]<<"_"<<n_events<<"ev_merged.root";
        cout << "Opening file " << simul_name.str() << endl;
 
        TFile *my_file = TFile::Open(simul_name.str().data());
        //TTree *tree_1e1a = (TTree*)my_file->Get("tree_rec_1e1a_source");
        //TTree *tree_1e1a = (TTree*)my_file->Get("tree_rec_1e1a_tracker");
        TChain *tree_1e1a;
        if(type_of_sel=="source") tree_1e1a = new TChain("tree_rec_1e1a_source");
        else if(type_of_sel=="tracker") tree_1e1a = new TChain("tree_rec_1e1a_tracker");
        else cerr << "ERROR: wrong selection type given as input" << endl;

        tree_1e1a->Add(simul_name.str().data());
        
        std::vector<double> *e_vtx = 0;
        double e_reco_energy;
        int a_n_geiger;
        double a_track_length, delta_t;
        
        tree_1e1a->SetBranchAddress("sd.primary_electrons_vertex",&e_vtx);
        tree_1e1a->SetBranchAddress("ptd.reconstructed_electrons_calo_energy",&e_reco_energy);
        tree_1e1a->SetBranchAddress("ptd.reconstructed_alphas_n_geiger_hits",&a_n_geiger);
        tree_1e1a->SetBranchAddress("ptd.reconstructed_alphas_track_length",&a_track_length);
        tree_1e1a->SetBranchAddress("ptd.reconstructed_1e1a_delta_t",&delta_t);

        Long64_t nentries = tree_1e1a->GetEntries();
        
        efficiencies[i_file] = (double)nentries/(0.973*(double)(n_events));
        cout << efficiencies[i_file] << endl;

        for (Long64_t i=0;i<nentries;i++) {
            tree_1e1a->GetEntry(i);
            
            
            if(i<10) {
                cout << delta_t << endl;
                //cout << e_vtx->size() << endl;
                //cout << e_vtx->at(0) <<" "<<e_vtx->at(1) <<" "<<e_vtx->at(2) <<endl;
                //cout << e_reco_energy << endl;
                //cout << "distance from source: " << abs((*e_vtx)[0]) << endl;;
            }
            
            double exposure_correction = 1./(mc_exposures[i_file]);
            
            h_a_track_length[i_file]->Fill(a_track_length,exposure_correction);
            h_a_n_geiger[i_file]->Fill(a_n_geiger,exposure_correction);
            h_e_calo_energy[i_file]->Fill(e_reco_energy,exposure_correction);
            h_delta_t[i_file]->Fill(delta_t,exposure_correction);
            
            h_cumulative_a_track_length->Fill(a_track_length,exposure_correction);
            h_cumulative_a_n_geiger->Fill(a_n_geiger,exposure_correction);
            h_cumulative_e_calo_energy->Fill(e_reco_energy,exposure_correction);
            h_cumulative_delta_t->Fill(delta_t,exposure_correction);

            h_vtx_top_view->Fill(((*e_vtx)[0]),((*e_vtx)[1]));
            h_vtx_side_view->Fill(((*e_vtx)[0]),((*e_vtx)[2]));

            hh_track_length_vs_source_dist->Fill(abs((*e_vtx)[0]),a_track_length);

        }

    }
    
    /// Plots ///
    
    TCanvas *c_vtx_top = new TCanvas("c_vtx_top","c_vtx_top",1961,344,700,502);
    h_vtx_top_view->Draw("colz");
    
    TCanvas *c_vtx_side = new TCanvas("c_vtx_side","c_vtx_side",1961,344,700,502);
    h_vtx_side_view->Draw("colz");
    
    TCanvas *c_length_vs_dist = new TCanvas("c_length_vs_dist","c_length_vs_dist",1961,344,700,502);
    c_length_vs_dist->SetLogz();
    //hh_track_length_vs_source_dist->SetMarkerStyle(7);
    hh_track_length_vs_source_dist->Draw("colz");

    
    TCanvas *c_a_track_length = new TCanvas("a_track_length","a_track_length",1961,344,700,502);
    TLegend *l_a_track_length = new TLegend(0.45,0.7,0.89,0.89);

    h_cumulative_a_track_length->GetXaxis()->SetTitle("#alpha track length (mm)");
    h_cumulative_a_track_length->GetYaxis()->SetTitle("(s^{-1})");
    h_cumulative_a_track_length->SetLineWidth(2);
    h_cumulative_a_track_length->SetLineColor(1);
    h_cumulative_a_track_length->SetFillStyle(0);
    h_cumulative_a_track_length->Draw("hist");
    
    l_a_track_length->AddEntry(h_cumulative_a_track_length,h_cumulative_a_track_length->GetTitle(),"lep");
    h_cumulative_a_track_length->SetTitle("");
    
    for(int i_file=first_file; i_file<n_files; i_file+=step) {
            
        h_a_track_length[i_file]->SetLineWidth(2);
        h_a_track_length[i_file]->SetLineColor(i_file+2);
        h_a_track_length[i_file]->SetFillColor(i_file+2);
        h_a_track_length[i_file]->SetFillStyle(3003);
        h_a_track_length[i_file]->Draw("histsames");
        
        l_a_track_length->AddEntry(h_a_track_length[i_file],h_a_track_length[i_file]->GetTitle(),"lep");
        h_a_track_length[i_file]->SetTitle("");
    }
    if(step == 1) h_a_track_length[1]->Draw("histsames");
    l_a_track_length->Draw();
    c_a_track_length->SaveAs("a_track_length.pdf");
    
    
    TCanvas *c_a_n_geiger = new TCanvas("a_n_geiger","a_n_geiger",1961,344,700,502);
    TLegend *l_a_n_geiger = new TLegend(0.45,0.7,0.89,0.89);
    
    h_cumulative_a_n_geiger->GetXaxis()->SetTitle("#alpha geiger hits");
    h_cumulative_a_n_geiger->GetYaxis()->SetTitle("(s^{-1})");
    h_cumulative_a_n_geiger->SetLineWidth(2);
    h_cumulative_a_n_geiger->SetLineColor(1);
    h_cumulative_a_n_geiger->SetFillStyle(0);
    h_cumulative_a_n_geiger->Draw("hist");
    
    l_a_n_geiger->AddEntry(h_cumulative_a_n_geiger,h_cumulative_a_n_geiger->GetTitle(),"lep");
    h_cumulative_a_n_geiger->SetTitle("");
    
    for(int i_file=first_file; i_file<n_files; i_file+=step) {
        
        h_a_n_geiger[i_file]->SetLineWidth(2);
        h_a_n_geiger[i_file]->SetLineColor(i_file+2);
        h_a_n_geiger[i_file]->SetFillColor(i_file+2);
        h_a_n_geiger[i_file]->SetFillStyle(3003);
        h_a_n_geiger[i_file]->Draw("histsames");
        
        l_a_n_geiger->AddEntry(h_a_n_geiger[i_file],h_a_n_geiger[i_file]->GetTitle(),"lep");
        h_a_n_geiger[i_file]->SetTitle("");
    }
    if(step == 1) h_a_n_geiger[1]->Draw("histsames");
    l_a_n_geiger->Draw();
    c_a_n_geiger->SaveAs("a_n_geiger.pdf");

    
    TCanvas *c_e_calo_energy = new TCanvas("e_calo_energy","e_calo_energy",1961,344,700,502);
    TLegend *l_e_calo_energy = new TLegend(0.45,0.7,0.89,0.89);
    
    h_cumulative_e_calo_energy->GetXaxis()->SetTitle("electron Edep in calo (MeV)");
    h_cumulative_e_calo_energy->GetYaxis()->SetTitle("(s^{-1})");
    h_cumulative_e_calo_energy->SetLineWidth(2);
    h_cumulative_e_calo_energy->SetLineColor(1);
    h_cumulative_e_calo_energy->SetFillStyle(0);
    h_cumulative_e_calo_energy->Draw("hist");
    
    l_e_calo_energy->AddEntry(h_cumulative_e_calo_energy,h_cumulative_e_calo_energy->GetTitle(),"lep");
    h_cumulative_e_calo_energy->SetTitle("");
    
    for(int i_file=first_file; i_file<n_files; i_file+=step) {
        
        h_e_calo_energy[i_file]->SetLineWidth(2);
        h_e_calo_energy[i_file]->SetLineColor(i_file+2);
        h_e_calo_energy[i_file]->SetFillColor(i_file+2);
        h_e_calo_energy[i_file]->SetFillStyle(3003);
        h_e_calo_energy[i_file]->Draw("histsames");
        
        l_e_calo_energy->AddEntry(h_e_calo_energy[i_file],h_e_calo_energy[i_file]->GetTitle(),"lep");
        h_e_calo_energy[i_file]->SetTitle("");
    }
    if(step == 1) h_e_calo_energy[1]->Draw("histsames");
    l_e_calo_energy->Draw();
    c_e_calo_energy->SaveAs("e_calo_energy.pdf");

    
    TCanvas *c_delta_t = new TCanvas("delta_t","delta_t",1961,344,700,502);
    TLegend *l_delta_t = new TLegend(0.45,0.7,0.89,0.89);
    
    c_delta_t->SetLogy();
    
    h_cumulative_delta_t->GetXaxis()->SetTitle("e-#alpha #Delta t (#mu s)");
    h_cumulative_delta_t->GetYaxis()->SetTitle("(s^{-1})");
    h_cumulative_delta_t->SetLineWidth(2);
    h_cumulative_delta_t->SetLineColor(1);
    h_cumulative_delta_t->SetFillStyle(0);
    h_cumulative_delta_t->Draw("hist");
    
    l_delta_t->AddEntry(h_cumulative_delta_t,h_cumulative_delta_t->GetTitle(),"lep");
    h_cumulative_delta_t->SetTitle("");
    
    for(int i_file=first_file; i_file<n_files; i_file+=step) {
        
        h_delta_t[i_file]->SetLineWidth(2);
        h_delta_t[i_file]->SetLineColor(i_file+2);
        h_delta_t[i_file]->SetFillColor(i_file+2);
        h_delta_t[i_file]->SetFillStyle(3003);
        h_delta_t[i_file]->Draw("histsames");
        
        l_delta_t->AddEntry(h_delta_t[i_file],h_delta_t[i_file]->GetTitle(),"lep");
        h_delta_t[i_file]->SetTitle("");
    }
    if(step == 1) h_delta_t[1]->Draw("histsames");
    l_delta_t->Draw();
    c_delta_t->SaveAs("delta_t.pdf");
    
    
    /// Statistical analysis ///
    
    TH1D *h_fitted_activities[n_files];
    for(int i_file=first_file; i_file<n_files; i_file+=step) {

        sprintf(name,"fitted_actvity_%s",simul_names[i_file]);

        h_fitted_activities[i_file] = new TH1D(name,name,(int)(n_pseudo*n_bin_coeff),fit_limits[i_file][0],fit_limits[i_file][1]);
    }
    
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(1);

    RooWorkspace w("w");
    
    RooRealVar track_length_obs("track_length_obs","track_length_obs",0,500);
    RooRealVar *fitted_activity[n_files];
    
    RooDataHist datahist("datahist","datahist",track_length_obs, h_cumulative_a_track_length);
    RooHistPdf pdf("pdf","pdf",track_length_obs,datahist,poly_order);
    
    RooDataHist *datahists[n_files];
    RooHistPdf *pdfs[n_files];
    
    for(int i_file=first_file; i_file<n_files; i_file+=step) {
        
        stringstream dh_name, pdf_name, var_name;
        dh_name << "data_hist_" << simul_names[i_file];
        pdf_name << "pdf_" << simul_names[i_file];
        var_name << "fitted_activity_" << simul_names[i_file];
        
        datahists[i_file] = new RooDataHist(dh_name.str().data(),dh_name.str().data(),
                                            track_length_obs, h_a_track_length[i_file]);
        pdfs[i_file] = new RooHistPdf(pdf_name.str().data(),pdf_name.str().data(),track_length_obs,
                                      *datahists[i_file],poly_order);
        //fitted_activity[i_file] = new RooRealVar(var_name.str().data(),var_name.str().data(),activities[i_file]*given_exposure*0.,activities[i_file]*given_exposure*10.);
        fitted_activity[i_file] = new RooRealVar(var_name.str().data(),var_name.str().data(),(given_exposure*efficiencies[i_file])*fit_limits[i_file][0],(given_exposure*efficiencies[i_file])*fit_limits[i_file][1]);
        //fitted_activity[i_file] = new RooRealVar(var_name.str().data(),var_name.str().data(),activities[i_file]*(given_exposure*efficiencies[i_file]),"s-1");
        fitted_activity[i_file]->setVal(activities[i_file]*(given_exposure*efficiencies[i_file]));

    }

    RooAddPdf model("model","model",RooArgSet(*pdfs[0],*pdfs[1],*pdfs[2]),RooArgList(*fitted_activity[0],*fitted_activity[1],*fitted_activity[2]));
    //RooAddPdf model("model","model",RooArgSet(*pdfs[0],*pdfs[2]),RooArgList(*fitted_activity[0],*fitted_activity[2]));

    for(int i_pseudo=0;i_pseudo<n_pseudo;i_pseudo++) {
        
        RooDataSet * data = pdf.generate( track_length_obs, h_cumulative_a_track_length->Integral()*given_exposure);
        
        data->SetName("data");
        model.fitTo(*data,PrintLevel(-1));

        for(int i_file=first_file; i_file<n_files; i_file+=step) {
            
            h_fitted_activities[i_file]->Fill( fitted_activity[i_file]->getValV()/(given_exposure*efficiencies[i_file]) );
        }

        if(verbose) {
        
            TCanvas *c_mock = new TCanvas("cacca","puzza",1961,344,700,502);
            TLegend *l_mock = new TLegend(0.5,0.65,0.89,0.89);
            
            //h_cumulative_a_track_length->SetTitle("");
            TH1 *h_cumulative_a_track_length_data = data->createHistogram("track_length_obs",100);
            h_cumulative_a_track_length_data->GetXaxis()->SetTitle("#alpha track length (mm)");
            //h_cumulative_a_track_length_data->GetYaxis()->SetTitle("(s^{-1})");
            h_cumulative_a_track_length_data->SetLineWidth(2);
            h_cumulative_a_track_length_data->SetLineColor(1);
            h_cumulative_a_track_length_data->SetFillStyle(0);
            h_cumulative_a_track_length_data->Draw();
            h_cumulative_a_track_length_data->SetTitle("mock data (30 days exposure)");
            l_mock->AddEntry(h_cumulative_a_track_length_data,h_cumulative_a_track_length_data->GetTitle(),"lep");
            h_cumulative_a_track_length_data->SetTitle("");

            TH1 *h_cumulative_a_track_length_mock = model.createHistogram("track_length_obs",100);
            //h_cumulative_a_track_length_mock->GetXaxis()->SetTitle("#alpha track length (mm)");
            //h_cumulative_a_track_length_mock->GetYaxis()->SetTitle("(s^{-1})");
            //h_cumulative_a_track_length_mock->Scale(fitted_activity[0]->getValV()+fitted_activity[1]->getValV()+fitted_activity[2]->getValV());
            h_cumulative_a_track_length_mock->Scale(5);
            h_cumulative_a_track_length_mock->SetLineWidth(2);
            h_cumulative_a_track_length_mock->SetLineColor(16);
            h_cumulative_a_track_length_mock->SetFillColor(16);
            h_cumulative_a_track_length_mock->SetFillStyle(3003);
            h_cumulative_a_track_length_mock->Draw("histsames");
            h_cumulative_a_track_length_mock->SetTitle("global PDF fitted to mock data");
            l_mock->AddEntry(h_cumulative_a_track_length_mock,h_cumulative_a_track_length_mock->GetTitle(),"l");
            h_cumulative_a_track_length_mock->SetTitle("");

            TH1 *h_a_track_length_mock[n_files];
            
            for(int i_file=first_file; i_file<n_files; i_file+=step) {

                h_a_track_length_mock[i_file] = pdfs[i_file]->createHistogram("track_length_obs",100);
                h_a_track_length_mock[i_file]->Scale(fitted_activity[i_file]->getValV());
                h_a_track_length_mock[i_file]->SetLineWidth(2);
                h_a_track_length_mock[i_file]->SetLineColor(i_file+2);
                h_a_track_length_mock[i_file]->SetFillColor(i_file+2);
                h_a_track_length_mock[i_file]->SetFillStyle(3003);
                h_a_track_length_mock[i_file]->Draw("histsames");
                h_a_track_length_mock[i_file]->SetTitle(pdf_explicit_names[i_file]);
                l_mock->AddEntry(h_a_track_length_mock[i_file],h_a_track_length_mock[i_file]->GetTitle(),"l");
                h_a_track_length_mock[i_file]->SetTitle("");
            }
            if(step == 1) h_a_track_length_mock[1]->Draw("histsames");
            l_mock->Draw();
            
            TCanvas *c_mock1 = new TCanvas("track_length_pseudo","track_length_pseudo",1961,344,700,502);
        
            RooPlot * pl = track_length_obs.frame(Title("pseudo track length distribution"));
            data->plotOn(pl);
            model.plotOn(pl);
            //pdfs[0]->plotOn(pl);
            //h_fitted_activities[0]->plotOn(pl);
            pl->SetTitle("");
            pl->Draw();

            for(int i_file=first_file; i_file<n_files; i_file+=step) {
                cout << activities[i_file]*(given_exposure*efficiencies[i_file])
                    <<" "<< fitted_activity[i_file]->getValV() << endl;
            }
        }
    }

    //IMPORTANT: fitted_activity[0]->getValV() is equivalent to h_a_track_length[0]->Integral()*given_exposure

    sprintf(name,"./fitted_ab_%s_sel.txt",type_of_sel.data());
    ofstream f_ab (name,ios::app);
    f_ab << days_of_exposure << " ";
        
    TCanvas *c_fa[n_files];
    for(int i_file=first_file;i_file<n_files;i_file+=step) {
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
        f_ab << (myfit->GetParameter(2)/myfit->GetParameter(1)) << " ";
    }
    f_ab << endl;
    f_ab.close();
    
    
}

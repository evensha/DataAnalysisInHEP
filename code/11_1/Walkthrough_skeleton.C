/*
 * Project:        Exercises 11.1-11.3
 * File:           Walkthrough_skeleton.C
 * Author:         Ivo van Vulpen, Aart Heijboer
 * Version (date): 1.0 (23.06.2013)
 *
 * Copyright (C) 2013, Ivo van Vulpen, Aart Heijboer
 * All rights reserved.
 *
 * Description:
 * A code skeleton for the searches part.
 *
 * This code is distributed with the solution manual to the book
 *
 * Data Analysis in High Energy Physics: A Practical Guide to Statistical Methods,
 * Wiley-VCH (2013),
 * O. Behnke, K. Kroeninger, G. Schott, Th. Schoerner-Sadenius (editors)
 */


#include "TBox.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/ProbFuncMathCore.h" // for ROOT::Math::gaussian_cdf

#include <iostream>
using namespace std;


//-------------------------------------------------------------------------------------------------------------------------
//-- full functions
TH1D * GetMassDistribution(int Itype = 1, double scalefactor = 1.00);
void MassPlot(int Irebin = 20);

//-- skeleton functions
void SideBandFit(int Irebin = 10, int print = 1);
double Get_TestStatistic(TH1D *h_mass_dataset, TH1D *h_template_bgr, TH1D *h_template_sig);
TH1D * GenerateToyDataSet(TH1D *h_mass_template);

//-- some fusefull functions
double IntegratePoissonFromRight(double mu, int N_obs);
double IntegrateFromRight(TH1D * h_X_bgr, double X_value);
vector<double> Get_Quantiles( TH1D* hist );
void AddText( Double_t txt_x = 0.50, Double_t txt_y = 0.50, const char * txt = "dummy", Double_t txt_size = 0.045, 
	      Double_t txt_angle = 0., const char * Alignment = "left", Int_t UseNormalizedSize = 1, Int_t txt_color =1 );

//-- my own functions 
void LumiScan(double max_lumi, double lumi_step); 
void MakePoissonHist(); 
//-------------------------------------------------------------------------------------------------------------------------



//========================================
// S O M E   F I N A L   F U N C T I O N S 
//========================================


//========================================================
TH1D * GetMassDistribution(int Itype, double scalefactor){
//========================================================
 //----------------------------------------------------------
 // Goal: return the histogram of the 4-lepton invariant mass
 //       for given type with an optional scale factor
 //
 //       Itype 1 = ZZ SM background
 //             2 = data
 //           125 = Higgs 125
 //           200 = Higgs 200
 //
 //      scalefactor: histograms will be scaled with this number
 //
 //  Note: Histograms have ~200 MeV bins, so need to rebin
 //---------------------------------------------------------

  //-- [1] Get histogram from the file
  TH1D *h_mass = 0;
  TDirectory* dir = gDirectory;   
  TFile *file = new TFile("Histograms_fake.root", "READ");
  dir->cd();

  //-- Higgs 125
  if(Itype == 125){
    h_mass  = (TH1D*) file->Get("h_m4l_Higgs125_fake")->Clone("h_mass");     
  }
  //-- Higgs 200
  if(Itype == 200){
    h_mass  = (TH1D*) file->Get("h_m4l_Higgs200_fake")->Clone("h_mass");     
  }
  //-- ZZ SM background
  if(Itype == 1){
    h_mass  = (TH1D*) file->Get("h_m4l_ZZ_fake")->Clone("h_mass");     
  }
  //-- data
  if(Itype == 2){
    h_mass  = (TH1D*) file->Get("h_m4l_data_fake")->Clone("h_mass");     
  }
  
  //-- [2] scale histograms
  int Nbins = h_mass->GetNbinsX();
  for (int i_bin = 1; i_bin < Nbins; i_bin++){
    double mu_bin = h_mass->GetBinContent(i_bin);
    h_mass -> SetBinContent( i_bin, scalefactor * mu_bin);
  }


  file->Close();
  //-- [3] return histogram   
  return h_mass;

  //===========================
} // end GetMassDistribution()
  //===========================




//========================
void MassPlot(int Irebin){
//========================
  // ------------------------------------------
  // Goal: produce SM+Higgs+data plot
  //       Note: rebinning is only for plotting
  // ------------------------------------------

  //------------------------------------
  //-- Standard stuff and prepare canvas
  //------------------------------------
  gROOT->Clear();
  gROOT->Delete();

  //-- Prepare canvas and plot histograms
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125); 
  canvas1->cd(); 

  //------------------------------------------------------------------
  //-- [1] Prepare histograms
  //--     o Get histograms from the files (signal, background and data)
  //--     o Make cumulative histograms (for signal and background)
  //------------------------------------------------------------------

  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  h_sig  = GetMassDistribution(125);
  h_bgr  = GetMassDistribution(1);
  h_data = GetMassDistribution(2);

  //-----------------------------------
  //-- [2] Plot histograms and make gif
  //--     o rebin histograms 
  //--     o prepare cumulative histogram
  //--     o make plot + opsmuk + gif
  //-----------------------------------

  //-- Rebin histograms (only for plotting)
  h_sig->Rebin(Irebin);
  h_bgr->Rebin(Irebin);
  h_data->Rebin(Irebin);

  //-- Prepare cumulative histogram for signal + background 
  TH1D *h_sig_plus_bgr = (TH1D* ) h_bgr->Clone("h_sig_plus_bgr");
  h_sig_plus_bgr->Reset();
  for (int i_bin = 1; i_bin < h_bgr->GetNbinsX(); i_bin++){
       h_sig_plus_bgr->SetBinContent( i_bin, h_sig->GetBinContent(i_bin) + h_bgr->GetBinContent(i_bin));
       //printf("  REBINNED HISTOGRAM:  bin %d, Ndata = %d\n",i_bin,(int)h_data->GetBinContent(i_bin));
  }

  //-- prepare histograms and plot them on canvas
  double Data_max = h_data->GetBinContent(h_data->GetMaximumBin());
  double Ymax_plot = 1.10* (Data_max + TMath::Sqrt(Data_max));
  h_sig_plus_bgr->SetFillColor(7); 
  h_sig_plus_bgr->SetAxisRange(0.,Ymax_plot,"Y");
  h_sig_plus_bgr->SetAxisRange(0.,400.,"X");
  h_bgr->SetFillColor(2); 
  h_sig_plus_bgr->Draw("hist");  
  h_bgr->Draw("same");  
  h_bgr->Draw("axis same");  
  h_data->Draw("e same");

  //-- some nice axes and add legend
  AddText( 0.900, 0.035, "4-muon invariant mass [GeV]",0.060, 0.,"right");                             // X-axis
  AddText( 0.040, 0.900, Form("Number of events / %3.1f GeV",h_bgr->GetBinWidth(1)) ,0.060,90.,"right"); // Y-axis
  TLegend *leg1 = new TLegend(0.65,0.65,0.90,0.85);
  leg1->SetBorderSize(0); leg1->SetFillColor(0);
  TLegendEntry *leg1a = leg1->AddEntry(h_bgr,          " SM(ZZ)", "f");  leg1a->SetTextSize(0.04);
  TLegendEntry *leg1b = leg1->AddEntry(h_sig_plus_bgr, " Higgs" , "f");  leg1b->SetTextSize(0.04);
  leg1->Draw();

  //-- prepare gif
  canvas1->Print(Form("./MassPlot_rebin%d.png",Irebin));
  canvas1->Close(); 
  
  return;

   //===============
 } // end MassPlot()
   //===============




//===============================================
// S O M E   S K E L E T O N    F U N C T I O N S 
//===============================================


//========================================================================================
double Significance_Optimization(double Lumi_scalefactor = 1.00, int print_and_plot = 0){
//========================================================================================

  printf("\n Significance_Optimization()\n\n");

  //------------------------------------------------------------------
  //-- [1] Prepare histograms
  //--     o Get histograms from the files (signal, background and data)
  //--     o scale to correct luminosity
  //------------------------------------------------------------------

  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  printf ("\n  INFO: Mass distribution in the 4 lepton channel\n");
  h_sig  = GetMassDistribution(125, Lumi_scalefactor);
  h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
  h_data = GetMassDistribution(2, Lumi_scalefactor);

  //-------------------------------------------------------
  //-- [2] Compute significance for various mass windows
  //--     o try various options for window (use histogram)
  //--     o compute expected and observed significance 
  //-------------------------------------------------------

  //-- Define histogram that defines mass windows to try
  TH1D *h_masswindow          = new TH1D("h_masswindow","",250,0.,25.);           // make a mass window - full width between 0 and 25 GeV
  TH1D *h_masswindow_expected = new TH1D("h_masswindow_expected","",250,0.,25.);  // histogram to hold results for expected
  TH1D *h_masswindow_observed = new TH1D("h_masswindow_observed","",250,0.,25.);  // histogram to hold results for observed

  double best_exp_masswindow = 0;
  double best_obs_masswindow = 0; 
  double best_exp_significance = 0; 
  double best_obs_significance = 0; 
  double n_sig_best_window = 0; 
  double n_bgr_best_window = 0;
  double n_data_best_window = 0; 

  //---------------------------------
  //-- Loop over various mass windows
  //---------------------------------
  for (int i_bin = 1; i_bin<h_masswindow->GetNbinsX(); i_bin++ ){
  
    //-- get full width of mass window (bin center of our histogram) and the number of events in mass window for each event type
    double masswindow_fullwidth = h_masswindow->GetBinCenter(i_bin);

    //-- [a] determine the number of events in the mass window for each event type
    //       Ndata_win, Nbgr_win and Nsig_win
    int bin_125 = h_data->GetXaxis()->FindBin(125);
    double Ndata_win = h_data->Integral(bin_125-masswindow_fullwidth/(2*0.2),bin_125+masswindow_fullwidth/(2*0.2));   
    double Nbgr_win = h_bgr->Integral(bin_125-masswindow_fullwidth/(2*0.2),bin_125+masswindow_fullwidth/(2*0.2)); 
    double Nsig_win = h_sig->Integral(bin_125-masswindow_fullwidth/(2*0.2),bin_125+masswindow_fullwidth/(2*0.2)); 

    if((Nbgr_win+Nsig_win)<1){continue;}  

    //-- [b] compute EXPECTED significance and save in histogram
    double pvalue_expected       = IntegratePoissonFromRight(Nbgr_win, Nbgr_win+Nsig_win); // you need to do this yourself
    double significance_expected = ROOT::Math::gaussian_quantile_c(pvalue_expected,1);
    h_masswindow_expected->SetBinContent(i_bin, significance_expected);    

    //-- [c] compute OBSERVED significance and save in histogram
    double pvalue_observed       = IntegratePoissonFromRight(Nbgr_win, Ndata_win); // you need to do this yourself
    double significance_observed = ROOT::Math::gaussian_quantile_c(pvalue_observed,1);
    h_masswindow_observed->SetBinContent(i_bin, significance_observed);

    if(significance_expected >= best_exp_significance){ 
      best_exp_significance = significance_expected; 
      best_exp_masswindow = masswindow_fullwidth;
      n_sig_best_window = Nsig_win; 
      n_bgr_best_window = Nbgr_win;
      n_data_best_window = Ndata_win; 
    } 
    if(significance_observed >= best_obs_significance){ best_obs_significance = significance_observed; best_obs_masswindow = masswindow_fullwidth;} 
  } // end loop over width mass window

  //-- print optimum to the screen
  if(print_and_plot == 1){
    cout << "-------------------------" << endl; 
    cout << "Best expected mass window: " << best_exp_masswindow << endl;
    cout << "Best expected significance: " <<  best_exp_significance << endl;
    cout << "Number of background events: " << n_bgr_best_window << endl; 
    cout << "Number of signal events: " << n_sig_best_window << endl; 
    cout << "Number of data events: " << n_data_best_window << endl; 
    cout << "Best observed mass window: " <<  best_obs_masswindow << endl; 
    cout << "Best observed significance: " << best_obs_significance << endl; 
    cout << "-------------------------" << endl;   

    //----------------------------------
    //-- [3] Plot histogram and make gif
    //----------------------------------
    TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
    canvas1->SetLeftMargin(0.125);
    canvas1->SetBottomMargin(0.125); 
    canvas1->cd(); 
 
    h_masswindow_expected->SetLineColor(1);
    h_masswindow_expected->SetLineWidth(2);
    h_masswindow_observed->SetLineColor(4);
    h_masswindow_observed->SetLineWidth(2);

    h_masswindow_expected->SetAxisRange(-1.,6.,"Y");
    h_masswindow_expected->Draw("l");
    if(fabs(Lumi_scalefactor-1.00)<0.01){
      h_masswindow_observed->Draw("l same");
    }
    //-- axes
    AddText( 0.900, 0.035, "Mass window (GeV)",0.060, 0.,"right"); // X-axis
    AddText( 0.040, 0.900, "Significance" ,0.060,90.,"right");   // Y-axis    
    //AddText( 0.225, 0.825, Form("Luminosity scalefactor = %5.1f",Lumi_scalefactor),0.050, 0.,"left");            

    AddText( 0.700, 0.200, "Expected significance",0.050, 0.,"right",1,1);                        
    if(fabs(Lumi_scalefactor-1.00)<0.01){
      AddText( 0.700, 0.300, "Observed significance",0.050, 0.,"right",1,4);                        
    }
    //-- prepare gif
    canvas1->Print(Form("./Significance_Optimization_lumiscalefactor%d.png",int(Lumi_scalefactor)));
    canvas1->Close(); 
  }

  //-- Prepare canvas and plot histograms
  TCanvas * canvas2 = new TCanvas( "canvas2","Standard Canvas",600,400);  
  canvas2->SetLeftMargin(0.125);
  canvas2->SetBottomMargin(0.125); 
  canvas2->cd(); 

  //-- Rebin histograms (only for plotting)
  h_sig->Rebin(20);
  h_bgr->Rebin(20);
  h_data->Rebin(20);

  //-- Prepare cumulative histogram for signal + background 
  TH1D *h_sig_plus_bgr = (TH1D* ) h_bgr->Clone("h_sig_plus_bgr");
  h_sig_plus_bgr->Reset();
  for (int i_bin = 1; i_bin < h_bgr->GetNbinsX(); i_bin++){
       h_sig_plus_bgr->SetBinContent( i_bin, h_sig->GetBinContent(i_bin) + h_bgr->GetBinContent(i_bin));
       //printf("  REBINNED HISTOGRAM:  bin %d, Ndata = %d\n",i_bin,(int)h_data->GetBinContent(i_bin));
  }

  //-- prepare histograms and plot them on canvas
  double Data_max = h_data->GetBinContent(h_data->GetMaximumBin());
  double Ymax_plot = 1.10* (Data_max + TMath::Sqrt(Data_max));
  h_sig_plus_bgr->SetFillColor(7); 
  h_sig_plus_bgr->SetAxisRange(0.,Ymax_plot,"Y");
  h_sig_plus_bgr->SetAxisRange(110.,140.,"X");
  h_bgr->SetFillColor(2); 
  h_sig_plus_bgr->Draw("hist");  
  h_bgr->Draw("same");  
  h_bgr->Draw("axis same");  
  h_data->Draw("e same");

  //-- some nice axes and add legend
  AddText( 0.900, 0.035, "4-muon invariant mass [GeV]",0.060, 0.,"right");                             // X-axis
  AddText( 0.040, 0.900, Form("Number of events / %3.1f GeV",h_bgr->GetBinWidth(1)) ,0.060,90.,"right"); // Y-axis
  TLegend *leg1 = new TLegend(0.65,0.65,0.90,0.85);
  leg1->SetBorderSize(0); leg1->SetFillColor(0);
  TLegendEntry *leg1a = leg1->AddEntry(h_bgr,          " SM(ZZ)", "f");  leg1a->SetTextSize(0.04);
  TLegendEntry *leg1b = leg1->AddEntry(h_sig_plus_bgr, " Higgs" , "f");  leg1b->SetTextSize(0.04);
  leg1->Draw();

  TLine *l1 = new TLine(125-3.575, 0, 125-3.575, Ymax_plot);  
  TLine *l2 = new TLine(125+3.575, 0, 125+3.575, Ymax_plot); 
  l1->Draw(); l2->Draw(); 
  //-- prepare gif
  canvas2->Print("./MassPlot_optimized_window.png");
  canvas2->Close(); 

  
  h_masswindow->Delete(); 
  h_masswindow_expected->Delete(); 
  h_masswindow_observed->Delete(); 

  return best_exp_significance;

  //================================
} // end Significance_Optimization()
  //================================






//===========================================
void SideBandFit(int Irebin, int print = 0){
//===========================================

  printf("\n SideBandFit()\n\n");

  //-------------------------
  //-- [1] Prepare histograms
  //-------------------------
  TH1D *h_bgr, *h_data;
  h_bgr  = GetMassDistribution(1);
  h_data = GetMassDistribution(2);
 
  //-- rebin histograms if necessary
  h_bgr->Rebin(Irebin);
  h_data->Rebin(Irebin);
  printf(" INFO: Rebinning the histograms with a factor %d. Binwidth is now %5.2f GeV\n", Irebin, h_data->GetBinWidth(1));

  //-----------------------------------------
  //-- [2] Loop over scale factor (alpha_bgr)
  //-----------------------------------------
  TH1D *h_scalefactor_bgr = new TH1D("h_scalefactor_bgr","",300,0.1,3.1); // you'll want to put more here I guess
  double scalefactor_bgr  = 1.00; 
  double min_loglik = 1E10; 
  double best_sf_bgr = 0; 

  int lower_bin = h_data->GetXaxis()->FindBin(150); 
  int upper_bin = h_data->GetXaxis()->FindBin(400); 
  //cout << lower_bin << endl; 
  //cout << upper_bin << endl; 

  //---------------------------------------------
  // [2a] Loop 1: loop over scalefactors in alpha
  //---------------------------------------------
  for (int i_bin_sf = 1; i_bin_sf <=h_scalefactor_bgr->GetNbinsX(); i_bin_sf++){

    //-- determine the scale factor for the background
    scalefactor_bgr = h_scalefactor_bgr->GetBinLowEdge(i_bin_sf);
    //printf(" Loop 1: I am now trying alpha = %5.2f\n",scalefactor_bgr);
  
    //-----------------------------------------------------------------------------------
    // [2b] Loop 2: loop over bins in the histogram, compute loglik and save in histogram
    //-----------------------------------------------------------------------------------
    double loglik = 0;
    for (int i_bin = lower_bin; i_bin <= upper_bin; i_bin++){
      //printf("        bin %d, m = %5.2f, Ndata = %5.2f, Nbgr = %5.2f\n",i_bin, h_bgr->GetBinCenter(i_bin),h_bgr->GetBinContent(i_bin),h_data->GetBinContent(i_bin));
      loglik +=TMath::Log(TMath::Poisson(h_data->GetBinContent(i_bin), h_bgr->GetBinContent(i_bin)*scalefactor_bgr)); // you'll have to do this yourself;  	
    } // end loop over bins

    if(-2*loglik < min_loglik){ min_loglik = -2*loglik; best_sf_bgr = scalefactor_bgr; }
    h_scalefactor_bgr->SetBinContent(i_bin_sf,-2.*loglik);   
  } // end loop over scale factors for the background

  //----------------------------------------------------
  //-- [3] Interpret the -2Log (Likelihood distribution)
  //----------------------------------------------------

  double lh_unc =  h_scalefactor_bgr->GetBinContent(h_scalefactor_bgr->FindBin(best_sf_bgr))+1; 

  int bin_min = h_scalefactor_bgr->GetXaxis()->FindBin(best_sf_bgr); 
  int bin_max = h_scalefactor_bgr->GetXaxis()->FindBin(best_sf_bgr); 
  double lh_min = 0;
  double lh_max = 0; 

  while(lh_min<lh_unc){
    lh_min = h_scalefactor_bgr->GetBinContent(bin_min);
    bin_min--; 
  }
  while(lh_max<lh_unc){
    lh_max = h_scalefactor_bgr->GetBinContent(bin_max);
    bin_max++; 
  }

  double delta_p = h_scalefactor_bgr->GetXaxis()->GetBinLowEdge(bin_max) - best_sf_bgr; 
  double delta_m = h_scalefactor_bgr->GetXaxis()->GetBinLowEdge(bin_min) - best_sf_bgr; 
  if(print == 1){
    cout << "-------------------------------------" << endl; 
    cout << "Best background scalefactor: " << best_sf_bgr << endl; 
    cout << "Lower uncertainty in alpha: " << delta_m  << endl; 
    cout << "Upper uncertainty in alpha: " << delta_p << endl;
  }

  TH1D *h_scalefactor_bgr_rescaled  = (TH1D*) h_scalefactor_bgr->Clone("h_scalefactor_bgr_rescaled");  
  h_scalefactor_bgr_rescaled->SetAxisRange(0.8, 1.4, "X"); 
  //h_scalefactor_bgr_rescaled->Reset();   
    
  //--Rescale and find +/- 1 sigma errors
  //double unscaled_events = (h_bgr->Integral(h_bgr->GetXaxis()->FindBin(125-5),h_bgr->GetXaxis()->FindBin(125+5)));
  double unscaled_events = h_bgr->Integral(h_bgr->GetXaxis()->FindBin(125-7.15/2.0),h_bgr->GetXaxis()->FindBin(125+7.15/2.0));
  double scaled_events = best_sf_bgr*unscaled_events;

  if(print == 1){
    cout << "Events in mass window: " << endl; 
    cout << "Unscaled background events: " << unscaled_events << endl; 
    cout << "Scaled background events: " << scaled_events << "  +" << scaled_events/best_sf_bgr*delta_p << "  -" << fabs(scaled_events/best_sf_bgr*delta_m) << endl; 
  }

  //-- print summary to screen

  //----------------------------------
  //-- Plot histogram and make png
  //----------------------------------
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125); 
  canvas1->cd(); 
  
  h_scalefactor_bgr_rescaled->SetLineColor(kBlue);
    
  //h_masswindow_expected->SetAxisRange(-1.,6.,"Y");
  h_scalefactor_bgr_rescaled->Draw("l");

  //-- axes
  AddText( 0.900, 0.035, "#alpha_{bgr}",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "-2ln#it{L}" ,0.060,90.,"right");   // Y-axis    
  if(fabs(delta_m) != delta_p){
    AddText( 0.325, 0.725, Form("#hat{#alpha}_{bgr} = %5.2f{}^{+%.2f}_{-%.2f}",best_sf_bgr,delta_p,fabs(delta_m)),0.080, 0.,"left");            
  }
  else{
    AddText( 0.325, 0.725, Form("#hat{#alpha}_{bgr} = %5.2f #pm %.2f",best_sf_bgr, delta_p),0.080, 0.,"left");            
  }

  
  //-- prepare gif
  canvas1->Print("./Loglikelihood_sf_bgr.png");
  canvas1->Close(); 
    

  return;

  //==================
} // end SideBandFit()
  //==================



//=========================================================================================
double Get_TestStatistic(TH1D *h_mass_dataset, TH1D *h_template_bgr, TH1D *h_template_sig){
//=========================================================================================
   
  double n_bgr = 0; 
  double n_sig = 0; 
  double n_data = 0; 
  double alpha_bgr = 1.0;
  double loglik_bgr = 0; 
  double loglik_sig_plus_bgr = 0; 
  int lower_bin = 1.0; 
  int upper_bin = h_template_bgr->GetNbinsX(); 

  //-- do likelihood fit
  for (int i_bin = lower_bin; i_bin <= upper_bin; i_bin++){
    n_bgr = h_template_bgr->GetBinContent(i_bin); 
    n_sig = h_template_sig->GetBinContent(i_bin); 
    n_data = h_mass_dataset->GetBinContent(i_bin); 
    loglik_bgr += TMath::Log(TMath::Poisson(n_data, n_bgr*alpha_bgr)); // mu=0 (b-only) 
    loglik_sig_plus_bgr += TMath::Log(TMath::Poisson(n_data, n_sig+n_bgr*alpha_bgr)); // mu=1 (s+b) 
  } // end loop over bins
 
  //-- compute test statistic
  double test_statistic  =  -2*loglik_sig_plus_bgr+2*loglik_bgr; 

  //-- return test_statistic
  return test_statistic;

} // end Get_TestStatistic()





//===============================================
TH1D * GenerateToyDataSet(TH1D *h_mass_template){
//===============================================
  //-------------------------------------------------
  // Goal: Generate Toy data set from input histogram
  // How:  dumb way -> draw random Poisson in each bin       
  //--------------------------------------------------
  TRandom3 *R = new TRandom3(0);

  //-- Create new histogram for the data-set
  TH1D *h_mass_toydataset = (TH1D*) h_mass_template->Clone("h_mass_toydataset"); h_mass_toydataset->Reset();
 
  //-- Loop over bins and draw Poisson number of event in each bin
  double expected; 
  double poisson; 
  for(int bin_i = 1; bin_i <= h_mass_template->GetNbinsX(); bin_i++){
    expected = h_mass_template->GetBinContent(bin_i); 
    poisson = R->Poisson(expected);
    h_mass_toydataset->SetBinContent(bin_i, poisson); 
  }

  R->Delete();   
  //-- return histogram of toy data-set
  return h_mass_toydataset;
  
} // end GenerateToyDataSet()












//================================================
// S O M E   U S E F U L L   F U N C T I O N S 
//================================================



//====================================================
double IntegratePoissonFromRight(double mu, int N_obs){
//====================================================
// --------------------------------------------------------------------
// Compute p-value for case zero background uncertainty, i.e. 
//         just integrate Poisson from the right from N_obs to infinity
// --------------------------------------------------------------------

  double integral = 1.; 
  for(int i_obs = 0; i_obs < N_obs; i_obs++){
    integral -= TMath::Poisson(i_obs,mu);
  }
  
  return integral;
  
} // end IntegratePoissonFromRight()


//========================================================
double IntegrateFromRight(TH1D * h_X_bgr, double X_value){
//========================================================
// --------------------------------------------------------------------
// Compute p-value: integrate number of events from X_value to infinity
// --------------------------------------------------------------------

  //-- Integrate distributions
  int Nbins = h_X_bgr->GetNbinsX();
  int X_bin = h_X_bgr->FindBin(X_value); 

  //-- Compute integral from X-value to infinity
  double pvalue = h_X_bgr->Integral(X_bin,Nbins) / h_X_bgr->Integral();
  
  return pvalue;

} // end IntegrateFrom Right()




//=========================================
vector<double> Get_Quantiles( TH1D* hist ){
//=========================================
// Quantiles returns a vector<double> with 5 entries.
// Entries 0 and 4 are the values on the histogram x-axis
// so that 95% of the content lies between these values.
// Entries 1 and 3 bracket 68% in the same way.
// Entry 2 is the median of the histogram.

  //-- define quantiles
  double fraction_1sigma = ROOT::Math::gaussian_cdf(-1.,1.,0.); // 15.8655 %
  double fraction_2sigma = ROOT::Math::gaussian_cdf(-2.,1.,0.); //  2.2750 %
  double probs[5] = {fraction_2sigma, fraction_1sigma, 0.50, 1.00-fraction_1sigma, 1.00-fraction_2sigma };

  //-- output of the quantiles
  double Xvalues[5];

  //-- extract quantiles
  hist->GetQuantiles( 5, Xvalues, probs );
  
  vector<double> Xvalues_output(5);
  for (int i=0; i<5; i++) 
    {
      Xvalues_output[i] = Xvalues[i];
    }

  return Xvalues_output;
} // end Get_Quantiles()





//=======================================================================================================================
void AddText( Double_t txt_x, Double_t txt_y, const char * txt, Double_t txt_size,
              Double_t txt_angle, const char * Alignment, Int_t UseNormalizedSize, Int_t txt_color)
//=======================================================================================================================
{
  Int_t txt_align = 12;
  if ( !strcmp(Alignment, "left"))   { txt_align = 12; } // left 
  if ( !strcmp(Alignment, "right"))  { txt_align = 32; } // right
  if ( !strcmp(Alignment, "center")) { txt_align = 22; } // center

  TLatex* t1 = new TLatex( txt_x, txt_y, txt);
  if(UseNormalizedSize) {t1->SetNDC(kTRUE);} // <- use NDC coordinate
  t1->SetTextSize(txt_size);
  t1->SetTextAlign(txt_align);
  t1->SetTextAngle(txt_angle);
  t1->SetTextColor(txt_color);
  t1->Draw();

} // end AddText()


//===============
void LumiScan(){
//===============

  double lumi;
  double sign;
  TH1D *h_lumi = new TH1D("h_lumi", "", 60, 1, 7);

  for(int i = 1; i<=h_lumi->GetNbinsX(); i++){
    lumi = h_lumi->GetXaxis()->GetBinLowEdge(i); 
    sign = Significance_Optimization(lumi);
    h_lumi->SetBinContent(i, sign); 
  }  

  cout << h_lumi->GetBinContent(h_lumi->FindFirstBinAbove(5)) << endl; 
  cout << h_lumi->GetBinLowEdge(h_lumi->FindFirstBinAbove(5)) << endl; 
  
  //----------------------------------
  //-- Plot histogram and make gif
  //----------------------------------
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125); 
  canvas1->cd(); 
 
  h_lumi->SetLineColor(kBlue);
  h_lumi->Draw("l");

  canvas1->Update(); 

  TLine *line = new TLine(1,5,7,5);
  line->SetLineColor(kRed);
  line->Draw(); 
  
  //-- axes
  AddText( 0.900, 0.035, "Integrated luminosity [fb{}^{-1}]",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Significance [#sigma]" ,0.060,90.,"right");   // Y-axis    
  //AddText( 0.325, 0.725, Form("#hat{#alpha}_{bgr} = %5.2f{}^{+ 0.07}_{- 0.06}",best_sf_bgr),0.080, 0.,"left");            

  //-- prepare gif
  canvas1->Print("./Significance_vs_luminosity.png");
  canvas1->Close(); 

  
  return; 
}


//=========================================================================================================================
void ExpectedSignificance_ToyMC(double exp_bgr=4.64, double exp_sig=5.41, double delta_bgr=0.0, int obs_data=13, int N=1E6)
//=========================================================================================================================
{ 

  TRandom3 *R = new TRandom3(0);

  TH1D *h_poisson_b = new TH1D("Poisson dist b-only", "", 25, 0, 25); 
  TH1D *h_poisson_sb = new TH1D("Poisson dist s+b", "", 25, 0, 25); 

  for(int i = 1; i<N; i++){
    
    h_poisson_b->Fill(R->Poisson(R->Gaus(exp_bgr,delta_bgr))); 
    h_poisson_sb->Fill(R->Poisson(R->Gaus(exp_bgr,delta_bgr)+exp_sig)); 

  }

  double expected_pvalue = (h_poisson_b->Integral(h_poisson_b->GetXaxis()->FindBin(exp_bgr+exp_sig), h_poisson_b->GetNbinsX()))/(h_poisson_b->Integral(1, h_poisson_b->GetNbinsX())); 
  double expected_significance = ROOT::Math::gaussian_quantile_c(expected_pvalue,1);
  double observed_pvalue = (h_poisson_b->Integral(h_poisson_b->GetXaxis()->FindBin(obs_data), h_poisson_b->GetNbinsX()))/(h_poisson_b->Integral(1, h_poisson_b->GetNbinsX())); 
  double observed_significance = ROOT::Math::gaussian_quantile_c(observed_pvalue,1); 
  cout << "---------------------------" << endl; 
  cout << "Expected p-value: " << expected_pvalue << endl; 
  cout << "Expected significance: " << expected_significance << endl; 
  cout << "Observed p-value: " << observed_pvalue << endl; 
  cout << "Observed significance: " << observed_significance << endl; 

  //----------------------------------
  //-- Plot histogram and make gif
  //----------------------------------
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125); 
  canvas1->cd(); 

  h_poisson_b->SetLineColor(kRed); 
  h_poisson_sb->SetLineColor(kBlue);

  TLegend *leg = new TLegend(0.70, 0.70, 0.88, 0.88); 
  leg->SetFillStyle(4000); 
  leg->SetFillColor(0); 
  leg->SetBorderSize(0); 
  leg->SetTextFont(22); 

  
  //h_masswindow_expected->SetAxisRange(-1.,6.,"Y");
  h_poisson_b->Draw("hist");
  h_poisson_sb->Draw("same hist"); 
  
  //-- axes
  AddText( 0.900, 0.035, "Number of events",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Number of experiments" ,0.060,90.,"right");   // Y-axis    
  AddText( 0.425, 0.725, "Observed events",0.050, 0.,"left");            
  
  leg->AddEntry(h_poisson_b, "b-only", "l"); 
  leg->AddEntry(h_poisson_sb, "s+b", "l"); 
  leg->Draw();  

  TArrow *ar1 = new TArrow(obs_data,0.05*h_poisson_b->GetMaximum(),obs_data,0.75*h_poisson_b->GetMaximum(), 0.02, "<|");
  ar1->Draw();
    
  //-- prepare gif
  canvas1->Print("./ToyExperiment_Nevents.png");
  canvas1->Close(); 

  h_poisson_b->Delete(); 
  h_poisson_sb->Delete(); 

  return; 
}

//============================================================================================
void Significance_LikelihoodRatio_ToyMC(int N, double signal_sf = 1.0, double lumi_sf=2.0 )
//============================================================================================
{

  TH1D *h_bgr, *h_data, *h_sig, *h_toy_bgr_only, *h_toy_sig_plus_bgr;
  h_bgr  = GetMassDistribution(1, lumi_sf); //h_bgr->Rebin(rebin); 
  h_data = GetMassDistribution(2, lumi_sf); //h_data->Rebin(rebin); 
  h_sig = GetMassDistribution(125, lumi_sf); //h_sig->Rebin(rebin); 
  h_sig->Scale(signal_sf); 

  TH1D *h_sig_plus_bgr = (TH1D* ) h_bgr->Clone("h_sig_plus_bgr");
  h_sig_plus_bgr->Reset(); 
  for (int i_bin = 1; i_bin < h_bgr->GetNbinsX(); i_bin++){
       h_sig_plus_bgr->SetBinContent( i_bin, h_sig->GetBinContent(i_bin) + h_bgr->GetBinContent(i_bin));
  }


  TH1D *h_test_bgr_only = new TH1D("test bgr only", "", 275, -35, 20);
  TH1D *h_test_sig_plus_bgr = new TH1D("test bgr only", "", 275, -35, 20);

  if(signal_sf != 1.0 || lumi_sf != 1.0 ){ 
    h_test_bgr_only->SetBins(1000, -100, 100); 
    h_test_sig_plus_bgr->SetBins(1000, -100, 100); 
  }

  double test_stat_bgr_only;
  double test_stat_sig_plus_bgr; 
  double test_stat_data = Get_TestStatistic(h_data, h_bgr, h_sig); 
  cout << "Data test statistic: " <<  test_stat_data << endl; 
  
  for(int i = 1; i<N; i++){
    h_toy_bgr_only = GenerateToyDataSet(h_bgr); 
    h_toy_sig_plus_bgr = GenerateToyDataSet(h_sig_plus_bgr); 
    test_stat_bgr_only = Get_TestStatistic(h_toy_bgr_only, h_bgr, h_sig);     
    test_stat_sig_plus_bgr = Get_TestStatistic(h_toy_sig_plus_bgr, h_bgr, h_sig);     
    h_test_bgr_only->Fill(test_stat_bgr_only); 
    h_test_sig_plus_bgr->Fill(test_stat_sig_plus_bgr); 
    if(i%1000 == 0){ cout << i << " toy experiments done!" << endl;} 
  }


  double mean_bgr_only = h_test_bgr_only->GetMean(); 
  double mean_sig_plus_bgr = h_test_sig_plus_bgr->GetMean(); 

  double x_bgr_only, x_sig_plus_bgr , q; 
  q = 0.5; 
  h_test_bgr_only->GetQuantiles(1, &x_bgr_only, &q);  
  h_test_sig_plus_bgr->GetQuantiles(1, &x_sig_plus_bgr, &q); 

  cout << "Median t-value of b-only exp: " << x_bgr_only << endl; 
  cout << "Median t-value of s+b exp: " << x_sig_plus_bgr << endl; 

  double CLb_bgr_only, CLb_sig_plus_bgr, CLb_data; 
  CLb_bgr_only = h_test_bgr_only->Integral(1,h_test_bgr_only->FindBin(x_bgr_only))/h_test_bgr_only->Integral(1,h_test_bgr_only->GetNbinsX()); 
  CLb_sig_plus_bgr = h_test_bgr_only->Integral(1,h_test_bgr_only->FindBin(x_sig_plus_bgr))/h_test_bgr_only->Integral(1,h_test_bgr_only->GetNbinsX()); 
  CLb_data = h_test_bgr_only->Integral(1,h_test_bgr_only->FindBin(test_stat_data))/h_test_bgr_only->Integral(1,h_test_bgr_only->GetNbinsX()); 

  cout << "---------------------------------------------" << endl; 
  cout << "Discovery analysis:" << endl; 
  cout << "1-CL_b for b-only: " << CLb_bgr_only << endl; 
  cout << "1-CL_b for s+b: " << CLb_sig_plus_bgr << endl; 
  cout << "1-CL_b for data: " << CLb_data << endl; 

  double Zb_bgr_only = ROOT::Math::gaussian_quantile_c(CLb_bgr_only,1);  
  double Zb_sig_plus_bgr = ROOT::Math::gaussian_quantile_c(CLb_sig_plus_bgr,1);  
  double Zb_data = ROOT::Math::gaussian_quantile_c(CLb_data,1);

  cout << "Significance b-only: " << Zb_bgr_only << endl; 
  cout << "Significance s+b: " << Zb_sig_plus_bgr << endl; 
  cout << "Significance data: " << Zb_data << endl; 

  double CLsb_bgr_only, CLsb_sig_plus_bgr, CLsb_data; 
  CLsb_bgr_only = h_test_sig_plus_bgr->Integral(h_test_sig_plus_bgr->FindBin(x_bgr_only),h_test_sig_plus_bgr->GetNbinsX())/h_test_sig_plus_bgr->Integral(1,h_test_sig_plus_bgr->GetNbinsX()); 
  CLsb_sig_plus_bgr = h_test_sig_plus_bgr->Integral(h_test_sig_plus_bgr->FindBin(x_sig_plus_bgr),h_test_sig_plus_bgr->GetNbinsX())/h_test_sig_plus_bgr->Integral(1,h_test_sig_plus_bgr->GetNbinsX()); 
  CLsb_data = h_test_sig_plus_bgr->Integral(h_test_sig_plus_bgr->FindBin(test_stat_data),h_test_sig_plus_bgr->GetNbinsX())/h_test_sig_plus_bgr->Integral(1,h_test_sig_plus_bgr->GetNbinsX()); 

  cout << "---------------------------------------------" << endl; 
  cout << "Exclusion analysis:" << endl; 
  cout << "1-CL_s+b for b-only: " << CLsb_bgr_only << endl; 
  cout << "1-CL_s+b for s+b: " << CLsb_sig_plus_bgr << endl; 
  cout << "1-CL_s+b for data: " << CLsb_data << endl; 


    
  //----------------------------------
  //-- Plot histogram and make gif
  //----------------------------------
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125); 
  canvas1->cd(); 

  h_test_bgr_only->SetLineColor(kRed);
  h_test_sig_plus_bgr->SetLineColor(kBlue); 
  
  TLegend *leg = new TLegend(0.15, 0.70, 0.33, 0.88); 
  leg->SetFillStyle(4000); 
  leg->SetFillColor(0); 
  leg->SetBorderSize(0); 
  leg->SetTextFont(22);   
  

  //h_masswindow_expected->SetAxisRange(-1.,6.,"Y");
  h_test_bgr_only->Draw("hist");
  h_test_sig_plus_bgr->Draw("same hist"); 

  //-- axes
  AddText( 0.900, 0.035, "t=-2ln#it{Q}",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Number of experiments" ,0.060,90.,"right");   // Y-axis    
  AddText( 0.375, 0.675, "t_{data}" ,0.080, 0.,"left");            

  leg->AddEntry(h_test_bgr_only, "b-only", "l"); 
  leg->AddEntry(h_test_sig_plus_bgr, "s+b", "l"); 
  leg->Draw();  

  canvas1->Update(); 
  double x_max = h_test_bgr_only->GetMaximum();

  TArrow *ar1 = new TArrow(-11.5,0.05*x_max,-11.5,0.75*x_max, 0.02, "<|");
  ar1->Draw();

  //-- prepare gif
  if( signal_sf == 1.0 && lumi_sf == 1.0){ 
    canvas1->Print("./ToyExperiments_test_stat.png");
  }
  canvas1->Close(); 

  h_test_bgr_only->Delete(); 
  h_test_sig_plus_bgr->Delete(); 

  return; 
}


//===================
void MuFit(int rebin)
//===================
{

  //-------------------------
  //-- [1] Prepare histograms
  //-------------------------
  TH1D *h_bgr, *h_sig, *h_data;
  h_bgr  = GetMassDistribution(1);
  h_data = GetMassDistribution(2);
  h_sig = GetMassDistribution(125); 

  //-- rebin histograms if necessary
  h_bgr->Rebin(rebin);
  h_sig->Rebin(rebin); 
  h_data->Rebin(rebin);

  cout << h_sig->Integral(1,h_sig->GetNbinsX()) << endl; 
  
  //-----------------------------------------
  //-- [2] Loop over scale factor (alpha_bgr)
  //-----------------------------------------
  //TH1D *h_scalefactor_bgr = new TH1D("h_scalefactor_bgr","",300,0.1,3.1); 
  //TH1D *h_scalefactor_sig = new TH1D("h_scalefactor_sig","",300,0.1,3.1); 
  TH2D *h_scalefactors = new TH2D("h_scalefactors", "", 300, 0.5, 2.0, 300, 0.5, 2.0); // x = background, y = signal 
  double scalefactor_bgr  = 1.00; 
  double scalefactor_sig = 1.00; 
  double min_loglik = 1E10; 
  double best_sf_bgr = 0; 
  double best_sf_sig = 0; 
  double n_obs = 0; 
  double n_exp = 0; 

  //---------------------------------------------
  // [2a] Loop 1: loop over scalefactors
  //---------------------------------------------
  for (int i = 1; i <=h_scalefactors->GetNbinsX(); i++){
    for(int j = 1; j <= h_scalefactors->GetNbinsY(); j++){

      //-- determine the scale factor for the background
      scalefactor_bgr = h_scalefactors->GetXaxis()->GetBinLowEdge(i);
      scalefactor_sig = h_scalefactors->GetYaxis()->GetBinLowEdge(j); 
      //printf(" Loop 1: I am now trying alpha = %5.2f\n",scalefactor_bgr);
  
      //-----------------------------------------------------------------------------------
      // [2b] Loop 2: loop over bins in the histogram, compute loglik and save in histogram
      //-----------------------------------------------------------------------------------
      double loglik = 0;
      for (int i_bin = 1; i_bin <= h_bgr->GetNbinsX(); i_bin++){
	n_obs = h_data->GetBinContent(i_bin); 
	n_exp = h_bgr->GetBinContent(i_bin)*scalefactor_bgr + h_sig->GetBinContent(i_bin)*scalefactor_sig; 
	loglik +=TMath::Log(TMath::Poisson(n_obs, n_exp));   	
      } // end loop over bins

      if(-2*loglik < min_loglik){ min_loglik = -2*loglik; best_sf_bgr = scalefactor_bgr; best_sf_sig = scalefactor_sig;}
      h_scalefactors->SetBinContent(i, j, -2.*loglik);
   
    }
  } // end loop over scale factors for the background
  

  cout << "Best background scalefactor: " << best_sf_bgr << endl; 
  cout << "Best signal scalefactor: " << best_sf_sig << endl; 
  double lh_unc =  h_scalefactors->GetBinContent(h_scalefactors->FindBin(best_sf_bgr, best_sf_sig))+1; 

  int x_bin = h_scalefactors->GetXaxis()->FindBin(best_sf_bgr);
  int y_bin_min = h_scalefactors->GetYaxis()->FindBin(best_sf_sig); 
  int y_bin_max = h_scalefactors->GetYaxis()->FindBin(best_sf_sig); 
  double lh_min = 0;
  double lh_max = 0; 

  while(lh_min<lh_unc){
    lh_min = h_scalefactors->GetBinContent(x_bin, y_bin_min);
    y_bin_min--; 
  }
  while(lh_max<lh_unc){
    lh_max = h_scalefactors->GetBinContent(x_bin, y_bin_max);
    y_bin_max++; 
  }

  cout << "Lower uncertainty in mu: " << h_scalefactors->GetYaxis()->GetBinLowEdge(y_bin_min) - best_sf_sig << endl; 
  cout << "Upper uncertainty in mu: " << h_scalefactors->GetYaxis()->GetBinLowEdge(y_bin_max) - best_sf_sig << endl;

  TH2D *h_scalefactors_rescaled  = (TH2D*) h_scalefactors->Clone("h_scalefactors_rescaled");  
  h_scalefactors_rescaled->SetAxisRange(1.0, 1.2, "X"); 
  h_scalefactors_rescaled->SetAxisRange(0.8, 1.8, "Y"); 


  //-------------------
  // Plot 2D likelihood 
  //-------------------

  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125); 
  canvas1->cd(); 

  h_scalefactors_rescaled->Draw("COLZ");

  //-- axes
  AddText( 0.900, 0.035, "#alpha_{bgr}",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "#mu_{s}" ,0.060,90.,"right");   // Y-axis    

  //-- prepare gif
  canvas1->Print("./Loglikelihood_sig_and_bgr_sf.png");
  canvas1->Close(); 


  h_scalefactors->Delete(); 
  return; 

}

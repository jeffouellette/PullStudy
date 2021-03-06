#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TF1.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <math.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

#include <GlobalParams.h>
#include <Utilities.h>

#include "header.h"

using namespace atlashi;
using namespace pullstudy;

TH1D* h_trk_pth_pull[3][3][6];
TH1D* h_trk_xhz_pull[3][3][6];

TH1D* h_trk_pth_yield[3][3];
TH1D* h_trk_xhz_yield[3][3];
TH2D* h2_trk_pth_cov[3][3];
TH2D* h2_trk_xhz_cov[3][3];

TGraph* g_pth_yield_pull[3][3][3][6];
TGraph* g_xhz_yield_pull[3][3][3][6];

TGraphErrors* g_pth_pullWidth[3][3];
TGraphErrors* g_xhz_pullWidth[3][3];

int main () {

  SetAtlasStyle();

  TFile* histFile = new TFile ("rootFiles/outFile.root", "update");

  for (short iPtZ : {2, 3, 4}) {
    for (short iCent : {1, 2, 3}) {

      h_trk_pth_yield[iPtZ-2][iCent-1] = (TH1D*) histFile->Get (Form ("h_trk_pth_yield_iPtZ%i_iCent%i", iPtZ, iCent));
      h_trk_xhz_yield[iPtZ-2][iCent-1] = (TH1D*) histFile->Get (Form ("h_trk_xhz_yield_iPtZ%i_iCent%i", iPtZ, iCent));
      h2_trk_pth_cov[iPtZ-2][iCent-1] = (TH2D*) histFile->Get (Form ("h2_trk_pth_cov_iPtZ%i_iCent%i", iPtZ, iCent));
      h2_trk_xhz_cov[iPtZ-2][iCent-1] = (TH2D*) histFile->Get (Form ("h2_trk_xhz_cov_iPtZ%i_iCent%i", iPtZ, iCent));

      for (short iX = 0; iX < 6; iX++) {
        h_trk_pth_pull[iPtZ-2][iCent-1][iX] = (TH1D*) histFile->Get (Form ("h_trk_pth_pull_iPtZ%i_iCent%i_iPth%i", iPtZ, iCent, iX));
        h_trk_xhz_pull[iPtZ-2][iCent-1][iX] = (TH1D*) histFile->Get (Form ("h_trk_xhz_pull_iPtZ%i_iCent%i_iXhZ%i", iPtZ, iCent, iX));

        g_pth_yield_pull[0][iPtZ-2][iCent-1][iX] = (TGraph*) histFile->Get (Form ("g_pth_yield_pull_ee_iPtZ%i_iCent%i_iPth%i", iPtZ, iCent, iX));
        g_xhz_yield_pull[0][iPtZ-2][iCent-1][iX] = (TGraph*) histFile->Get (Form ("g_xhz_yield_pull_ee_iPtZ%i_iCent%i_iXhZ%i", iPtZ, iCent, iX));
        g_pth_yield_pull[1][iPtZ-2][iCent-1][iX] = (TGraph*) histFile->Get (Form ("g_pth_yield_pull_mumu_iPtZ%i_iCent%i_iPth%i", iPtZ, iCent, iX));
        g_xhz_yield_pull[1][iPtZ-2][iCent-1][iX] = (TGraph*) histFile->Get (Form ("g_xhz_yield_pull_mumu_iPtZ%i_iCent%i_iXhZ%i", iPtZ, iCent, iX));
        g_pth_yield_pull[2][iPtZ-2][iCent-1][iX] = (TGraph*) histFile->Get (Form ("g_pth_yield_pull_comb_iPtZ%i_iCent%i_iPth%i", iPtZ, iCent, iX));
        g_xhz_yield_pull[2][iPtZ-2][iCent-1][iX] = (TGraph*) histFile->Get (Form ("g_xhz_yield_pull_comb_iPtZ%i_iCent%i_iXhZ%i", iPtZ, iCent, iX));
      } // end loop over iX
    } // end loop over iCent
  } // end loop over iPtZ


  TCanvas* c_pull_pth_iPtZ2 = new TCanvas ("c_pull_pth_iPtZ2", "", 800, 600);
  c_pull_pth_iPtZ2->Divide (4, 3);
  TCanvas* c_pull_xhz_iPtZ2 = new TCanvas ("c_pull_xhz_iPtZ2", "", 800, 600);
  c_pull_xhz_iPtZ2->Divide (4, 3);
  TCanvas* c_pull_pth_iPtZ3 = new TCanvas ("c_pull_pth_iPtZ3", "", 1000, 600);
  c_pull_pth_iPtZ3->Divide (5, 3);
  TCanvas* c_pull_xhz_iPtZ3 = new TCanvas ("c_pull_xhz_iPtZ3", "", 1000, 600);
  c_pull_xhz_iPtZ3->Divide (5, 3);
  TCanvas* c_pull_pth_iPtZ4 = new TCanvas ("c_pull_pth_iPtZ4", "", 1200, 600);
  c_pull_pth_iPtZ4->Divide (6, 3);
  TCanvas* c_pull_xhz_iPtZ4 = new TCanvas ("c_pull_xhz_iPtZ4", "", 1200, 600);
  c_pull_xhz_iPtZ4->Divide (6, 3);

  float pthSigmaValues[3][3][6]; // iPtZ, iCent, iX
  float pthSigmaErrors[3][3][6]; // iPtZ, iCent, iX
  float xhzSigmaValues[3][3][6]; // iPtZ, iCent, iX
  float xhzSigmaErrors[3][3][6]; // iPtZ, iCent, iX

  for (short iX = 0; iX < 4; iX++) {
    for (short iCent : {1, 2, 3}) {
      c_pull_pth_iPtZ2->cd (4*(iCent-1)+iX+1);
      h_trk_pth_pull[0][iCent-1][iX]->Draw ("hist");

      TF1* f = new TF1 (Form ("f_trk_pth_pull_iPtZ2_iCent%i_iX%i", iCent, iX), "gaus(0)", -4, 4);
      h_trk_pth_pull[0][iCent-1][iX]->Fit (f, "RN0Q");
      f->SetLineColor (kAzure-1);
      f->SetLineWidth (1);
      f->Draw ("same");

      float chi2 = f->GetChisquare ();
      int ndf = f->GetNDF ();
      pthSigmaValues[0][iCent-1][iX] = f->GetParameter (2);
      pthSigmaErrors[0][iCent-1][iX] = f->GetParError (2);
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f #pm %.2f", f->GetParameter (2), f->GetParError (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);
      myText (0.2, 0.25, kBlack, Form ("%s < #it{p}_{T}^{ch} < %s GeV", GetMinTrkStr (iX, true).c_str (), GetMaxTrkStr (iX, true).c_str ()), 0.06);

      c_pull_xhz_iPtZ2->cd (4*(iCent-1)+iX+1);
      h_trk_xhz_pull[0][iCent-1][iX+2]->Draw ("hist");

      f = new TF1 (Form ("f_trk_xhz_pull_iPtZ2_iCent%i_iX%i", iCent, iX+2), "gaus(0)", -4, 4);
      h_trk_xhz_pull[0][iCent-1][iX+2]->Fit (f, "RN0Q");
      f->SetLineColor (kAzure-1);
      f->SetLineWidth (1);
      f->Draw ("same");

      chi2 = f->GetChisquare ();
      ndf = f->GetNDF ();
      xhzSigmaValues[0][iCent-1][iX] = f->GetParameter (2);
      xhzSigmaErrors[0][iCent-1][iX] = f->GetParError (2);
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f #pm %.2f", f->GetParameter (2), f->GetParError (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);
      myText (0.2, 0.25, kBlack, Form ("%s < #it{x}_{hZ} < %s", GetMinTrkStr (iX+2, false).c_str (), GetMaxTrkStr (iX+2, false).c_str ()), 0.06);
    } // end loop over iCent
  } // end loop over iX

  for (short iX = 0; iX < 5; iX++) {
    for (short iCent : {1, 2, 3}) {
      c_pull_pth_iPtZ3->cd (5*(iCent-1)+iX+1);
      h_trk_pth_pull[1][iCent-1][iX]->Draw ("hist");

      TF1* f = new TF1 (Form ("f_trk_pth_pull_iPtZ3_iCent%i_iX%i", iCent, iX), "gaus(0)", -4, 4);
      h_trk_pth_pull[1][iCent-1][iX]->Fit (f, "RN0Q");
      f->SetLineColor (kAzure-1);
      f->SetLineWidth (1);
      f->Draw ("same");

      float chi2 = f->GetChisquare ();
      int ndf = f->GetNDF ();
      pthSigmaValues[1][iCent-1][iX] = f->GetParameter (2);
      pthSigmaErrors[1][iCent-1][iX] = f->GetParError (2);
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f #pm %.2f", f->GetParameter (2), f->GetParError (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);
      myText (0.2, 0.25, kBlack, Form ("%s < #it{p}_{T}^{ch} < %s GeV", GetMinTrkStr (iX, true).c_str (), GetMaxTrkStr (iX, true).c_str ()), 0.06);

      c_pull_xhz_iPtZ3->cd (5*(iCent-1)+iX+1);
      h_trk_xhz_pull[1][iCent-1][iX+1]->Draw ("hist");

      f = new TF1 (Form ("f_trk_xhz_pull_iPtZ3_iCent%i_iX%i", iCent, iX+1), "gaus(0)", -4, 4);
      h_trk_xhz_pull[1][iCent-1][iX+1]->Fit (f, "RN0Q");
      f->SetLineColor (kAzure-1);
      f->SetLineWidth (1);
      f->Draw ("same");

      chi2 = f->GetChisquare ();
      ndf = f->GetNDF ();
      xhzSigmaValues[1][iCent-1][iX] = f->GetParameter (2);
      xhzSigmaErrors[1][iCent-1][iX] = f->GetParError (2);
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f #pm %.2f", f->GetParameter (2), f->GetParError (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);
      myText (0.2, 0.25, kBlack, Form ("%s < #it{x}_{hZ} < %s", GetMinTrkStr (iX+1, false).c_str (), GetMaxTrkStr (iX+1, false).c_str ()), 0.06);
    } // end loop over iCent
  } // end loop over iX

  for (short iX = 0; iX < 6; iX++) {
    for (short iCent : {1, 2, 3}) {
      c_pull_pth_iPtZ4->cd (6*(iCent-1)+iX+1);
      h_trk_pth_pull[2][iCent-1][iX]->Draw ("hist");

      TF1* f = new TF1 (Form ("f_trk_pth_pull_iPtZ4_iCent%i_iX%i", iCent, iX), "gaus(0)", -4, 4);
      h_trk_pth_pull[2][iCent-1][iX]->Fit (f, "RN0Q");
      f->SetLineColor (kAzure-1);
      f->SetLineWidth (1);
      f->Draw ("same");

      float chi2 = f->GetChisquare ();
      int ndf = f->GetNDF ();
      pthSigmaValues[2][iCent-1][iX] = f->GetParameter (2);
      pthSigmaErrors[2][iCent-1][iX] = f->GetParError  (2);
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f #pm %.2f", f->GetParameter (2), f->GetParError (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);
      myText (0.2, 0.25, kBlack, Form ("%s < #it{p}_{T}^{ch} < %s GeV", GetMinTrkStr (iX, true).c_str (), GetMaxTrkStr (iX, true).c_str ()), 0.06);

      c_pull_xhz_iPtZ4->cd (6*(iCent-1)+iX+1);
      h_trk_xhz_pull[2][iCent-1][iX]->Draw ("hist");

      f = new TF1 (Form ("f_trk_xhz_pull_iPtZ4_iCent%i_iX%i", iCent, iX), "gaus(0)", -4, 4);
      h_trk_xhz_pull[2][iCent-1][iX]->Fit (f, "RN0Q");
      f->SetLineColor (kAzure-1);
      f->SetLineWidth (1);
      f->Draw ("same");

      chi2 = f->GetChisquare ();
      ndf = f->GetNDF ();
      xhzSigmaValues[2][iCent-1][iX] = f->GetParameter (2);
      xhzSigmaErrors[2][iCent-1][iX] = f->GetParError  (2);
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f #pm %.2f", f->GetParameter (2), f->GetParError (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);
      myText (0.2, 0.25, kBlack, Form ("%s < #it{x}_{hZ} < %s", GetMinTrkStr (iX, false).c_str (), GetMaxTrkStr (iX, false).c_str ()), 0.06);
    } // end loop over iCent
  } // end loop over iX

  c_pull_pth_iPtZ2->SaveAs ("Plots/pull_pth_iPtZ2.pdf");
  c_pull_xhz_iPtZ2->SaveAs ("Plots/pull_xhz_iPtZ2.pdf");
  c_pull_pth_iPtZ3->SaveAs ("Plots/pull_pth_iPtZ3.pdf");
  c_pull_xhz_iPtZ3->SaveAs ("Plots/pull_xhz_iPtZ3.pdf");
  c_pull_pth_iPtZ4->SaveAs ("Plots/pull_pth_iPtZ4.pdf");
  c_pull_xhz_iPtZ4->SaveAs ("Plots/pull_xhz_iPtZ4.pdf");


  TCanvas* c_yield_pull = new TCanvas ("c_yield_pull", "", 600, 600);
  for (short iCent : {1, 2, 3}) {
    short iPtZ = 2;
    for (short iX = 0; iX < 4; iX++) {
      c_yield_pull->Clear ();

      TGraph* g = g_pth_yield_pull[2][iPtZ-2][iCent-1][iX];
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerSize (0.5);
      g->SetMarkerColor (kBlack);
      g->GetXaxis ()->SetTitle ("Per-event Yield");
      g->GetYaxis ()->SetTitle ("Channel pull");
      g->Draw ("AP");

      myText (0.22, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
      myText (0.22, 0.81, kBlack, Form ("%i-%i%%", GetMaxCent (iCent), GetMinCent (iCent)), 0.04);
      if (iPtZ == 4)  myText (0.22, 0.76, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", GetMinZPt (iPtZ)), 0.04);
      else            myText (0.22, 0.76, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", GetMinZPt (iPtZ), GetMaxZPt (iPtZ)), 0.04);
      myText (0.22, 0.71, kBlack, Form ("%s < #it{p}_{T}^{ch} < %s GeV", GetMinTrkStr (iX, true).c_str (), GetMaxTrkStr (iX, true).c_str ()), 0.04);

      c_yield_pull->SaveAs (Form ("Plots/2d/yield_vs_pull_iPtZ%i_iCent%i_iPth%i.pdf", iPtZ, iCent, iX));

      c_yield_pull->Clear ();

      g = g_xhz_yield_pull[2][iPtZ-2][iCent-1][iX+2];
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerSize (0.5);
      g->SetMarkerColor (kBlack);
      g->GetXaxis ()->SetTitle ("Per-event Yield");
      g->GetYaxis ()->SetTitle ("Channel pull");
      g->Draw ("AP");

      myText (0.22, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
      myText (0.22, 0.81, kBlack, Form ("%i-%i%%", GetMaxCent (iCent), GetMinCent (iCent)), 0.04);
      if (iPtZ == 4)  myText (0.22, 0.76, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", GetMinZPt (iPtZ)), 0.04);
      else            myText (0.22, 0.76, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", GetMinZPt (iPtZ), GetMaxZPt (iPtZ)), 0.04);
      myText (0.22, 0.71, kBlack, Form ("%s < #it{x}_{hZ} < %s", GetMinTrkStr (iX+2, false).c_str (), GetMaxTrkStr (iX+2, false).c_str ()), 0.04);

      c_yield_pull->SaveAs (Form ("Plots/2d/yield_vs_pull_iPtZ%i_iCent%i_iXhZ%i.pdf", iPtZ, iCent, iX+2));
    } // end loop over iX

    iPtZ = 3;
    for (short iX = 0; iX < 5; iX++) {
      c_yield_pull->Clear ();

      TGraph* g = g_pth_yield_pull[2][iPtZ-2][iCent-1][iX];
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerSize (0.5);
      g->SetMarkerColor (kBlack);
      g->GetXaxis ()->SetTitle ("Per-event Yield");
      g->GetYaxis ()->SetTitle ("Channel pull");
      g->Draw ("AP");

      myText (0.22, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
      myText (0.22, 0.81, kBlack, Form ("%i-%i%%", GetMaxCent (iCent), GetMinCent (iCent)), 0.04);
      if (iPtZ == 4)  myText (0.22, 0.76, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", GetMinZPt (iPtZ)), 0.04);
      else            myText (0.22, 0.76, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", GetMinZPt (iPtZ), GetMaxZPt (iPtZ)), 0.04);
      myText (0.22, 0.71, kBlack, Form ("%s < #it{p}_{T}^{ch} < %s GeV", GetMinTrkStr (iX, true).c_str (), GetMaxTrkStr (iX, true).c_str ()), 0.04);

      c_yield_pull->SaveAs (Form ("Plots/2d/yield_vs_pull_iPtZ%i_iCent%i_iPth%i.pdf", iPtZ, iCent, iX));

      c_yield_pull->Clear ();

      g = g_xhz_yield_pull[2][iPtZ-2][iCent-1][iX+1];
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerSize (0.5);
      g->SetMarkerColor (kBlack);
      g->GetXaxis ()->SetTitle ("Per-event Yield");
      g->GetYaxis ()->SetTitle ("Channel pull");
      g->Draw ("AP");

      myText (0.22, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
      myText (0.22, 0.81, kBlack, Form ("%i-%i%%", GetMaxCent (iCent), GetMinCent (iCent)), 0.04);
      if (iPtZ == 4)  myText (0.22, 0.76, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", GetMinZPt (iPtZ)), 0.04);
      else            myText (0.22, 0.76, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", GetMinZPt (iPtZ), GetMaxZPt (iPtZ)), 0.04);
      myText (0.22, 0.71, kBlack, Form ("%s < #it{x}_{hZ} < %s", GetMinTrkStr (iX+1, false).c_str (), GetMaxTrkStr (iX+1, false).c_str ()), 0.04);

      c_yield_pull->SaveAs (Form ("Plots/2d/yield_vs_pull_iPtZ%i_iCent%i_iXhZ%i.pdf", iPtZ, iCent, iX+1));
    } // end loop over iX

    iPtZ = 4;
    for (short iX = 0; iX < 6; iX++) {
      c_yield_pull->Clear ();

      TGraph* g = g_pth_yield_pull[2][iPtZ-2][iCent-1][iX];
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerSize (0.5);
      g->SetMarkerColor (kBlack);
      g->GetXaxis ()->SetTitle ("Per-event Yield");
      g->GetYaxis ()->SetTitle ("Channel pull");
      g->Draw ("AP");

      myText (0.22, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
      myText (0.22, 0.81, kBlack, Form ("%i-%i%%", GetMaxCent (iCent), GetMinCent (iCent)), 0.04);
      if (iPtZ == 4)  myText (0.22, 0.76, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", GetMinZPt (iPtZ)), 0.04);
      else            myText (0.22, 0.76, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", GetMinZPt (iPtZ), GetMaxZPt (iPtZ)), 0.04);
      myText (0.22, 0.71, kBlack, Form ("%s < #it{p}_{T}^{ch} < %s GeV", GetMinTrkStr (iX, true).c_str (), GetMaxTrkStr (iX, true).c_str ()), 0.04);

      c_yield_pull->SaveAs (Form ("Plots/2d/yield_vs_pull_iPtZ%i_iCent%i_iPth%i.pdf", iPtZ, iCent, iX));

      c_yield_pull->Clear ();

      g = g_xhz_yield_pull[2][iPtZ-2][iCent-1][iX];
      g->SetMarkerStyle (kOpenCircle);
      g->SetMarkerSize (0.5);
      g->SetMarkerColor (kBlack);
      g->GetXaxis ()->SetTitle ("Per-event Yield");
      g->GetYaxis ()->SetTitle ("Channel pull");
      g->Draw ("AP");

      myText (0.22, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
      myText (0.22, 0.81, kBlack, Form ("%i-%i%%", GetMaxCent (iCent), GetMinCent (iCent)), 0.04);
      if (iPtZ == 4)  myText (0.22, 0.76, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", GetMinZPt (iPtZ)), 0.04);
      else            myText (0.22, 0.76, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", GetMinZPt (iPtZ), GetMaxZPt (iPtZ)), 0.04);
      myText (0.22, 0.71, kBlack, Form ("%s < #it{x}_{hZ} < %s", GetMinTrkStr (iX, false).c_str (), GetMaxTrkStr (iX, false).c_str ()), 0.04);

      c_yield_pull->SaveAs (Form ("Plots/2d/yield_vs_pull_iPtZ%i_iCent%i_iXhZ%i.pdf", iPtZ, iCent, iX));
    } // end loop over iX
  } // end loop over iCent


  for (short iCent : {1, 2, 3}) {
    short iPtZ = 2;
    g_pth_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_pth_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_pth_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    g_xhz_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_xhz_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_xhz_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    for (int iX = 0; iX < 4; iX++) { 
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxPth (iX)+GetMinPth (iX)), pthSigmaValues[iPtZ-2][iCent-1][iX]);
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxPth (iX)-GetMinPth (iX)), pthSigmaErrors[iPtZ-2][iCent-1][iX]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_xhz_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxXhZ (iX+2)+GetMinXhZ (iX+2)), xhzSigmaValues[iPtZ-2][iCent-1][iX+2]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_xhz_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxXhZ (iX+2)-GetMinXhZ (iX+2)), xhzSigmaErrors[iPtZ-2][iCent-1][iX+2]);
    }
    iPtZ = 3;
    g_pth_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_pth_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_pth_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    g_xhz_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_xhz_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_xhz_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    for (int iX = 0; iX < 5; iX++) {
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxPth (iX)+GetMinPth (iX)), pthSigmaValues[iPtZ-2][iCent-1][iX]);
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxPth (iX)-GetMinPth (iX)), pthSigmaErrors[iPtZ-2][iCent-1][iX]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_xhz_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxXhZ (iX+1)+GetMinXhZ (iX+1)), xhzSigmaValues[iPtZ-2][iCent-1][iX+1]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_xhz_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxXhZ (iX+1)-GetMinXhZ (iX+1)), xhzSigmaErrors[iPtZ-2][iCent-1][iX+1]);
    }
    iPtZ = 4;
    g_pth_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_pth_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_pth_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    g_xhz_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_xhz_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_xhz_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    for (int iX = 0; iX < 6; iX++) {
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxPth (iX)+GetMinPth (iX)), pthSigmaValues[iPtZ-2][iCent-1][iX]);
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxPth (iX)-GetMinPth (iX)), pthSigmaErrors[iPtZ-2][iCent-1][iX]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_xhz_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxXhZ (iX)+GetMinXhZ (iX)), xhzSigmaValues[iPtZ-2][iCent-1][iX]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_xhz_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxXhZ (iX)-GetMinXhZ (iX)), xhzSigmaErrors[iPtZ-2][iCent-1][iX]);
    } // end loop over iX
  } // end loop over iCent

  TCanvas* c_pullWidth = new TCanvas ("c_pullWidth", "", 600, 600);
  const Style_t markerStyles[9] = {kOpenCircle, kOpenSquare, kOpenCross, kOpenTriangleUp, kOpenDiamond, kOpenCrossX, kOpenTriangleDown, kOpenDoubleDiamond, kOpenFourTrianglesPlus};
  const Color_t markerColors[9] = {kRed+1, kAzure-1, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1};

  c_pullWidth->Clear ();
  gPad->SetLogx ();
  bool axesDrawn = false;
  for (short iCent : {1, 2, 3}) {
    for (short iPtZ : {4, 3, 2}) {
      TGraphErrors* g = g_pth_pullWidth[iPtZ-2][iCent-1];
      g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      g->GetYaxis ()->SetTitle ("#sigma (pull)");

      g->GetYaxis ()->SetRangeUser (0.3, 1.8);

      g->SetMarkerStyle (markerStyles[iPtZ-2+3*(iCent-1)]);
      g->SetMarkerColor (markerColors[iPtZ-2+3*(iCent-1)]);
      g->SetLineColor (markerColors[iPtZ-2+3*(iCent-1)]);
      g->SetMarkerSize (1.0);

      g->Draw (axesDrawn ? "P" : "AP");
      axesDrawn = true;
    } // end loop over iPtZ
  } // end loop over iCent
  myText (0.25, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  for (short iCent : {1, 2, 3}) {
    for (short iPtZ : {4, 3, 2}) {
      if (iPtZ == 4)  myMarkerTextNoLine (0.25, 0.84-(iPtZ-2)*0.03-(iCent-1)*0.09, markerColors[iPtZ-2+3*(iCent-1)], markerStyles[iPtZ-2+3*(iCent-1)], Form ("%i-%i%%, #it{p}_{T}^{Z} > %g GeV", GetMaxCent (iCent), GetMinCent (iCent), GetMinZPt (iPtZ)), 1.00, 0.03);
      else            myMarkerTextNoLine (0.25, 0.84-(iPtZ-2)*0.03-(iCent-1)*0.09, markerColors[iPtZ-2+3*(iCent-1)], markerStyles[iPtZ-2+3*(iCent-1)], Form ("%i-%i%%, %g < #it{p}_{T}^{Z} < %g GeV", GetMaxCent (iCent), GetMinCent (iCent), GetMinZPt (iPtZ), GetMaxZPt (iPtZ)), 1.00, 0.03);
    } // end loop over iX
  } // end loop over iCent
  c_pullWidth->SaveAs ("Plots/pullWidth_pth.pdf");

  c_pullWidth->Clear ();
  gPad->SetLogx ();
  axesDrawn = false;
  for (short iCent : {1, 2, 3}) {
    for (short iPtZ : {4, 3, 2}) {
      TGraphErrors* g = g_xhz_pullWidth[iPtZ-2][iCent-1];
      g->GetXaxis ()->SetTitle ("#it{x}_{hZ}");
      g->GetYaxis ()->SetTitle ("#sigma (pull)"); 

      g->GetYaxis ()->SetRangeUser (0.3, 1.8);

      g->SetMarkerStyle (markerStyles[iPtZ-2+3*(iCent-1)]);
      g->SetMarkerColor (markerColors[iPtZ-2+3*(iCent-1)]);
      g->SetLineColor (markerColors[iPtZ-2+3*(iCent-1)]);
      g->SetMarkerSize (1.00);

      g->Draw (axesDrawn ? "P" : "AP");
      axesDrawn = true;
    } // end loop over iPtZ
  } // end loop over iCent
  myText (0.25, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  for (short iCent : {1, 2, 3}) {
    for (short iPtZ : {4, 3, 2}) {
      if (iPtZ == 4)  myMarkerTextNoLine (0.25, 0.84-(iPtZ-2)*0.03-(iCent-1)*0.09, markerColors[iPtZ-2+3*(iCent-1)], markerStyles[iPtZ-2+3*(iCent-1)], Form ("%i-%i%%, #it{p}_{T}^{Z} > %g GeV", GetMaxCent (iCent), GetMinCent (iCent), GetMinZPt (iPtZ)), 1.00, 0.03);
      else            myMarkerTextNoLine (0.25, 0.84-(iPtZ-2)*0.03-(iCent-1)*0.09, markerColors[iPtZ-2+3*(iCent-1)], markerStyles[iPtZ-2+3*(iCent-1)], Form ("%i-%i%%, %g < #it{p}_{T}^{Z} < %g GeV", GetMaxCent (iCent), GetMinCent (iCent), GetMinZPt (iPtZ), GetMaxZPt (iPtZ)), 1.00, 0.03);
    } // end loop over iPtZ
  } // end loop over iCent
  c_pullWidth->SaveAs ("Plots/pullWidth_xhz.pdf");


  // now save histograms to a rootfile
  for (short iPtZ : {2, 3, 4}) {
    for (short iCent : {1, 2, 3}) {
      g_pth_pullWidth[iPtZ-2][iCent-1]->Write ("", TObject::kOverwrite);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->Write ("", TObject::kOverwrite);
    } // end loop over iCent
  } // end loop over iPtZ
  histFile->Close ();


  // now save extracted pull values to a file
  ofstream correctionsFile;
  correctionsFile.open ("pullCorrections.dat");
  for (short iPtZ : {2, 3, 4}) {
    for (short iCent : {1, 2, 3}) {
      correctionsFile << "iPtZ = " << iPtZ << ", iCent = " << iCent << endl;
      for (short iX = 0; iX < 6; iX++) {
        correctionsFile << pthSigmaValues[iPtZ-2][iCent-1][iX] << "\t" << pthSigmaErrors[iPtZ-2][iCent-1][iX] << "\t";
      }
      correctionsFile << endl;
      for (short iX = 0; iX < 6; iX++) {
        correctionsFile << xhzSigmaValues[iPtZ-2][iCent-1][iX] << "\t" << xhzSigmaErrors[iPtZ-2][iCent-1][iX] << "\t";
      }
      correctionsFile << endl;
    }
  }
  correctionsFile.close ();

  return 0;
}

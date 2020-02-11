#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLorentzVector.h>
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


int main () {

  SetAtlasStyle();

  const int nSeeds = 2001;
  TFile* inFile;
  TTree* inTree;

  //int code = 0;
  //int id1 = 0;
  //int id2 = 0;
  //float x1pdf = 0;
  //float x2pdf = 0;
  //float Q = 0;
  //bool isValence1 = false;
  //bool isValence2 = false;

  float z_pt = 0;
  float z_eta = 0;
  float z_phi = 0;
  float z_m = 0;
  int part_n = 0;
  float part_pt[10000];
  float part_eta[10000];
  float part_phi[10000];


  float trkPtYields[2][6];
  float trkPtYieldWeightsSq[2][6];
  float trkXYields[2][6];
  float trkXYieldWeightsSq[2][6];

  float trkPtYieldAvg[2][6];
  float trkPtYieldSumSq[2][6];
  float trkPtYieldCov[6];
  float trkXYieldAvg[2][6];
  float trkXYieldSumSq[2][6];
  float trkXYieldCov[6];

  TH1D* h_trk_pth_pull[3][3][6];
  TH1D* h_trk_xhz_pull[3][3][6];

  TH1D* h_trk_pth_yield[3][3];
  TH1D* h_trk_xhz_yield[3][3];
  TH2D* h2_trk_pth_cov[3][3];
  TH2D* h2_trk_xhz_cov[3][3];

  TGraph* g_pth_yield_pull[3][3][3][6];
  TGraph* g_xhz_yield_pull[3][3][3][6];

  const float pthBins[7] = {1, 2, 4, 8, 15, 30, 60};
  const short nPthBins = sizeof (pthBins) / sizeof (pthBins[0]) - 1;
  const float xhzBins[7] = {1./60., 1./30., 1./15., 1./8., 1./4., 1./2., 1.};
  const short nXhZBins = sizeof (xhzBins) / sizeof (xhzBins[0]) - 1;

  float trk_counts[2][6] = {{}, {}};

  for (short iPtZ : {2, 3, 4}) {
    for (short iCent : {1, 2, 3}) {

      h_trk_pth_yield[iPtZ-2][iCent-1] = new TH1D (Form ("h_trk_pth_yield_iPtZ%i_iCent%i", iPtZ, iCent), ";#it{p}_{T}^{ch} [GeV]", nPthBins, pthBins);
      h_trk_pth_yield[iPtZ-2][iCent-1]->Sumw2 ();
      h_trk_xhz_yield[iPtZ-2][iCent-1] = new TH1D (Form ("h_trk_xhz_yield_iPtZ%i_iCent%i", iPtZ, iCent), ";#it{x}_{hZ}", nXhZBins, xhzBins);
      h_trk_xhz_yield[iPtZ-2][iCent-1]->Sumw2 ();
      h2_trk_pth_cov[iPtZ-2][iCent-1] = new TH2D (Form ("h2_trk_pth_cov_iPtZ%i_iCent%i", iPtZ, iCent), ";#it{p}_{T}^{ch} [GeV];#it{p}_{T}^{ch} [GeV]", nPthBins, pthBins, nPthBins, pthBins);
      h2_trk_pth_cov[iPtZ-2][iCent-1]->Sumw2 ();
      h2_trk_xhz_cov[iPtZ-2][iCent-1] = new TH2D (Form ("h2_trk_xhz_cov_iPtZ%i_iCent%i", iPtZ, iCent), ";#it{x}_{hZ};#it{x}_{hZ}", nXhZBins, xhzBins, nXhZBins, xhzBins);
      h2_trk_xhz_cov[iPtZ-2][iCent-1]->Sumw2 ();

      for (short iX = 0; iX < 6; iX++) {
        h_trk_pth_pull[iPtZ-2][iCent-1][iX] = new TH1D (Form ("h_trk_pth_pull_iPtZ%i_iCent%i_iPth%i", iPtZ, iCent, iX), ";Pull;Counts", 24, -4, 4);
        h_trk_xhz_pull[iPtZ-2][iCent-1][iX] = new TH1D (Form ("h_trk_xhz_pull_iPtZ%i_iCent%i_iXhZ%i", iPtZ, iCent, iX), ";Pull;Counts", 24, -4, 4);
        h_trk_pth_pull[iPtZ-2][iCent-1][iX]->Sumw2 ();
        h_trk_xhz_pull[iPtZ-2][iCent-1][iX]->Sumw2 ();

        g_pth_yield_pull[0][iPtZ-2][iCent-1][iX] = new TGraph ();
        g_pth_yield_pull[0][iPtZ-2][iCent-1][iX]->SetName (Form ("g_pth_yield_pull_ee_iPtZ%i_iCent%i_iPth%i", iPtZ, iCent, iX));
        g_xhz_yield_pull[0][iPtZ-2][iCent-1][iX] = new TGraph ();
        g_xhz_yield_pull[0][iPtZ-2][iCent-1][iX]->SetName (Form ("g_xhz_yield_pull_ee_iPtZ%i_iCent%i_iXhZ%i", iPtZ, iCent, iX));
        g_pth_yield_pull[1][iPtZ-2][iCent-1][iX] = new TGraph ();
        g_pth_yield_pull[1][iPtZ-2][iCent-1][iX]->SetName (Form ("g_pth_yield_pull_mumu_iPtZ%i_iCent%i_iPth%i", iPtZ, iCent, iX));
        g_xhz_yield_pull[1][iPtZ-2][iCent-1][iX] = new TGraph ();
        g_xhz_yield_pull[1][iPtZ-2][iCent-1][iX]->SetName (Form ("g_xhz_yield_pull_mumu_iPtZ%i_iCent%i_iXhZ%i", iPtZ, iCent, iX));
        g_pth_yield_pull[2][iPtZ-2][iCent-1][iX] = new TGraph ();
        g_pth_yield_pull[2][iPtZ-2][iCent-1][iX]->SetName (Form ("g_pth_yield_pull_comb_iPtZ%i_iCent%i_iPth%i", iPtZ, iCent, iX));
        g_xhz_yield_pull[2][iPtZ-2][iCent-1][iX] = new TGraph ();
        g_xhz_yield_pull[2][iPtZ-2][iCent-1][iX]->SetName (Form ("g_xhz_yield_pull_comb_iPtZ%i_iCent%i_iXhZ%i", iPtZ, iCent, iX));
      }

      const int nGroup1 = GetNumInGroup1 (iPtZ, iCent);
      const int nGroup2 = GetNumInGroup2 (iPtZ, iCent);

      for (int iSeed = 0; iSeed < nSeeds; iSeed++) {

        inFile = new TFile (Form ("rootFiles/iSeed%i_iPtZ%i_iCent%i.root", iSeed, iPtZ, iCent), "read");
        inTree = (TTree*) inFile->Get ("tree");

        //inTree->SetBranchAddress ("code",       &code);
        //inTree->SetBranchAddress ("id1",        &id1);
        //inTree->SetBranchAddress ("id2",        &id2);
        //inTree->SetBranchAddress ("x1pdf",      &x1pdf);
        //inTree->SetBranchAddress ("x2pdf",      &x2pdf);
        //inTree->SetBranchAddress ("Q",          &Q);
        //inTree->SetBranchAddress ("isValence1", &isValence1);
        //inTree->SetBranchAddress ("isValence2", &isValence2);
        inTree->SetBranchAddress ("z_pt",       &z_pt);
        inTree->SetBranchAddress ("z_eta",      &z_eta);
        inTree->SetBranchAddress ("z_phi",      &z_phi);
        inTree->SetBranchAddress ("z_m",        &z_m);
        inTree->SetBranchAddress ("part_n",     &part_n);
        inTree->SetBranchAddress ("part_pt",    &part_pt);
        inTree->SetBranchAddress ("part_eta",   &part_eta);
        inTree->SetBranchAddress ("part_phi",   &part_phi);

        const int nEvents = inTree->GetEntries ();
        assert (nEvents == nGroup1+nGroup2);

        for (int iEvent = 0; iEvent < nEvents; iEvent++) {
          inTree->GetEntry (iEvent);

          const int iGroup = (iEvent < nGroup1 ? 0 : 1); // 0 corresponds to group 1, 1 to group 2

          for (int iPart = 0; iPart < part_n; iPart++) {
            if (fabs (part_eta[iPart]) > 2.5)
              continue;
            if (DeltaPhi (part_phi[iPart], z_phi) < 3*pi/4)
              continue;

            const float trkpt = part_pt[iPart];
            const float xhz = trkpt / z_pt;

            if (trkpt < 60) {
              if (30 <= trkpt)      trk_counts[0][5] += 1.;
              else if (15 <= trkpt) trk_counts[0][4] += 1.;
              else if (8 <= trkpt)  trk_counts[0][3] += 1.;
              else if (4 <= trkpt)  trk_counts[0][2] += 1.;
              else if (2 <= trkpt)  trk_counts[0][1] += 1.;
              else if (1 <= trkpt)  trk_counts[0][0] += 1.;
            }

            if (1./60. <= xhz) {
              if (xhz <= 1./30.)      trk_counts[1][0] += 1.;
              else if (xhz <= 1./15.) trk_counts[1][1] += 1.;
              else if (xhz <= 1./8.)  trk_counts[1][2] += 1.;
              else if (xhz <= 1./4.)  trk_counts[1][3] += 1.;
              else if (xhz <= 1./2.)  trk_counts[1][4] += 1.;
              else if (xhz <= 1.)     trk_counts[1][5] += 1.;
            }
          } // end loop over iPart

          for (short iX = 0; iX < nPthBins; iX++)
            for (short iY = 0; iY < nPthBins; iY++)
              h2_trk_pth_cov[iPtZ-2][iCent-1]->SetBinContent (iX+1, iY+1, h2_trk_pth_cov[iPtZ-2][iCent-1]->GetBinContent (iX+1, iY+1) + (trk_counts[0][iX])*(trk_counts[0][iY]));
          for (short iX = 0; iX < nXhZBins; iX++)
            for (short iY = 0; iY < nXhZBins; iY++)
              h2_trk_xhz_cov[iPtZ-2][iCent-1]->SetBinContent (iX+1, iY+1, h2_trk_xhz_cov[iPtZ-2][iCent-1]->GetBinContent (iX+1, iY+1) + (trk_counts[1][iX])*(trk_counts[1][iY]));

          for (int i = 0; i < 6; i++) {
            trkPtYields[iGroup][i] += trk_counts[0][i];
            trkPtYieldWeightsSq[iGroup][i] += pow (trk_counts[0][i], 2);
            trkXYields[iGroup][i] += trk_counts[1][i];
            trkXYieldWeightsSq[iGroup][i] += pow (trk_counts[1][i], 2);
          }

          for (int i = 0; i < 6; i++) {
            trk_counts[0][i] = 0;
            trk_counts[1][i] = 0;
          }
        } // end loop over iEvent

        inFile->Close ();
        SaferDelete (inFile);

        // calculate pull on per-Z yields
        float yield1, yield2, sigma1, sigma2;
        for (short iX = 0; iX < 6; iX++) {
          yield1 = trkPtYields[0][iX] / nGroup1;
          sigma1 = sqrt (trkPtYields[0][iX]) / nGroup1;
          //sigma1 = sqrt (trkPtYieldWeightsSq[0][iX]) / nGroup1;
          yield2 = trkPtYields[1][iX] / nGroup2;
          sigma2 = sqrt (trkPtYields[1][iX]) / nGroup2;
          //sigma2 = sqrt (trkPtYieldWeightsSq[1][iX]) / nGroup2;

          h_trk_pth_yield[iPtZ-2][iCent-1]->SetBinContent (iX+1, h_trk_pth_yield[iPtZ-2][iCent-1]->GetBinContent (iX+1) + (nGroup1*yield1+nGroup2*yield2));

          trkPtYieldAvg[0][iX] += yield1;
          trkPtYieldAvg[1][iX] += yield2;
          trkPtYieldSumSq[0][iX] += yield1*yield1;
          trkPtYieldSumSq[1][iX] += yield2*yield2;
          trkPtYieldCov[iX] += yield1 * yield2;

          const float pullPth = (yield1 - yield2) / sqrt (sigma1*sigma1 + sigma2*sigma2);
          h_trk_pth_pull[iPtZ-2][iCent-1][iX]->Fill (pullPth);
          g_pth_yield_pull[0][iPtZ-2][iCent-1][iX]->SetPoint (g_pth_yield_pull[0][iPtZ-2][iCent-1][iX]->GetN (), yield1, pullPth);
          g_pth_yield_pull[1][iPtZ-2][iCent-1][iX]->SetPoint (g_pth_yield_pull[1][iPtZ-2][iCent-1][iX]->GetN (), yield2, pullPth);
          g_pth_yield_pull[2][iPtZ-2][iCent-1][iX]->SetPoint (g_pth_yield_pull[2][iPtZ-2][iCent-1][iX]->GetN (), (nGroup1*yield1+nGroup2*yield2) / (nGroup1+nGroup2), pullPth);

          yield1 = trkXYields[0][iX] / nGroup1;
          sigma1 = sqrt (trkXYields[0][iX]) / nGroup1;
          //sigma1 = sqrt (trkXYieldWeightsSq[0][iX]) / nGroup1;
          yield2 = trkXYields[1][iX] / nGroup2;
          sigma2 = sqrt (trkXYields[1][iX]) / nGroup2;
          //sigma2 = sqrt (trkXYieldWeightsSq[1][iX]) / nGroup2;

          h_trk_xhz_yield[iPtZ-2][iCent-1]->SetBinContent (iX+1, h_trk_xhz_yield[iPtZ-2][iCent-1]->GetBinContent (iX+1) + (nGroup1*yield1+nGroup2*yield2));

          trkXYieldAvg[0][iX] += yield1;
          trkXYieldAvg[1][iX] += yield2;
          trkXYieldSumSq[0][iX] += yield1*yield1;
          trkXYieldSumSq[1][iX] += yield2*yield2;
          trkXYieldCov[iX] += yield1 * yield2;

          const float pullXhZ = (yield1 - yield2) / sqrt (sigma1*sigma1 + sigma2*sigma2);
          h_trk_xhz_pull[iPtZ-2][iCent-1][iX]->Fill (pullXhZ);
          g_xhz_yield_pull[0][iPtZ-2][iCent-1][iX]->SetPoint (g_xhz_yield_pull[0][iPtZ-2][iCent-1][iX]->GetN (), yield1, pullXhZ);
          g_xhz_yield_pull[1][iPtZ-2][iCent-1][iX]->SetPoint (g_xhz_yield_pull[1][iPtZ-2][iCent-1][iX]->GetN (), yield2, pullXhZ);
          g_xhz_yield_pull[2][iPtZ-2][iCent-1][iX]->SetPoint (g_xhz_yield_pull[2][iPtZ-2][iCent-1][iX]->GetN (), (nGroup1*yield1+nGroup2*yield2) / (nGroup1+nGroup2), pullXhZ);
        }

        for (short iX = 0; iX < 6; iX++) {
          trkPtYields[0][iX]  = 0.;
          trkPtYields[1][iX]  = 0.;
          trkPtYieldWeightsSq[0][iX]  = 0.;
          trkPtYieldWeightsSq[1][iX]  = 0.;
          trkXYields[0][iX]   = 0.;
          trkXYields[1][iX]   = 0.;
          trkXYieldWeightsSq[0][iX]   = 0.;
          trkXYieldWeightsSq[1][iX]   = 0.;
        } // end loop over iX

      } // end loop over seeds


      {
        h_trk_pth_yield[iPtZ-2][iCent-1]->Scale (1./(nSeeds*(nGroup1+nGroup2)), "width");
        h_trk_xhz_yield[iPtZ-2][iCent-1]->Scale (1./(nSeeds*(nGroup1+nGroup2)), "width");

        h2_trk_pth_cov[iPtZ-2][iCent-1]->Scale (1., "width");
        h2_trk_xhz_cov[iPtZ-2][iCent-1]->Scale (1., "width");

        TH2D* h2 = h2_trk_pth_cov[iPtZ-2][iCent-1];
        TH1D* h = h_trk_pth_yield[iPtZ-2][iCent-1];
        for (short iX = 0; iX < h2->GetNbinsX (); iX++)
          for (short iY = 0; iY < h2->GetNbinsY (); iY++)
            h2->SetBinContent (iX+1, iY+1, h2->GetBinContent (iX+1, iY+1) - (nSeeds*(nGroup1+nGroup2))*(h->GetBinContent (iX+1))*(h->GetBinContent (iY+1)));
        h2 = h2_trk_xhz_cov[iPtZ-2][iCent-1];
        h = h_trk_xhz_yield[iPtZ-2][iCent-1];
        for (short iX = 0; iX < h2->GetNbinsX (); iX++)
          for (short iY = 0; iY < h2->GetNbinsY (); iY++)
            h2->SetBinContent (iX+1, iY+1, h2->GetBinContent (iX+1, iY+1) - (nSeeds*(nGroup1+nGroup2))*(h->GetBinContent (iX+1))*(h->GetBinContent (iY+1)));

        h2_trk_pth_cov[iPtZ-2][iCent-1]->Scale (1./(nSeeds*(nGroup1+nGroup2)));
        h2_trk_xhz_cov[iPtZ-2][iCent-1]->Scale (1./(nSeeds*(nGroup1+nGroup2)));
      }


      //cout << "iPtZ = " << iPtZ << ", iCent = " << iCent << endl;
      //float sigma1, sigma2;
      //for (short iX = 0; iX < 6; iX++) {
      //  trkPtYieldAvg[0][iX] = trkPtYieldAvg[0][iX] / nSeeds;
      //  trkPtYieldAvg[1][iX] = trkPtYieldAvg[1][iX] / nSeeds;
      //  sigma1 = sqrt (trkPtYieldSumSq[0][iX] / nSeeds - trkPtYieldAvg[0][iX]);
      //  sigma2 = sqrt (trkPtYieldSumSq[1][iX] / nSeeds - trkPtYieldAvg[1][iX]);
      //  trkPtYieldCov[iX] = (trkPtYieldCov[iX] / nSeeds) - trkPtYieldAvg[0][iX]*trkPtYieldAvg[1][iX];
      //  cout << "Pt correlation coeff.: " << trkPtYieldCov[iX]/(sigma1*sigma2) << endl;

      //  trkXYieldAvg[0][iX] = trkXYieldAvg[0][iX] / nSeeds;
      //  trkXYieldAvg[1][iX] = trkXYieldAvg[1][iX] / nSeeds;
      //  sigma1 = sqrt (trkXYieldSumSq[0][iX] / nSeeds - trkXYieldAvg[0][iX]);
      //  sigma2 = sqrt (trkXYieldSumSq[1][iX] / nSeeds - trkXYieldAvg[1][iX]);
      //  trkXYieldCov[iX] = (trkXYieldCov[iX] / nSeeds) - trkXYieldAvg[0][iX]*trkXYieldAvg[1][iX];

      //  //cout << "X covariance: " << trkPtYieldCov[iX] << endl;
      //} // end loop over iX

      for (short iX = 0; iX < 6; iX++) {
        trkPtYieldAvg[0][iX] = 0.;
        trkPtYieldAvg[1][iX] = 0.;
        trkPtYieldSumSq[0][iX] = 0.;
        trkPtYieldSumSq[1][iX] = 0.;
        trkPtYieldCov[iX] = 0.;
        trkXYieldAvg[0][iX] = 0.;
        trkXYieldAvg[1][iX] = 0.;
        trkXYieldSumSq[0][iX] = 0.;
        trkXYieldSumSq[1][iX] = 0.;
        trkXYieldCov[iX] = 0.;
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


  TGraphErrors* g_pth_pullWidth[3][3];
  TGraphErrors* g_xhz_pullWidth[3][3];
  for (short iCent : {1, 2, 3}) {
    short iPtZ = 2;
    g_pth_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_pth_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_pth_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    g_xhz_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_xhz_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_xhz_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    for (int iX = 0; iX < 4; iX++) { 
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxPth (iX)+GetMinPth (iX)), pthSigmaValues[iPtZ-2][iCent-1][iX]);
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxPth (iX)-GetMinPth (iX)), pthSigmaErrors[iPtZ-2][iCent-1][iX]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxXhZ (iX+2)+GetMinXhZ (iX+2)), xhzSigmaValues[iPtZ-2][iCent-1][iX+2]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_xhz_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxXhZ (iX+2)-GetMinXhZ (iX+2)), pthSigmaErrors[iPtZ-2][iCent-1][iX+2]);
    }
    iPtZ = 3;
    g_pth_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_pth_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_pth_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    g_xhz_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_xhz_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_xhz_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    for (int iX = 0; iX < 5; iX++) {
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxPth (iX)+GetMinPth (iX)), pthSigmaValues[iPtZ-2][iCent-1][iX]);
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxPth (iX)-GetMinPth (iX)), pthSigmaErrors[iPtZ-2][iCent-1][iX]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxXhZ (iX+1)+GetMinXhZ (iX+1)), xhzSigmaValues[iPtZ-2][iCent-1][iX+1]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_xhz_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxXhZ (iX+1)-GetMinXhZ (iX+1)), pthSigmaErrors[iPtZ-2][iCent-1][iX+1]);
    }
    iPtZ = 4;
    g_pth_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_pth_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_pth_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    g_xhz_pullWidth[iPtZ-2][iCent-1] = new TGraphErrors ();
    g_xhz_pullWidth[iPtZ-2][iCent-1]->SetName (Form ("g_xhz_pullWidth_iPtZ%i_iCent%i", iPtZ, iCent));
    for (int iX = 0; iX < 6; iX++) {
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxPth (iX)+GetMinPth (iX)), pthSigmaValues[iPtZ-2][iCent-1][iX]);
      g_pth_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxPth (iX)-GetMinPth (iX)), pthSigmaErrors[iPtZ-2][iCent-1][iX]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPoint (g_pth_pullWidth[iPtZ-2][iCent-1]->GetN (), 0.5*(GetMaxXhZ (iX)+GetMinXhZ (iX)), xhzSigmaValues[iPtZ-2][iCent-1][iX]);
      g_xhz_pullWidth[iPtZ-2][iCent-1]->SetPointError (g_xhz_pullWidth[iPtZ-2][iCent-1]->GetN ()-1, 0.5*(GetMaxXhZ (iX)-GetMinXhZ (iX)), pthSigmaErrors[iPtZ-2][iCent-1][iX]);
    }
  }

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
      g->GetYaxis ()->SetTitle ("#sigma");

      g->GetYaxis ()->SetRangeUser (0.3, 1.8);

      g->SetMarkerStyle (markerStyles[iPtZ-2+3*(iCent-1)]);
      g->SetMarkerColor (markerColors[iPtZ-2+3*(iCent-1)]);
      g->SetLineColor (markerColors[iPtZ-2+3*(iCent-1)]);
      g->SetMarkerSize (1.0);

      g->Draw (axesDrawn ? "P" : "AP");
      axesDrawn = true;
    }
  }
  myText (0.55, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  for (short iCent : {1, 2, 3}) {
    for (short iPtZ : {4, 3, 2}) {
      if (iPtZ == 4)  myMarkerTextNoLine (0.25, 0.84-(iPtZ-2)*0.03-(iCent-1)*0.09, markerColors[iPtZ-2+3*(iCent-1)], markerStyles[iPtZ-2+3*(iCent-1)], Form ("%i-%i%%, #it{p}_{T}^{Z} > %g GeV", GetMaxCent (iCent), GetMinCent (iCent), GetMinZPt (iPtZ)), 1.00, 0.03);
      else            myMarkerTextNoLine (0.25, 0.84-(iPtZ-2)*0.03-(iCent-1)*0.09, markerColors[iPtZ-2+3*(iCent-1)], markerStyles[iPtZ-2+3*(iCent-1)], Form ("%i-%i%%, %g < #it{p}_{T}^{Z} < %g GeV", GetMaxCent (iCent), GetMinCent (iCent), GetMinZPt (iPtZ), GetMaxZPt (iPtZ)), 1.00, 0.03);
    }
  }
  c_pullWidth->SaveAs ("Plots/pullWidth_pth.pdf");

  c_pullWidth->Clear ();
  gPad->SetLogx ();
  axesDrawn = false;
  for (short iCent : {1, 2, 3}) {
    for (short iPtZ : {4, 3, 2}) {
      TGraphErrors* g = g_pth_pullWidth[iPtZ-2][iCent-1];
      g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      g->GetYaxis ()->SetTitle ("#sigma"); 

      g->GetYaxis ()->SetRangeUser (0.3, 1.8);

      g->SetMarkerStyle (markerStyles[iPtZ-2+3*(iCent-1)]);
      g->SetMarkerColor (markerColors[iPtZ-2+3*(iCent-1)]);
      g->SetLineColor (markerColors[iPtZ-2+3*(iCent-1)]);
      g->SetMarkerSize (1.00);

      g->Draw (axesDrawn ? "P" : "AP");
      axesDrawn = true;
    }
  }
  myText (0.55, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  for (short iCent : {1, 2, 3}) {
    for (short iPtZ : {4, 3, 2}) {
      if (iPtZ == 4)  myMarkerTextNoLine (0.25, 0.84-(iPtZ-2)*0.03-(iCent-1)*0.09, markerColors[iPtZ-2+3*(iCent-1)], markerStyles[iPtZ-2+3*(iCent-1)], Form ("%i-%i%%, #it{p}_{T}^{Z} > %g GeV", GetMaxCent (iCent), GetMinCent (iCent), GetMinZPt (iPtZ)), 1.00, 0.03);
      else            myMarkerTextNoLine (0.25, 0.84-(iPtZ-2)*0.03-(iCent-1)*0.09, markerColors[iPtZ-2+3*(iCent-1)], markerStyles[iPtZ-2+3*(iCent-1)], Form ("%i-%i%%, %g < #it{p}_{T}^{Z} < %g GeV", GetMaxCent (iCent), GetMinCent (iCent), GetMinZPt (iPtZ), GetMaxZPt (iPtZ)), 1.00, 0.03);
    }
  }
  c_pullWidth->SaveAs ("Plots/pullWidth_xhz.pdf");


  // now save histograms to a rootfile
  TFile* outFile = new TFile ("outFile.root", "recreate");
  for (short iPtZ : {2, 3, 4}) {
    for (short iCent : {1, 2, 3}) {
      for (short iX = 0; iX < 6; iX++) {
        h_trk_pth_pull[iPtZ-2][iCent-1][iX]->Write ();
        h_trk_xhz_pull[iPtZ-2][iCent-1][iX]->Write ();
        for (short iChannel : {0, 1, 2}) {
          g_pth_yield_pull[iChannel][iPtZ-2][iCent-1][iX]->Write ();
          g_xhz_yield_pull[iChannel][iPtZ-2][iCent-1][iX]->Write ();
        }
      }
      h_trk_pth_yield[iPtZ-2][iCent-1]->Write ();
      h_trk_xhz_yield[iPtZ-2][iCent-1]->Write ();
      h2_trk_pth_cov[iPtZ-2][iCent-1]->Write ();
      h2_trk_xhz_cov[iPtZ-2][iCent-1]->Write ();
      g_pth_pullWidth[iPtZ-2][iCent-1]->Write ();
      g_xhz_pullWidth[iPtZ-2][iCent-1]->Write ();
    }
  }
  outFile->Close ();


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


void analyze () {
  main ();
}

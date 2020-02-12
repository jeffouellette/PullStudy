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
#include <TRandom3.h>

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

//const float bkgExponent = -5.;
//TRandom3* rndm = new TRandom3 ();
//
//float SamplePowerLaw (float x) {
//  return pow (1 + (pow (60, bkgExponent+1) - 1)*x, 1./(bkgExponent+1.));
//}

int main (int argc, char** argv) {

  if (argc < 3) {
    cout << "Usage: analyze iPtZ iCent" << endl;
    return 0;
  }

  const short iPtZ = atoi (argv[1]);
  const short iCent = atoi (argv[2]);

  SetAtlasStyle();

  const int nSeeds = 2001;

  //TFile* nchFile = new TFile ("nch_ztagged_PbPb.root", "read");
  //TH1D* h_nch_ztagged_PbPb[3];
  //h_nch_ztagged_PbPb[0] = (TH1D*) nchFile->Get ("h_nch_ztagged_PbPb_cent30_80");
  //h_nch_ztagged_PbPb[1] = (TH1D*) nchFile->Get ("h_nch_ztagged_PbPb_cent10_30");
  //h_nch_ztagged_PbPb[2] = (TH1D*) nchFile->Get ("h_nch_ztagged_PbPb_cent0_10");

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

  cout << "Working on " << GetMinZPt (iPtZ) << " - " << GetMaxZPt (iPtZ) << " GeV Z's in " << GetMaxCent (iCent) << "-" << GetMinCent (iCent) << "% central events." << endl;

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
    if (nSeeds > 100 && iSeed % (nSeeds / 100) == 0)
      cout << iSeed / (nSeeds / 100) << "\% of seeds done...\r" << flush;


    inFile = new TFile (Form ("rootFiles/iPtZ%i/iCent%i/iSeed%i_iPtZ%i_iCent%i.root", iPtZ, iCent, iSeed, iPtZ, iCent), "read");
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

      // first loop over the particles in the recorded event
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


      /*
      if (iCent != 0) {
        // now in the case of Pb+Pb sample a falling background to embed this Z in
        {
          const int nchBkg = h_nch_ztagged_PbPb[iCent-1]->GetRandom ();
          for (int iPart = 0; iPart < nchBkg; iPart++) {
            const float trkpt = SamplePowerLaw (rndm->Rndm ());
            const float xhz = trkpt / z_pt;
            if (DeltaPhi (2*pi*rndm->Rndm (), z_phi) < 3*pi/4)
              continue;

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
        }

        // then subtract the same falling background
        for (int iMixedEvent = 0; iMixedEvent < 20; iMixedEvent++)
        {
          const int nchBkg = h_nch_ztagged_PbPb[iCent-1]->GetRandom ();
          for (int iPart = 0; iPart < nchBkg; iPart++) {
            const float trkpt = SamplePowerLaw (rndm->Rndm ());
            const float xhz = trkpt / z_pt;
            if (DeltaPhi (2*pi*rndm->Rndm (), z_phi) < 3*pi/4)
              continue;

            if (trkpt < 60) {
              if (30 <= trkpt)      trk_counts[0][5] -= 1./20.;
              else if (15 <= trkpt) trk_counts[0][4] -= 1./20.;
              else if (8 <= trkpt)  trk_counts[0][3] -= 1./20.;
              else if (4 <= trkpt)  trk_counts[0][2] -= 1./20.;
              else if (2 <= trkpt)  trk_counts[0][1] -= 1./20.;
              else if (1 <= trkpt)  trk_counts[0][0] -= 1./20.;
            }

            if (1./60. <= xhz) {
              if (xhz <= 1./30.)      trk_counts[1][0] -= 1./20.;
              else if (xhz <= 1./15.) trk_counts[1][1] -= 1./20.;
              else if (xhz <= 1./8.)  trk_counts[1][2] -= 1./20.;
              else if (xhz <= 1./4.)  trk_counts[1][3] -= 1./20.;
              else if (xhz <= 1./2.)  trk_counts[1][4] -= 1./20.;
              else if (xhz <= 1.)     trk_counts[1][5] -= 1./20.;
            }
          } // end loop over iPart
        }
      }
      */

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
      //sigma1 = sqrt (trkPtYields[0][iX]) / nGroup1;
      sigma1 = sqrt (trkPtYieldWeightsSq[0][iX] / pow (nGroup1, 2) - pow (yield1, 2) / nGroup1);
      yield2 = trkPtYields[1][iX] / nGroup2;
      //sigma2 = sqrt (trkPtYields[1][iX]) / nGroup2;
      sigma2 = sqrt (trkPtYieldWeightsSq[1][iX] / pow (nGroup2, 2) - pow (yield2, 2) / nGroup2);

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
      //sigma1 = sqrt (trkXYields[0][iX]) / nGroup1;
      sigma1 = sqrt (trkXYieldWeightsSq[0][iX] / pow (nGroup1, 2) - pow (yield1, 2) / nGroup1);
      yield2 = trkXYields[1][iX] / nGroup2;
      //sigma2 = sqrt (trkXYields[1][iX]) / nGroup2;
      sigma2 = sqrt (trkXYieldWeightsSq[1][iX] / pow (nGroup2, 2) - pow (yield2, 2) / nGroup2);

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

  cout << endl << endl;


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
        

  // now save histograms to a rootfile
  TFile* outFile = new TFile (Form ("rootFiles/outFile_iPtZ%i_iCent%i.root", iPtZ, iCent), "recreate");
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
  
  
  outFile->Close ();

  return 0;
}

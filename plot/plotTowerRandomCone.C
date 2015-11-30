#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

//#include "/Users/mverweij/wrk/macros/plotUtils.C"
//#include "/Users/mverweij/wrk/macros/style.C"

//#include "ForestPFs.h"

Double_t deltaR(const Double_t phi1, const Double_t phi2, const Double_t eta1, const Double_t eta2);

void plotTowerRandomCone(TString str = "forest.root", Int_t nRCPerEvent = 1) {

  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  double minEta = -1.5;
  double maxEta = 1.5;
  double minPhi = -TMath::Pi();
  double maxPhi = TMath::Pi();
  double radiusRC = 0.4;
  
  TChain *fChain = new TChain("hiEvtAnalyzer/HiTree");
  TChain *pfTree = new TChain("rechitanalyzer/tower");
  TChain *skimTree = new TChain("skimanalysis/HltTree");
  TChain *hltTree = new TChain("hltanalysis/HltTree");
  //  TFile *f = TFile::Open(str.Data());
  fChain->Add(str.Data());
  pfTree->Add(str.Data());
  skimTree->Add(str.Data());
  hltTree->Add(str.Data());
  fChain->AddFriend(pfTree);
  fChain->AddFriend(skimTree);
  fChain->AddFriend(hltTree);
  
  if(!fChain) {
    Printf("Couldn't find pfTree. Aborting!");
    return;
  }

  Int_t nPFpart;
  Int_t pfId[32768];
  Float_t pfPt[32768];
  Float_t pfEta[32768];
  Float_t pfPhi[32768];
  Float_t pfPtVS[32768];
  
  Int_t MinBiasTriggerBit;
  Int_t phfCoincFilter;

  Int_t           hiBin;
  TBranch        *b_hiBin;   //!
  fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
  //fChain->SetBranchAddress("HLT_HIL1MinimumBiasHF1AND_v1",&MinBiasTriggerBit);
  fChain->SetBranchAddress("HLT_HIL1MinimumBiasHF1ANDExpress_v1",&MinBiasTriggerBit);
  fChain->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);

  fChain->SetBranchAddress("n", &nPFpart);
  fChain->SetBranchAddress("et", pfPt);
  fChain->SetBranchAddress("eta", pfEta);
  fChain->SetBranchAddress("phi", pfPhi);
  fChain->SetBranchAddress("vsPt", pfPtVS);

  Printf("nentries: %d",(Int_t)fChain->GetEntries());

  TList *fOutput =  new TList();
  TH1::SetDefaultSumw2();
  
  TH3D *h2CentPtRCEta = new TH3D("h2CentPtRCEta","h2CentPtRCEta;centrality;p_{T,RC};#eta",100,0,100,250,-50.,200.,60,-6,6);
  fOutput->Add(h2CentPtRCEta);
  TH3D *h2CentPtRCEtaVS = new TH3D("h2CentPtRCEtaVS","h2CentPtRCEtaVS;centrality;p_{T,RC}^{VS};#eta",100,0,100,250,-50.,200.,60,-6,6);
  fOutput->Add(h2CentPtRCEtaVS);

  TH3D *h2MultPtRCEta = new TH3D("h2MultPtRCEta","h2MultPtRCEta;centrality;p_{T,RC};#eta",3000,0,6000,150,-50.,100.,60,-6,6);
  fOutput->Add(h2MultPtRCEta);
  TH3D *h2MultPtRCEtaVS = new TH3D("h2MultPtRCEtaVS","h2MultPtRCEtaVS;centrality;p_{T,RC}^{VS};#eta",3000,0,6000,150,-50.,100.,60,-6,6);
  fOutput->Add(h2MultPtRCEtaVS);
  Printf("histos defined");

  Int_t startEntry = 0;
  Int_t lastEntry = fChain->GetEntries();//100;
  
  for (int j=startEntry; j<lastEntry; j++) {
    fChain->GetEntry(j);

    if(!MinBiasTriggerBit) continue;
    if(!phfCoincFilter) continue;

    double etaRC = gRandom->Rndm() * (maxEta - minEta) + minEta;
    double phiRC = gRandom->Rndm() * (maxPhi - minPhi) + minPhi;

    double ptRC = 0.;
    double ptRCVS = 0.;
    Int_t pfCount = 0;
    for(Int_t i = 0; i<nPFpart; i++) {
      double pt = pfPt[i];
      double ptVS = pfPtVS[i];
      double phi = pfPhi[i];
      double eta = pfEta[i];

      double dr = deltaR(phi,phiRC,eta,etaRC);
      if(dr<radiusRC) {
        ptRC+=pt;
        ptRCVS+=ptVS;
      }

      if(abs(eta)<2.) pfCount++;
      
    }
    Double_t cent = (Double_t)hiBin/2.;
    h2CentPtRCEta->Fill(cent,ptRC,etaRC);
    h2CentPtRCEtaVS->Fill(cent,ptRCVS,etaRC);

    h2MultPtRCEta->Fill(pfCount,ptRC,etaRC);
    h2MultPtRCEtaVS->Fill(pfCount,ptRCVS,etaRC);
  }

  TFile *fout = new TFile("RandomCones.root","RECREATE");
  fOutput->Write();
  fout->Write();
  fout->Close();
}

Double_t deltaR(const Double_t phi1, const Double_t phi2, const Double_t eta1, const Double_t eta2) {
  //calculate distance
  Double_t dPhi = phi1 - phi2;
  Double_t dEta = eta1 - eta2;
  dPhi = TVector2::Phi_mpi_pi(dPhi);
  return TMath::Sqrt(dPhi * dPhi + dEta * dEta);
}

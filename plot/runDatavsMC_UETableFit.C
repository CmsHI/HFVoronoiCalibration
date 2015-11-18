// macro to check the difference between Data and MC for Run1 calibrations and apply that to the Run2 MC calibrations - to get a sort of mock data UE calibrations.
// the indices (for the fit function in the root file): predictor_index = 7, i = 0, y = 0,  (for all m, n, k), m = 0, 1; n = 0; k = 0, 1, 2
// method = generate histogram from the data fit function for all the above values, fit it with the corresponding MC function = a*F(b*x); F(x) = MC function.
// where HF_Energy is the x axis and a, b are parameters
//
// now that we have the code working, lets run over all the eta bins to check if the fit closes
// take the value for x scaling which is [11] from index=7, m = 0, n = 0, k = 1 and 
// 
// Author: Raghav Kunnawalam Elayavalli
//         Rutgers, @CERN, Nov17 2015
//

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include "TCanvas.h"
#include <TMarker.h>
#include <TString.h>
#include <TVirtualFitter.h>

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

static const int bound_m = 2;
static const int bound_n = 1;
static const int bound_k = 3;
static const int bound_index = 15;
static const int array_index[bound_index] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
static const int array_m[bound_m] = {0, 1}; // ?
static const int array_n[bound_n] = {0}; // ?
static const int array_k[bound_k] = {0, 1, 2}; // track, ecal, hcal

// make a histogram from TF1 function
TH1F *functionHist(TF1 *f, TH1F* h,char *fHistname)
{
  TH1F *hF = (TH1F*)h->Clone(fHistname);
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      double var = f->Integral(h->GetBinLowEdge(i),h->GetBinLowEdge(i+1))/h->GetBinWidth(i);
      hF->SetBinContent(i,var);
      hF->SetBinError(i,0);
    }
  return hF;
}

using namespace std;


void runDatavsMC_UETableFit(int isCalo = 0){

  // get the input files
  TFile * fRun1Data = new TFile(Form("VSTraining_HFEnergy_Run1_fitfunctions_%disCalorimetric_1isData.root",isCalo),"r");
  TFile * fRun1MC   = new TFile(Form("VSTraining_HFEnergy_Run1_fitfunctions_%disCalorimetric_0isData.root",isCalo),"r");
  
  // define the functions and get them from the files.
  // these are all response for the mid eta region: predictor index = 7, and v0 - v0 which is the HF energy distribution 
  TF1 * FRun1Data_HFResponse[bound_index][bound_m][bound_n][bound_k];
  TF1 * FRun1MC_HFResponse[bound_index][bound_m][bound_n][bound_k];
  TF1 * FRun1Data_MCFit[bound_index][bound_m][bound_n][bound_k];
  
  TH1F * hData_HFResponse[bound_index][bound_m][bound_n][bound_k];
  TH1F * htest = new TH1F("htest","",1000, 0, 10000);

  for(int index = 0; index<bound_index; ++index){
    cout<<"Eta bins index = "<<index<<endl;
    for(int m = 0; m<bound_m; ++m){
      for(int n = 0; n<bound_n; ++n){
	for(int k = 0; k<bound_k; ++k){
	  FRun1Data_HFResponse[index][m][n][k] = (TF1*)fRun1Data->Get(Form("HF_energyResponseFunction_PredictorIndix%d_i0_y0_m%d_n%d_k%d", array_index[index], array_m[m], array_n[n], array_k[k]));
	  hData_HFResponse[index][m][n][k] = (TH1F*)functionHist(FRun1Data_HFResponse[index][m][n][k], htest, Form("Data_HFEnergyResponse_index%d_m%d_n%d_k%d", array_index[index], array_m[m], array_n[n], array_k[k]));
	  htest->Reset();
	  FRun1MC_HFResponse[index][m][n][k] = (TF1*)fRun1MC->Get(Form("HF_energyResponseFunction_PredictorIndix%d_i0_y0_m%d_n%d_k%d", array_index[index], array_m[m], array_n[n], array_k[k]));
	  
	  FRun1Data_MCFit[index][m][n][k] = new TF1(Form("FRun1Data_MCFit_index%d_m%d_n%d_k%d", array_index[index], array_m[m], array_n[n], array_k[k]),"[12]*([0]+[1]*([11]*x*[10])+exp(-([11]*x*[10])^2)*((-3.913998411780905*([11]*x*[10])+2.6093322745206033*([11]*x*[10])^3)*[2]+(4.931174490213579*([11]*x*[10])-6.574899320284771*([11]*x*[10])^3+1.3149798640569543*([11]*x*[10])^5)*[3]+(-5.773117374387059*([11]*x*[10])+11.546234748774118*([11]*x*[10])^3-4.618493899509647*([11]*x*[10])^5+0.43985656185806166*([11]*x*[10])^7)*[4]+(6.507479403136423*([11]*x*[10])-17.353278408363792*([11]*x*[10])^3+10.411967045018276*([11]*x*[10])^5-1.9832318180987192*([11]*x*[10])^7+0.11017954544992885*([11]*x*[10])^9)*[5]+(-7.167191940825306*([11]*x*[10])+23.89063980275102*([11]*x*[10])^3-19.112511842200817*([11]*x*[10])^5+5.460717669200234*([11]*x*[10])^7-0.6067464076889149*([11]*x*[10])^9+0.02206350573414236*([11]*x*[10])^11)*[6]+(7.771206704387521*([11]*x*[10])-31.084826817550084*([11]*x*[10])^3+31.084826817550084*([11]*x*[10])^5-11.841838787638126*([11]*x*[10])^7+1.9736397979396878*([11]*x*[10])^9-0.14353743985015913*([11]*x*[10])^11+0.0036804471756451056*([11]*x*[10])^13)*[7]+(-8.331608118589472*([11]*x*[10])+38.88083788675087*([11]*x*[10])^3-46.65700546410104*([11]*x*[10])^5+22.217621649571925*([11]*x*[10])^7-4.9372492554604275*([11]*x*[10])^9+0.5386090096865921*([11]*x*[10])^11-0.027620974855722673*([11]*x*[10])^13+0.00052611380677567*([11]*x*[10])^15)*[8]+(8.856659222944476*([11]*x*[10])-47.23551585570387*([11]*x*[10])^3+66.12972219798543*([11]*x*[10])^5-37.7884126845631*([11]*x*[10])^7+10.496781301267527*([11]*x*[10])^9-1.5268045529116403*([11]*x*[10])^11+0.11744650407012618*([11]*x*[10])^13-0.004474152536004807*([11]*x*[10])^15+0.0000657963608236001*([11]*x*[10])^17)*[9]))",0, 10000);
	  //FRun1Data_MCFit[m][n][k] = new TF1(Form("FRun1Data_MCFit_m%d_n%d_k%d", array_m[m], array_n[n], array_k[k]),Form("[0]*HF_energyResponseFunction_PredictorIndix7_i0_y0_m%d_n%d_k%d + [1]*x", array_m[m], array_n[n], array_k[k]),0, 10000);
	  for(int l = 0; l<11; l++)
	    FRun1Data_MCFit[index][m][n][k]->FixParameter(l, FRun1MC_HFResponse[index][m][n][k]->GetParameter(l));
	  FRun1Data_MCFit[index][m][n][k]->SetParameter(11,1.0);
	  FRun1Data_MCFit[index][m][n][k]->SetParameter(12,1.0);
	  hData_HFResponse[index][m][n][k]->Fit(Form("FRun1Data_MCFit_index%d_m%d_n%d_k%d", array_index[index], array_m[m], array_n[n], array_k[k]),"LL","", 0, 10000);
	}//k
      }//n
    }//m
  }// index

  TFile * fout = new TFile(Form("HFEnergy_DatavsMC_Run1_%disCalorimetric.root",isCalo),"RECREATE");
  fout->cd();

  for(int index = 0; index<bound_index; ++index){
    for(int m = 0; m<bound_m; ++m){
      for(int n = 0; n<bound_n; ++n){
	for(int k = 0; k<bound_k; ++k){
	  hData_HFResponse[index][m][n][k]->Write();
	  FRun1Data_MCFit[index][m][n][k]->Write();
	}//k
      }//n
    }//m
  }//index
}

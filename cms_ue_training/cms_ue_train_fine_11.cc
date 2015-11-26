#include <cmath>
#include <cstdlib>
#include <tr1/tuple>

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TProfile.h>

size_t nevent;
static const size_t nfourier = 1;
static const bool   selMBTrigger = true;//false;

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
	 Int_t iflag)
{
}

size_t pf_id_reduce(const Int_t pf_id)
{
  // Particle::pdgId_ PFCandidate::particleId_
  // PFCandidate::ParticleType Particle
  // 0           0  X          unknown, or dummy 
  // +211, -211  1  h          charged hadron 
  // +11, -11    2  e          electron 
  // +13, -13    3  mu         muon 
  // 22          4  gamma      photon 
  // 130         5  h0         neutral hadron 
  // 130         6  h_HF       hadronic energy in an HF tower 
  // 22          7  egamma_HF  electromagnetic energy in an HF tower

  if (pf_id == 4) {
    return 1;
  }
  else if (pf_id >= 5 && pf_id <= 7) {
    return 2;
  }

  return 0;
}

int main(int argc, char *argv[])
{
  if (argc < 2) {
    exit(EXIT_FAILURE);
  }

  bool calorimetric = true;//false;

  const char *root_tree_name = calorimetric ?
    "rechitanalyzer/tower" : "pfcandAnalyzer/pfTree";
  const char *hlt_tree_name  = "hltanalysis/HltTree";
  const char *skim_tree_name  = "skimanalysis/HltTree"; //phfConcFilter3
  static const size_t nreduced_id = 3;

  const unsigned int plot_n = 0;

  fprintf(stderr, "%s:%d:\n", __FILE__, __LINE__);

  static const size_t nedge_pseudorapidity = 15 + 1;
  double edge_pseudorapidity[nedge_pseudorapidity] = {
    -5.191, -2.650, -2.043, -1.740, -1.479, -1.131, -0.783, -0.522,
    0.522, 0.783, 1.131, 1.479, 1.740, 2.043, 2.650, 5.191
  };

  static const double cms_hcal_edge_pseudorapidity[83] = {
    -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013,
    -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853,
    -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830,
    -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
    -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609,
    -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,
    0.000,
    0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,
    0.696,  0.783,  0.879,  0.957,  1.044,  1.131,  1.218,
    1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,
    1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  2.853,
    2.964,  3.139,  3.314,  3.489,  3.664,  3.839,  4.013,
    4.191,  4.363,  4.538,  4.716,  4.889,  5.191
  };

  const size_t cms_nonhcal_nedge_pseudorapidity = 344 + 1;
  double cms_nonhcal_edge_pseudorapidity[cms_nonhcal_nedge_pseudorapidity];

  for (size_t i = 0; i < cms_nonhcal_nedge_pseudorapidity; i++) {
    cms_nonhcal_edge_pseudorapidity[i] =
      i * (2 * 2.9928 / (cms_nonhcal_nedge_pseudorapidity - 1)) -
      2.9928;
  }

  fprintf(stderr, "%s:%d:\n", __FILE__, __LINE__);

  std::vector<std::string> filename_list;

  for (int index_arg = 1; index_arg < argc; index_arg++) {
    filename_list.push_back(argv[index_arg]);
  }

  gROOT->SetStyle("Plain");
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat(0);

  TApplication application("", &argc, argv);
  TCanvas canvas("canvas", "", 1440 + 4, 810 + 28);

#if 0
  if (plot_n == 0) {
    canvas.SetLogy();
  }
#endif
  //canvas.SetLogz();

  std::vector<TH1D> root_histogram_1d;

  for (size_t plot_reduced_id = 0; plot_reduced_id < 3; plot_reduced_id++) {
    for (size_t plot_edge_source = 0; plot_edge_source < nedge_pseudorapidity - 1; plot_edge_source++) {

      for (size_t i = 0; i < 5; i++) {
	char buf[BUFSIZ];

	snprintf(buf, BUFSIZ, "root_histogram_%lu_%lu_%lu", plot_reduced_id, plot_edge_source, i);

	if (plot_reduced_id != 2) {
	  if (i < 5) {
	    root_histogram_1d.push_back(TH1D(
					     buf, "", 344, -2.9928, 2.9928));
	  }
	  else {
	  }
	}
	else {
	  if (i < 5) {
	    root_histogram_1d.push_back(TH1D(
					     buf, "", 82, -5.191, 5.191));
	    root_histogram_1d.back().GetXaxis()->Set(
						     82, cms_hcal_edge_pseudorapidity);
	  }
	  else {
	  }
	}
	root_histogram_1d.back().Sumw2();
	root_histogram_1d.back().SetMarkerStyle(20);
	root_histogram_1d.back().SetLineColor(kBlack);
      }
    }
  }

#if 0
  TProfile *profile = NULL;
#endif

  int iaccept = 0;
  for (std::vector<std::string>::const_iterator iterator_filename =
	 filename_list.begin();
       iterator_filename != filename_list.end(); iterator_filename++) {

    if(iaccept>3000) continue;    
    
    TFile f(iterator_filename->c_str());
    TTree *root_tree = (TTree *)gDirectory->Get(root_tree_name);
    TTree *hlt_tree  = (TTree *)gDirectory->Get(hlt_tree_name);
    TTree *skim_tree  = (TTree *)gDirectory->Get(skim_tree_name);
    root_tree->AddFriend(hlt_tree);
    root_tree->AddFriend(skim_tree);
    
    fprintf(stderr, "%s:%d:\n", __FILE__, __LINE__);

    Int_t nPFpart;
    Int_t pfId[32768];
    Float_t pfPt[32768];
    Float_t pfEta[32768];
    Float_t pfPhi[32768];

    Int_t MinBiasTriggerBit;
    Int_t phfCoincFilter;
    
    if (calorimetric) {
      root_tree->SetBranchAddress("n", &nPFpart);
      root_tree->SetBranchAddress("et", pfPt);
      root_tree->SetBranchAddress("eta", pfEta);
      root_tree->SetBranchAddress("phi", pfPhi);
    }
    else {
      root_tree->SetBranchAddress("nPFpart", &nPFpart);
      root_tree->SetBranchAddress("pfId", pfId);
      root_tree->SetBranchAddress("pfPt", pfPt);
      root_tree->SetBranchAddress("pfEta", pfEta);
      root_tree->SetBranchAddress("pfPhi", pfPhi);
    }

    root_tree->SetBranchAddress("HLT_HIL1MinimumBiasHF1AND_v1",&MinBiasTriggerBit);
    //    root_tree->SetBranchAddress("HLT_L1MinimumBiasHF1_OR_part1_v1",&MinBiasTriggerBit);
    root_tree->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
    root_tree->AddFriend(hlt_tree);
    root_tree->AddFriend(skim_tree);

    size_t nentries = root_tree->GetEntries();

    for (size_t i = 0; i < nentries; i++) {
      if(iaccept>3000) continue;
      root_tree->GetEntry(i);
      
      if(selMBTrigger && !MinBiasTriggerBit) continue;
      if(selMBTrigger && !phfCoincFilter) continue; 
      iaccept++;
      
      if (i % 10 == 0) {
	fprintf(stderr, "%s:%d: %ld %d\n", __FILE__, __LINE__, i, nPFpart);
      }

      for (size_t plot_reduced_id = 0; plot_reduced_id < 3; plot_reduced_id++) {
	for (size_t plot_edge_source = 0; plot_edge_source < nedge_pseudorapidity - 1; plot_edge_source++) {
	  const size_t index_histogram = (plot_reduced_id * (nedge_pseudorapidity - 1) + plot_edge_source) * 5;

	  double pt_fourier[nedge_pseudorapidity - 1][nreduced_id][nfourier][2];
	  double pt_fourier_fine[344 + 1 - 1][nreduced_id][nfourier][2];

	  memset(pt_fourier, 0,
		 (nedge_pseudorapidity - 1) * nreduced_id * nfourier * 2 *
		 sizeof(double));
	  memset(pt_fourier_fine, 0,
		 (344 + 1 - 1) * nreduced_id * nfourier * 2 *
		 sizeof(double));

	  for (Int_t j = 0; j < nPFpart; j++) {
	    size_t reduced_id = calorimetric ?
	      ((pfEta[j] >= -2.650 && pfEta[j] < 2.650) ? 2 : 2) :
	      pf_id_reduce(pfId[j]);

	    for (size_t k = 1; k < nedge_pseudorapidity; k++) {
	      if (pfEta[j] >= edge_pseudorapidity[k - 1] &&
		  pfEta[j] < edge_pseudorapidity[k]) {
		for (size_t l = 0; l < nfourier; l++) {
		  pt_fourier[k - 1][reduced_id][l][0] +=
		    pfPt[j] * cos(l * pfPhi[j]);
		  pt_fourier[k - 1][reduced_id][l][1] +=
		    pfPt[j] * sin(l * pfPhi[j]);
		}
	      }
	    }
	  }

	  const size_t nedge_pseudorapidity_fine =
	    plot_reduced_id == 2 ? 83 : cms_nonhcal_nedge_pseudorapidity;
	  const double *edge_pseudorapidity_fine =
	    plot_reduced_id == 2 ? cms_hcal_edge_pseudorapidity : cms_nonhcal_edge_pseudorapidity;

	  for (Int_t j = 0; j < nPFpart; j++) {
	    size_t reduced_id = calorimetric ?
	      ((pfEta[j] >= -2.650 && pfEta[j] < 2.650) ? 2 : 2) :
	      pf_id_reduce(pfId[j]);

	    for (size_t k = 1; k < nedge_pseudorapidity_fine; k++) {
	      if (pfEta[j] >= edge_pseudorapidity_fine[k - 1] &&
		  pfEta[j] < edge_pseudorapidity_fine[k]) {
		for (size_t l = 0; l < nfourier; l++) {
		  const double phi0 =
		    atan2(pt_fourier[plot_edge_source][reduced_id][l][1],
			  pt_fourier[plot_edge_source][reduced_id][l][0]);

		  pt_fourier_fine[k - 1][reduced_id][l][0] +=
		    pfPt[j] * cos(l * pfPhi[j] - phi0);
		  pt_fourier_fine[k - 1][reduced_id][l][1] +=
		    pfPt[j] * sin(l * pfPhi[j] - phi0);
		}
	      }
	    }
	  }

	  // Normalization

	  for (size_t j = 0; j < nedge_pseudorapidity - 1; j++) {
	    for (size_t k = 0; k < nreduced_id; k++) {
	      pt_fourier[j][k][0][0] /= 2 * M_PI;
	      for (size_t l = 1; l < nfourier; l++) {
		pt_fourier[j][k][l][0] /= M_PI;
		pt_fourier[j][k][l][1] /= M_PI;
	      }
	    }
	  }
	  for (size_t j = 0; j < nedge_pseudorapidity_fine - 1; j++) {
	    for (size_t k = 0; k < nreduced_id; k++) {
	      pt_fourier_fine[j][k][0][0] /= 2 * M_PI;
	      for (size_t l = 1; l < nfourier; l++) {
		pt_fourier_fine[j][k][l][0] /= M_PI;
		pt_fourier_fine[j][k][l][1] /= M_PI;
	      }
	    }
	  }

#if 0
	  const double bin_size = (20 * M_PI) /
	    (root_histogram_1d[0].GetNbinsX() *
	     root_histogram_1d[0].GetNbinsY());
#endif

	  double e_hf = 0;

	  for (size_t m = 0; m < 3; m++) {
	    e_hf +=
	      (pt_fourier[0                       ][m][0][0] +
	       pt_fourier[nedge_pseudorapidity - 2][m][0][0]);
	  }


	  for (size_t j = 0; j < nedge_pseudorapidity_fine - 1; j++) {
	    const double numerator_cos =
	      pt_fourier_fine[j][plot_reduced_id][plot_n][0];
	    const double numerator_sin =
	      pt_fourier_fine[j][plot_reduced_id][plot_n][1];
	    const double denominator_cos =
	      pt_fourier[plot_edge_source][plot_reduced_id][plot_n][0];
	    const double denominator_sin =
	      pt_fourier[plot_edge_source][plot_reduced_id][plot_n][1];
	    const double denominator = sqrt(denominator_cos * denominator_cos + denominator_sin * denominator_sin);

	    if (denominator != 0) {
#if 0
	      const double area_fine = 1.0 / (nedge_pseudorapidity_fine - 1);
	      const double area_coarse = 1.0 / (nedge_pseudorapidity - 1);
#else
	      const double area_fine =
		2 * M_PI * (edge_pseudorapidity_fine[j + 1] -
			    edge_pseudorapidity_fine[j]);
	      const double area_coarse =
		2 * M_PI * (edge_pseudorapidity[plot_edge_source + 1] -
			    edge_pseudorapidity[plot_edge_source]);
#endif
	      const double x = 0.5 *
		(edge_pseudorapidity_fine[j] +
		 edge_pseudorapidity_fine[j + 1]);

	      root_histogram_1d[index_histogram + 0].Fill(
							  x,
							  (numerator_cos * area_coarse) /
							  (denominator * area_fine) * e_hf);
	      root_histogram_1d[index_histogram + 1].Fill(
							  x,
							  (numerator_sin * area_coarse) /
							  (denominator * area_fine) * e_hf);
	      root_histogram_1d[index_histogram + 2].Fill(x, e_hf);
	    }
	  }
	}
      }

      nevent = i + 1;
#if 0
      if (i == 200 || i % 1000 == 0) {
	fprintf(stderr, "\rReading ... %5.1f%%", i * (100.0 / nentries));
	for (int plot_reduced_id = 2; plot_reduced_id >= 0; plot_reduced_id--) {
	  for (size_t plot_edge_source = 0; plot_edge_source < nedge_pseudorapidity - 1; plot_edge_source++) {
	    const size_t index_histogram = (plot_reduced_id * (nedge_pseudorapidity - 1) + plot_edge_source) * 5;

	    root_histogram_1d[index_histogram + 3].Divide(&root_histogram_1d[index_histogram + 0], &root_histogram_1d[index_histogram + 2], 1, 1);
	    root_histogram_1d[index_histogram + 4].Divide(&root_histogram_1d[index_histogram + 1], &root_histogram_1d[index_histogram + 2], 1, 1);
	    root_histogram_1d[index_histogram + 3].SetMarkerStyle(20);
	    root_histogram_1d[index_histogram + 4].SetMarkerStyle(24);
	    root_histogram_1d[index_histogram + 3].Draw(plot_reduced_id == 2 && plot_edge_source == 0 ? "e1x0" : "e1x0same");
	    root_histogram_1d[index_histogram + 4].Draw("e1x0same");
	  }
	}
	canvas.Update();
      }
#endif
    }
    f.Close();
  }
  fprintf(stderr, "\r                                             "
	  "                                  \n");
  fprintf(stderr, "nevent = %lu  iaccept %d\n", nevent,iaccept);

#if 1
  for (int plot_reduced_id = 0; plot_reduced_id < 3; plot_reduced_id++) {
    for (size_t plot_edge_source = 0; plot_edge_source < nedge_pseudorapidity - 1; plot_edge_source++) {
      const size_t index_histogram = (plot_reduced_id * (nedge_pseudorapidity - 1) + plot_edge_source) * 5;

      root_histogram_1d[index_histogram + 3].Divide(&root_histogram_1d[index_histogram + 0], &root_histogram_1d[index_histogram + 2], 1, 1);
      for (int i = 0; i < root_histogram_1d[index_histogram + 3].GetNbinsX(); i++) {
	fprintf(stdout, "%lu %d %lf\n", plot_edge_source, i, root_histogram_1d[index_histogram + 3].GetBinContent(i + 1));
      }
    }
  }
#endif

  return EXIT_SUCCESS;
}

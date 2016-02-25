#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include <mkl_lapack.h>

static const size_t nfourier = 5;
static const size_t norder_emulate = 9;
static const size_t norder = 9;
static const size_t nfeature = 2 * nfourier - 1;
static const size_t ncolumn = nfeature * norder + 1;

double hermite_h_normalized(const size_t n, const double x)
{
	double y;

	switch (n) {
	case 3: y = -3.913998411780905*x + 2.6093322745206033*std::pow(x,3); break;
	case 5: y = 4.931174490213579*x - 6.574899320284771*std::pow(x,3) + 1.3149798640569543*std::pow(x,5); break;
	case 7: y = -5.773117374387059*x + 11.546234748774118*std::pow(x,3) - 4.618493899509647*std::pow(x,5) + 0.43985656185806166*std::pow(x,7); break;
	case 9: y = 6.507479403136423*x - 17.353278408363792*std::pow(x,3) + 10.411967045018276*std::pow(x,5) - 1.9832318180987192*std::pow(x,7) + 0.11017954544992885*std::pow(x,9); break;
	case 11: y = -7.167191940825306*x + 23.89063980275102*std::pow(x,3) - 19.112511842200817*std::pow(x,5) + 5.460717669200234*std::pow(x,7) - 0.6067464076889149*std::pow(x,9) + 0.02206350573414236*std::pow(x,11); break;
	case 13: y = 7.771206704387521*x - 31.084826817550084*std::pow(x,3) + 31.084826817550084*std::pow(x,5) - 11.841838787638126*std::pow(x,7) + 1.9736397979396878*std::pow(x,9) - 0.14353743985015913*std::pow(x,11) + 0.0036804471756451056*std::pow(x,13); break;
	case 15: y = -8.331608118589472*x + 38.88083788675087*std::pow(x,3) - 46.65700546410104*std::pow(x,5) + 22.217621649571925*std::pow(x,7) - 4.9372492554604275*std::pow(x,9) + 0.5386090096865921*std::pow(x,11) - 0.027620974855722673*std::pow(x,13) + 0.00052611380677567*std::pow(x,15); break;
	case 17: y = 8.856659222944476*x - 47.23551585570387*std::pow(x,3) + 66.12972219798543*std::pow(x,5) - 37.7884126845631*std::pow(x,7) + 10.496781301267527*std::pow(x,9) - 1.5268045529116403*std::pow(x,11) + 0.11744650407012618*std::pow(x,13) - 0.004474152536004807*std::pow(x,15) + 0.0000657963608236001*std::pow(x,17); break;
	default: fprintf(stderr, "%s:%d: n = %lu\n", __FILE__, __LINE__, n); exit(1); break;
	}

	return y;
}

int main(int argc, char *argv[])
{
	setenv("KMP_DUPLICATE_LIB_OK", "TRUE", 1);

	if (argc < 2) {
		exit(EXIT_FAILURE);
	}

	bool calorimetric = false;

	for (int i = 1; i < argc; i++) {
		if (strncmp(argv[i], "-c", 2) == 0 ||
			strncmp(argv[i], "-c", 2) == 0) {
			calorimetric = true;
			break;
		}
	}

	const char *root_tree_name = calorimetric ?
		"rechitanalyzer/tower" : "pfcandAnalyzer/pfTree";
	const char *hlt_tree_name = "hltanalysis/HltTree";
	const char *skim_tree_name = "skimanalysis/HltTree"; //phfConcFilter3
	static const size_t nreduced_id = 3;
	size_t nevent = 0;

	Int_t mb_trigger_bit;
	Int_t hf_coincidence_bit;

	for (int index_file = 0; index_file < argc - 1; index_file++) {
		if (argv[index_file + 1][0] == '-') {
			continue;
		}

		fprintf(stderr, "%s:%d: Count %s\n", __FILE__, __LINE__, argv[index_file + 1]);

		TFile *f = TFile::Open(argv[index_file + 1]);
		TTree *root_tree = reinterpret_cast<TTree *>(gDirectory->Get(root_tree_name));
		TTree *hlt_tree = reinterpret_cast<TTree *>(gDirectory->Get(hlt_tree_name));
		TTree *skim_tree = reinterpret_cast<TTree *>(gDirectory->Get(skim_tree_name));
		size_t nevent_file = root_tree->GetEntries();

		nevent_file = std::min(static_cast<size_t>(20000), nevent_file);

		bool is_mc = false;

		if (hlt_tree->GetBranch("HLT_HIL1MinimumBiasHF1AND_v1")) {
			hlt_tree->SetBranchAddress(
				"HLT_HIL1MinimumBiasHF1AND_v1", &mb_trigger_bit);
		}
		else {
			fprintf(stderr, "%s:%d: warning: MB trigger bit not found (ignore if this is MC)\n", __FILE__, __LINE__);
			is_mc = true;
		}
		if (skim_tree->GetBranch("phfCoincFilter3")) {
			skim_tree->SetBranchAddress(
				"phfCoincFilter3", &hf_coincidence_bit);
		}
		else {
			fprintf(stderr, "%s:%d: warning: HF coincidence skim bit not found (ignore if this is MC)\n", __FILE__, __LINE__);
			is_mc = true;
		}

		for (size_t i = 0; i < nevent_file; i++) {
			hlt_tree->GetEntry(i);
			skim_tree->GetEntry(i);

			if (!(is_mc || (mb_trigger_bit && hf_coincidence_bit))) {
				continue;
			}

			nevent++;

			if (i % 100 == 0) {
				fprintf(stderr, "%s:%d: i = %lu, nevent = %lu\n", __FILE__, __LINE__, i, nevent);
			}
		}

		f->Close();
	}

	fprintf(stderr, "%s:%d: nevent = %lu\n", __FILE__, __LINE__, nevent);

	//nevent = std::min(30000UL, nevent);

	if (nevent < 1) {
		// No point to continue, and LAPACK would also fail in this
		// case.
		exit(EXIT_FAILURE);
	}

	size_t nvariable = 0;

#if 1
	static const size_t nedge = 15 + 1;
	double edge_pseudorapidity[nedge] = {
		-5.191, -2.650, -2.043, -1.740, -1.479, -1.131, -0.783, -0.522,
		0.522, 0.783, 1.131, 1.479, 1.740, 2.043, 2.650, 5.191
	};
#else
	static const size_t nedge = 82 + 1;
	static const double edge_pseudorapidity[nedge] = {
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
#endif

	const size_t nregularization = nfeature * (norder - 1);
	const size_t feature_size = nevent + nregularization;

	static const double default_value = 0;

	// Use Fortran ordering (row = event number) for LAPACK
	std::vector<double> feature(ncolumn * feature_size, default_value);
	std::vector<double> feature_copy(ncolumn * feature_size, default_value);

	const size_t target_size =
		std::max(std::max(feature_size, ncolumn),
				 static_cast<size_t>(1));
	const size_t nrhs = (nedge - 1) * nreduced_id * nfeature;

	std::vector<double> target(nrhs * target_size, default_value);
	std::vector<double> target_copy(nrhs * target_size, default_value);

	for (size_t k = 1; k < nedge; k++) {
		const size_t l_min = k == 1 || k == nedge - 1 ? 2 : 0;

		for (size_t l = l_min; l < 3; l++) {
			nvariable += 2 * nfourier - 1;
		}
	}

	size_t nevent_read = 0;

	std::vector<double> density_table;
	std::vector<std::vector<double> > scale_check_table(nfourier, std::vector<double>());
	std::vector<double> weight_table;

	// This is valid for nfourier < 18 (where interference behavior
	// with the HF geometry starts to appear)

	std::vector<double> scale(nfourier, 1.0 / 200.0);

	if (nfourier >= 1) {
		scale[0] = 1.0 / 5400.0;
	}
	if (nfourier >= 2) {
		scale[1] = 1.0 / 130.0;
	}
	if (nfourier >= 3) {
		scale[2] = 1.0 / 220.0;
	}

	for (int index_file = 0; index_file < argc - 1; index_file++) {
		if (argv[index_file + 1][0] == '-') {
			continue;
		}
		TFile *f = TFile::Open(argv[index_file + 1]);
		TTree *root_tree = static_cast<TTree *>(gDirectory->Get(root_tree_name));
		TTree *hlt_tree = reinterpret_cast<TTree *>(gDirectory->Get(hlt_tree_name));
		TTree *skim_tree = reinterpret_cast<TTree *>(gDirectory->Get(skim_tree_name));
		size_t nevent_file = root_tree->GetEntries();

		nevent_file = std::min(static_cast<size_t>(20000), nevent_file);

		bool is_mc = false;

		if (hlt_tree->GetBranch("HLT_HIL1MinimumBiasHF1AND_v1")) {
			hlt_tree->SetBranchAddress(
				"HLT_HIL1MinimumBiasHF1AND_v1", &mb_trigger_bit);
		}
		else {
			fprintf(stderr, "%s:%d: warning: MB trigger bit not found (ignore if this is MC)\n", __FILE__, __LINE__);
			is_mc = true;
		}
		if (skim_tree->GetBranch("phfCoincFilter3")) {
			skim_tree->SetBranchAddress(
				"phfCoincFilter3", &hf_coincidence_bit);
		}
		else {
			fprintf(stderr, "%s:%d: warning: HF coincidence skim bit not found (ignore if this is MC)\n", __FILE__, __LINE__);
			is_mc = true;
		}

		Int_t nPFpart;
		Int_t pfId[32768];
		Float_t pfPt[32768];
		Float_t pfEta[32768];
		Float_t pfPhi[32768];

                std::vector<int>           *pfIdVec  = 0;
                std::vector<float>         *pfPtVec  = 0;
                std::vector<float>         *pfEtaVec = 0;
                std::vector<float>         *pfPhiVec = 0;

		if (calorimetric) {
			root_tree->SetBranchAddress("n", &nPFpart);
			root_tree->SetBranchAddress("et", pfPt);
			root_tree->SetBranchAddress("eta", pfEta);
			root_tree->SetBranchAddress("phi", pfPhi);
		}
		else {
                  root_tree->SetBranchAddress("nPFpart", &nPFpart);
                  root_tree->SetBranchAddress("pfId", &pfIdVec);
                  root_tree->SetBranchAddress("pfPt", &pfPtVec);
                  root_tree->SetBranchAddress("pfEta", &pfEtaVec);
                  root_tree->SetBranchAddress("pfPhi", &pfPhiVec);
		}

		fprintf(stderr, "%s:%d: %s is_mc = %d\n", __FILE__, __LINE__, argv[index_file + 1], static_cast<int>(is_mc));

	for (size_t i = 0; i < nevent_file; i++) {
		root_tree->GetEntry(i);
		hlt_tree->GetEntry(i);
		skim_tree->GetEntry(i);

		if (!(is_mc || (mb_trigger_bit && hf_coincidence_bit))) {
			continue;
		}

		double pt_fourier[nedge - 1][nreduced_id][nfourier][2];

		for (size_t k = 1; k < nedge; k++) {
			for (size_t l = 0; l < nreduced_id; l++) {
				for (size_t m = 0; m < nfourier; m++) {
					for (size_t re_or_im = 0; re_or_im < 2;
						 re_or_im++) {
						pt_fourier[k - 1][l][m][re_or_im] = 0;
					}
				}
			}
		}
		for (Int_t j = 0; j < nPFpart; j++) {
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

                  if (!calorimetric) {
                    pfId[j] = pfIdVec->at(j);
                    pfPt[j] = pfPtVec->at(j);
                    pfEta[j] = pfEtaVec->at(j);
                    pfPhi[j] = pfPhiVec->at(j);
                  }

			size_t reduced_id = 0;

			if (calorimetric) {
				reduced_id = (pfEta[j] >= -2.650 && pfEta[j] < 2.650) ? 2 : 2;
			}
			else {
				if (pfId[j] == 4) {
					reduced_id = 1;
				}
				else if (pfId[j] >= 5 && pfId[j] <= 7) {
					reduced_id = 2;
				}
			}

			for (size_t k = 1; k < nedge; k++) {
				if (pfEta[j] >= edge_pseudorapidity[k - 1] &&
					pfEta[j] < edge_pseudorapidity[k]) {
					for (size_t l = 0; l < nfourier; l++) {
						pt_fourier[k - 1][reduced_id][l][0] +=
							pfPt[j] * cos(l * pfPhi[j]);
						pt_fourier[k - 1][reduced_id][l][1] +=
							pfPt[j] * sin(l * pfPhi[j]);
					}
					break;
				}
			}
		}

		double eta_1_sum = 0;

		for (Int_t j = 0; j < nPFpart; j++) {
			if (pfEta[j] >= -1 && pfEta[j] < 1) {
				eta_1_sum += pfPt[j];
			}
		}

		const size_t row = nevent_read;

		if (row >= nevent) {
			break;
		}

		feature[(0) * feature_size + row] = 1;
		feature[(1) * feature_size + row] = 0;
		for (size_t m = 0; m < 3; m++) {
		feature[(1) * feature_size + row] += scale[0] *
			(pt_fourier[0        ][m][0][0] +
			 pt_fourier[nedge - 2][m][0][0]);
		}

		density_table.push_back(
			feature[(1) * feature_size + row] / scale[0]);
		scale_check_table[0].push_back(
			feature[(1) * feature_size + row] / scale[0]);

		// static const double weight = 1;

		if (i % 1000 == 0) {
			fprintf(stderr, "%s:%d: %lu %f\n", __FILE__, __LINE__, i, feature[feature_size + row]);
		}
#if 0
		for (size_t l = 1; l < norder; l++) {
			feature[(l + 1) * feature_size + row] =
				feature[(l) * feature_size + row] *
				feature[(1) * feature_size + row];
		}
#else
		// Hermite
		for (size_t l = 2; l < norder + 1; l++) {
			feature[(l) * feature_size + row] =
				hermite_h_normalized(2 * l - 1,
					 feature[(1) * feature_size + row]) *
				exp(-feature[(1) * feature_size + row] *
					 feature[(1) * feature_size + row]);
		}
#endif
		for (size_t k = 1; k < nfourier; k++) {
			feature[((2 * k - 1) * norder + 1) * feature_size + row] = 0;
			for (size_t m = 0; m < 3; m++) {
			feature[((2 * k - 1) * norder + 1) * feature_size + row] += scale[k] *
				(pt_fourier[0        ][m][k][0] +
				 pt_fourier[nedge - 2][m][k][0]);
			}
#if 0
			for (size_t l = 1; l < norder; l++) {
				feature[((2 * k - 1) * norder + l + 1) * feature_size + row] =
					feature[((2 * k - 1) * norder + l) * feature_size + row] *
					feature[((2 * k - 1) * norder + 1) * feature_size + row];
			}
#else
			for (size_t l = 2; l < norder + 1; l++) {
				feature[((2 * k - 1) * norder + l) * feature_size + row] =
					hermite_h_normalized(2 * l - 1,
						 feature[((2 * k - 1) * norder + 1) * feature_size + row]) *
					exp(-feature[((2 * k - 1) * norder + 1) * feature_size + row] *
						 feature[((2 * k - 1) * norder + 1) * feature_size + row]);
			}
#endif
			feature[((2 * k    ) * norder + 1) * feature_size + row] = 0;
			for (size_t m = 0; m < 3; m++) {
			feature[((2 * k    ) * norder + 1) * feature_size + row] += scale[k] *
				(pt_fourier[0        ][m][k][1] +
				 pt_fourier[nedge - 2][m][k][1]);
			}
#if 0
			for (size_t l = 1; l < norder; l++) {
				feature[((2 * k    ) * norder + l + 1) * feature_size + row] =
					feature[((2 * k    ) * norder + l) * feature_size + row] *
					feature[((2 * k    ) * norder + 1) * feature_size + row];
			}
#else
			for (size_t l = 2; l < norder + 1; l++) {
				feature[((2 * k    ) * norder + l) * feature_size + row] =
					hermite_h_normalized(2 * l - 1,
						 feature[((2 * k    ) * norder + 1) * feature_size + row]) *
					exp(-feature[((2 * k    ) * norder + 1) * feature_size + row] *
						 feature[((2 * k    ) * norder + 1) * feature_size + row]);
			}
#endif
			scale_check_table[k].push_back(std::max(feature[((2 * k - 1) * norder + 1) * feature_size + row], feature[((2 * k    ) * norder + 1) * feature_size + row]) / scale[k]);
		}

		for (size_t j = 1; j < nedge; j++) {
			for (size_t k = 0; k < nreduced_id; k++) {
				const size_t column_0 = ((j - 1) * nreduced_id + k) * nfeature;

				target[column_0 * target_size + row] = pt_fourier[j - 1][k][0][0];

				for (size_t l = 1; l < nfourier; l++) {
					for (size_t re_or_im = 0; re_or_im < 2;
						 re_or_im++) {
						const size_t column = column_0 + (2 * l - 1 + re_or_im);

						target[column * target_size + row] = pt_fourier[j - 1][k][l][re_or_im];
					}
				}
			}
		}

		for (size_t l = 0; l < nfeature; l++) {
			if (!std::isfinite(feature[(l) * feature_size + row])) {
				for (size_t k = 0; k < nfeature; k++) {
					fprintf(stderr, "%g ", feature[(k) * feature_size + row]);
				}
				fprintf(stderr, "\n");
				break;
			}
		}

#if 0
		if (i % 1000 == 0) {
			fprintf(stderr, "\rReading ... %5.1f%%", i * (100.0 / nevent));
		}
#endif
		nevent_read++;
	}
		f->Close();
	}
#if 0
	fprintf(stderr, "\r                                             "
			"                                  \n");
#endif
	fprintf(stderr, "%s:%d: nevent = %lu\n", __FILE__, __LINE__, nevent);

	// Statistics scaling

	for (size_t i = 0; i < nfourier; i++) {
		std::sort(scale_check_table[i].begin(), scale_check_table[i].end());

		const size_t index_95 = static_cast<size_t>(rint((scale_check_table[i].size() - 1) * 0.95));
		const size_t index_99 = static_cast<size_t>(rint((scale_check_table[i].size() - 1) * 0.99));

		fprintf(stderr, "%s:%d: %lu %f %f\n", __FILE__, __LINE__, i, scale_check_table[i][index_95], scale_check_table[i][index_99]);
	}

	std::sort(density_table.begin(), density_table.end());

	for (size_t i = 0; i < nevent; i++) {
		const double e_hf = feature[(1) * feature_size + i] / scale[0];
		const std::vector<double>::const_iterator density_table_iterator = std::lower_bound(density_table.begin(), density_table.end(), e_hf);
		double weight = 0;

		if (density_table_iterator != density_table.end()) {
			static const int naverage = 2 * 16 + 1;
			const int index = density_table_iterator - density_table.begin();
			const int index_min = std::max(0, index - (naverage - 1) / 2);
			const int index_max = std::min(static_cast<int>(density_table.size() - 1), index + (naverage - 1) / 2);
			const double range_min = density_table[index_min];
			const double range_max = density_table[index_max];

			const double inverse_density = (range_max - range_min) * nevent / (index_max - index_min);

			weight = sqrt(std::max(0.0, inverse_density));

			if (i % 1000 == 0) {
				fprintf(stderr, "%s:%d: %lu %f %f\n", __FILE__, __LINE__, i, e_hf, weight);
			}
		}

#if 0
		// This is probably wrong
		double max_abs_feature = 0;

		for (size_t j = 1; j < ncolumn; j += norder) {
			max_abs_feature = std::max(max_abs_feature, fabs(feature[j * feature_size + i]));
		}
		weight = 1.0 / (1.0 + exp(12.0 * (max_abs_feature - 0.8)));
#endif

		for (size_t j = 0; j < ncolumn; j++) {
			feature_copy[j * feature_size + i] = feature[j * feature_size + i];
			feature[j * feature_size + i] *= weight;
		}
		for (size_t j = 0; j < nrhs; j++) {
			target_copy[j * target_size + i] = target[j * target_size + i];
			target[j * target_size + i] *= weight;
		}

		weight_table.push_back(weight);
	}

	std::sort(weight_table.begin(), weight_table.end());

	// Regularization

	//const double tau = 1e-3 * nevent;
	size_t row = nevent;

	for (size_t i = 0; i < nfeature; i++) {
		const size_t index_95 = static_cast<size_t>(rint((weight_table.size() - 1) * 0.95));

		const double tau = (i == 0 ? 1e-6 : 1e-4) * nevent * weight_table[index_95];

		for (size_t j = 1; j < norder; j++) {
			const size_t column = i * norder + j + 1;

			feature[column * feature_size + row] = tau;
			//target[column * feature_size + row] = tau;

			for (size_t k = 1; k < nrhs; k++) {
				target[k * target_size + row] = 0;
			}

			row++;
		}
	}

#if 0
	fprintf(stdout, "feature =\n");
	for (size_t j = 0; j < ncolumn; j++) {
		for (size_t i = 0; i < feature_size; i++) {
			fprintf(stdout, "%g ", feature[j * feature_size + i]);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "target =\n");
	for (size_t j = 0; j < nrhs; j++) {
		for (size_t i = 0; i < target_size; i++) {
			fprintf(stdout, "%g ", target[j * target_size + i]);
		}
		fprintf(stdout, "\n");
	}
#endif

	std::vector<double> singular_value(
		std::min(feature_size, ncolumn),
		0);
	std::vector<double> working_array(
		3 * std::min(feature_size, ncolumn) +
		std::max(std::max(2 * std::min(feature_size, ncolumn),
						  std::max(feature_size, ncolumn)),
				 nrhs),
		0);
	double rank_condition = 0;

	int feature_size_int = feature_size;
	int ncolumn_int = ncolumn;
	int nrhs_int = nrhs;
	int target_size_int = target_size;
	int working_array_size = working_array.size();
	int rank;
	int info;

	dgelss_(// m, n, nrhs
			&feature_size_int, &ncolumn_int, &nrhs_int,
			// a, lda
			&feature[0], &feature_size_int,
			// b, ldb
			&target[0], &target_size_int,
			// s, rcond, rank
			&singular_value[0], &rank_condition, &rank,
			// work, lwork, info
			&working_array[0], &working_array_size, &info);

	fprintf(stderr, "%s:%d: DGELSS info = %d, rank = %d\n", __FILE__, __LINE__, info, rank);

	if (info < 0) {
		// info = -5 if LDA.LT.MAX( 1, M )
		return EXIT_FAILURE;
	}

	// Format:
	// pf_id pseudorapidity_index fourier_order re_or_im feature_index
	//
	// pseudorapidity_index is the bin defined by edge_pseudorapidity
	// re_or_im == 0 is cos(fourier_order * phi)
	// re_or_im == 1 is sin(fourier_order * phi)
	//
	// feature_index is packed HF Fourier decomposition wrt. the
	// following (unnormalized) basis functions: constant, cos(phi),
	// sin(phi), cos(2 phi), sin(2 phi) ...

	for (size_t k = 0; k < nreduced_id; k++) {
		for (size_t j = 1; j < nedge; j++) {
			const size_t column_0 = ((j - 1) * nreduced_id + k) * nfeature;
			size_t i_count = 0;

			for (size_t i = 0; i < ncolumn; i++) {
				fprintf(stdout, "%lu %lu %d %d %lu %.8e\n",
						k, j - 1, 0, 0, i_count, target[column_0 * target_size_int + i]);
				i_count++;
				// if i = norder + 1, 2 * norder + 1, ..., fill in the
				// remaining zeros
				if (i >= norder && (i - 1) % norder == norder - 1) {
					for (size_t m = norder; m < norder_emulate; m++) {
						fprintf(stdout, "%lu %lu %d %d %lu 0\n",
								k, j - 1, 0, 0, i_count);
						i_count++;
					}
				}
			}

			for (size_t l = 1; l < nfourier; l++) {
				for (size_t re_or_im = 0; re_or_im < 2;
					 re_or_im++) {
					size_t i_count = 0;

					for (size_t i = 0; i < ncolumn; i++) {
						const size_t column = column_0 + (2 * l - 1 + re_or_im);

						fprintf(stdout, "%lu %lu %lu %lu %lu %.8e\n",
								k, j - 1, l, re_or_im, i_count, target[column * target_size_int + i]);
						i_count++;
						if (i >= norder && (i - 1) % norder == norder - 1) {
							for (size_t m = norder; m < norder_emulate; m++) {
								fprintf(stdout, "%lu %lu %lu %lu %lu 0\n",
										k, j - 1, l, re_or_im, i_count);
								i_count++;
							}
						}
					}
				}
			}
		}
	}


	if (false) {
		for (size_t j = 1; j < 16; j++) {
			for (size_t k = 0; k < 3; k++) {
				const size_t column_0 = ((j - 1) * nreduced_id + k) * nfeature;

				// for (size_t l = 1; l < nfourier; l++) {
				// 	for (size_t re_or_im = 0; re_or_im < 2;
				// 		 re_or_im++) {
				// const size_t column = column_0 + (2 * l - 1 + re_or_im);

				for (size_t i = 0; i < nevent; i++) {
					double s = 0;

					for (size_t l = 0; l < ncolumn; l++) {
						s += target[column_0 * target_size_int + l] * feature_copy[l * feature_size + i];
					}
					fprintf(stderr, "%s:%d: %lu %lu %lu %.8e %.8e %.8e\n", __FILE__, __LINE__, i, j, k, feature_copy[nevent + i] / scale[0], target_copy[column_0 * target_size + i], s);

					s = 0;

					for (size_t l = 0; l < ncolumn; l++) {
						s += target[column_0 * target_size_int + l] * feature_copy[l * feature_size + i] * 2;
					}
					fprintf(stderr, "%s:%d: %lu %lu %lu %.8e %.8e %.8e\n", __FILE__, __LINE__, i, j, k, feature_copy[nevent + i] / scale[0] * 2, target_copy[column_0 * target_size + i] * 2, s);
				}
			}
		}
	}

	return EXIT_SUCCESS;
}

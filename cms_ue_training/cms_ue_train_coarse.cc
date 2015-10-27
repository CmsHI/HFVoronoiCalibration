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

int main(int argc, char *argv[])
{
	setenv("KMP_DUPLICATE_LIB_OK", "TRUE", 1);

	if (argc < 2) {
		exit(EXIT_FAILURE);
	}

	bool calorimetric = true;

	const char *root_tree_name = calorimetric ?
		"rechitanalyzer/tower" : "pfcandAnalyzer/pfTree";
	static const size_t nreduced_id = 3;
	size_t nevent = 0;

	for (int index_file = 0; index_file < argc - 1; index_file++) {
		TFile f(argv[index_file + 1]);
		TTree *root_tree = reinterpret_cast<TTree *>(gDirectory->Get(root_tree_name));
		size_t nevent_file = root_tree->GetEntries();

		// nevent_file = std::max(static_cast<size_t>(1000), nevent_file);

		nevent += nevent_file;
		f.Close();
	}

	if (nevent < 1) {
		// No point to continue, and LAPACK would also fail in this
		// case.
		exit(EXIT_FAILURE);
	}

	size_t nvariable = 0;

#if 1
	static const size_t nedge = 8;
	static const double edge_pseudorapidity[nedge] = {
		-5.191, -2.650, -1.479, -0.522, 0.522, 1.479, 2.650, 5.191
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
	// Use Fortran ordering (row = event number) for LAPACK
	std::vector<double> feature(ncolumn * nevent, NAN);
	std::vector<double> feature_copy(ncolumn * nevent, NAN);

	size_t target_size =
		std::max(std::max(nevent, ncolumn), static_cast<size_t>(1));

	const size_t nrhs = (nedge - 1) * nreduced_id * nfeature;

	std::vector<double> target(nrhs * target_size, NAN);
	std::vector<double> target_copy(nrhs * target_size, NAN);

	for (size_t k = 1; k < nedge; k++) {
		const size_t l_min = k == 1 || k == nedge - 1 ? 2 : 0;

		for (size_t l = l_min; l < 3; l++) {
			nvariable += 2 * nfourier - 1;
		}
	}

	size_t nevent_read = 0;

	std::vector<double> density_table;
	std::vector<std::vector<double> > scale_check_table(nfourier, std::vector<double>());

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
		TFile f(argv[index_file + 1]);
		TTree *root_tree = static_cast<TTree *>(gDirectory->Get(root_tree_name));
		size_t nevent_file = root_tree->GetEntries();

		// nevent_file = std::max(static_cast<size_t>(1000), nevent_file);

		Int_t nPFpart;
		Int_t pfId[32768];
		Float_t pfPt[32768];
		Float_t pfEta[32768];
		Float_t pfPhi[32768];

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

	fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, argv[index_file + 1]);

#if 1
	for (size_t i = 0; i < nevent_file; i++) {
		root_tree->GetEntry(i);

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

		feature[(0) * nevent + (nevent_read + i)] = 1;
		feature[(1) * nevent + (nevent_read + i)] = 0;
		for (size_t m = 0; m < 3; m++) {
		feature[(1) * nevent + (nevent_read + i)] += scale[0] *
			(pt_fourier[0        ][m][0][0] +
			 pt_fourier[nedge - 2][m][0][0]);
		}

		density_table.push_back(
			feature[(1) * nevent + (nevent_read + i)] / scale[0]);
		scale_check_table[0].push_back(
			feature[(1) * nevent + (nevent_read + i)] / scale[0]);

		// static const double weight = 1;

		if (i % 1000 == 0) {
			fprintf(stderr, "%s:%d: %lu %f\n", __FILE__, __LINE__, i, feature[nevent + i]);
		}
		for (size_t l = 1; l < norder; l++) {
			feature[(l + 1) * nevent + (nevent_read + i)] =
				feature[(l) * nevent + (nevent_read + i)] *
				feature[(1) * nevent + (nevent_read + i)];
		}
		for (size_t k = 1; k < nfourier; k++) {
			feature[((2 * k - 1) * norder + 1) * nevent + (nevent_read + i)] = 0;
			for (size_t m = 0; m < 3; m++) {
			feature[((2 * k - 1) * norder + 1) * nevent + (nevent_read + i)] += scale[k] *
				(pt_fourier[0        ][m][k][0] +
				 pt_fourier[nedge - 2][m][k][0]);
			}
			for (size_t l = 1; l < norder; l++) {
				feature[((2 * k - 1) * norder + l + 1) * nevent + (nevent_read + i)] =
					feature[((2 * k - 1) * norder + l) * nevent + (nevent_read + i)] *
					feature[((2 * k - 1) * norder + 1) * nevent + (nevent_read + i)];
			}
			feature[((2 * k    ) * norder + 1) * nevent + (nevent_read + i)] = 0;
			for (size_t m = 0; m < 3; m++) {
			feature[((2 * k    ) * norder + 1) * nevent + (nevent_read + i)] += scale[k] *
				(pt_fourier[0        ][m][k][1] +
				 pt_fourier[nedge - 2][m][k][1]);
			}
			for (size_t l = 1; l < norder; l++) {
				feature[((2 * k    ) * norder + l + 1) * nevent + (nevent_read + i)] =
					feature[((2 * k    ) * norder + l) * nevent + (nevent_read + i)] *
					feature[((2 * k    ) * norder + 1) * nevent + (nevent_read + i)];
			}
			scale_check_table[k].push_back(std::max(feature[((2 * k - 1) * norder + 1) * nevent + (nevent_read + i)], feature[((2 * k    ) * norder + 1) * nevent + (nevent_read + i)]) / scale[k]);
		}

		for (size_t j = 1; j < nedge; j++) {
			for (size_t k = 0; k < nreduced_id; k++) {
				const size_t column_0 = ((j - 1) * nreduced_id + k) * nfeature;

				target[column_0 * target_size + (nevent_read + i)] = pt_fourier[j - 1][k][0][0];

				for (size_t l = 1; l < nfourier; l++) {
					for (size_t re_or_im = 0; re_or_im < 2;
						 re_or_im++) {
						const size_t column = column_0 + (2 * l - 1 + re_or_im);

						target[column * target_size + (nevent_read + i)] = pt_fourier[j - 1][k][l][re_or_im];
					}
				}
			}
		}

		for (size_t l = 0; l < nfeature; l++) {
			if (!std::isfinite(feature[(l) * nevent + (nevent_read + i)])) {
				for (size_t k = 0; k < nfeature; k++) {
					fprintf(stderr, "%g ", feature[(k) * nevent + (nevent_read + i)]);
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
#endif
	}
		nevent_read += nevent_file;
		f.Close();
	}
#if 0
	fprintf(stderr, "\r                                             "
			"                                  \n");
#endif
	fprintf(stderr, "%s:%d: nevent = %lu\n", __FILE__, __LINE__, nevent);

	for (size_t i = 0; i < nfourier; i++) {
		std::sort(scale_check_table[i].begin(), scale_check_table[i].end());

		const size_t index_95 = static_cast<size_t>(rint((scale_check_table[i].size() - 1) * 0.95));
		const size_t index_99 = static_cast<size_t>(rint((scale_check_table[i].size() - 1) * 0.99));

		fprintf(stderr, "%s:%d: %lu %f %f\n", __FILE__, __LINE__, i, scale_check_table[i][index_95], scale_check_table[i][index_99]);
	}


#if 1
	std::sort(density_table.begin(), density_table.end());

	for (size_t i = 0; i < nevent; i++) {
		const double e_hf = feature[(1) * nevent + i] / scale[0];
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

		for (size_t j = 0; j < ncolumn; j++) {
			feature_copy[j * nevent + i] = feature[j * nevent + i];
			feature[j * nevent + i] *= weight;
		}
		for (size_t j = 0; j < nrhs; j++) {
			target_copy[j * target_size + i] = target[j * target_size + i];
			target[j * target_size + i] *= weight;
		}
	}
#endif

	std::vector<double> singular_value(
		std::min(nevent, ncolumn),
		0);
	std::vector<double> working_array(
		3 * std::min(nevent, ncolumn) +
		std::max(std::max(2 * std::min(nevent, ncolumn),
						  std::max(nevent, ncolumn)),
				 nrhs),
		0);
	double rank_condition = 0;

	int nevent_int = nevent;
	int ncolumn_int = ncolumn;
	int nrhs_int = nrhs;
	int target_size_int = target_size;
	int working_array_size = working_array.size();
	int rank;
	int info;

	dgelss_(// m, n, nrhs
			&nevent_int, &ncolumn_int, &nrhs_int,
			// a, lda
			&feature[0], &nevent_int,
			// b, ldb
			&target[0], &target_size_int,
			// s, rcond, rank
			&singular_value[0], &rank_condition, &rank,
			// work, lwork, info
			&working_array[0], &working_array_size, &info);

	fprintf(stderr, "%s:%d: DGELSS info = %d, rank = %d\n", __FILE__, __LINE__, info, rank);

	if (info != 0) {
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


	if (true) {
		for (size_t j = 2; j < 2 + 1; j++) {
			for (size_t k = 0; k < 1; k++) {
				const size_t column_0 = ((j - 1) * nreduced_id + k) * nfeature;

				// for (size_t l = 1; l < nfourier; l++) {
				// 	for (size_t re_or_im = 0; re_or_im < 2;
				// 		 re_or_im++) {
				// const size_t column = column_0 + (2 * l - 1 + re_or_im);

				for (size_t i = 0; i < nevent; i++) {
					double s = 0;

					for (size_t l = 0; l < ncolumn; l++) {
						s += target[column_0 * target_size_int + l] * feature_copy[l * nevent + i];
					}
					fprintf(stderr, "%s:%d: %lu %.8e %.8e %.8e\n", __FILE__, __LINE__, i, feature_copy[nevent + i] / scale[0], target_copy[column_0 * target_size + i], s);

					s = 0;

					for (size_t l = 0; l < ncolumn; l++) {
						s += target[column_0 * target_size_int + l] * feature_copy[l * nevent + i] * 2;
					}
					fprintf(stderr, "%s:%d: %lu %.8e %.8e %.8e\n", __FILE__, __LINE__, i, feature_copy[nevent + i] / scale[0] * 2, target_copy[column_0 * target_size + i] * 2, s);
				}
			}
		}
	}

	return EXIT_SUCCESS;
}

#include <fstream>
#include <sstream>
#include <algorithm>

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TMath.h>
#include <TRandom3.h>

#define NOBJECT_MAX 16384

double ue_predictor_pf[3][15][5][2][82];
double ue_interpolation_pf0[15][344];
double ue_interpolation_pf1[15][344];
double ue_interpolation_pf2[15][82];

namespace {

	double angular_range_reduce(const double x)
	{
		if (!std::isfinite(x))
			return x;

		static const double cody_waite_x_max = 1608.4954386379741381;
		static const double two_pi_0 = 6.2831853071795649157;
		static const double two_pi_1 = 2.1561211432631314669e-14;
		static const double two_pi_2 = 1.1615423895917441336e-27;
		double ret;

		if (x >= -cody_waite_x_max && x <= cody_waite_x_max) {
			static const double inverse_two_pi =
				0.15915494309189534197;
			const double k = rint(x * inverse_two_pi);
			ret = ((x - (k * two_pi_0)) - k * two_pi_1) -
				k * two_pi_2;
		}
		else {
			long double sin_x;
			long double cos_x;

			sincosl(x, &sin_x, &cos_x);
			ret = (double)atan2l(sin_x, cos_x);
		}
		if (ret == -M_PI) {
			ret = M_PI;
		}

		return ret;
	}

size_t pf_id_reduce(const Int_t pf_id)
{
	return 0;
}

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

	typedef CGAL::Delaunay_triangulation_2<
		CGAL::Exact_predicates_inexact_constructions_kernel>
	delaunay_triangulation_t;
	typedef delaunay_triangulation_t::Point point_2d_t;
	typedef CGAL::Voronoi_diagram_2<
		delaunay_triangulation_t,
		CGAL::Delaunay_triangulation_adaptation_traits_2<
			delaunay_triangulation_t>,
		CGAL::
Delaunay_triangulation_caching_degeneracy_removal_policy_2<
			delaunay_triangulation_t> > voronoi_diagram_t;
	typedef CGAL::Polygon_2<
		CGAL::Exact_predicates_inexact_constructions_kernel>
		polygon_t;
	typedef CGAL::Polygon_with_holes_2<
		CGAL::Exact_predicates_inexact_constructions_kernel>
		polygon_hole_t;

	// Obtain the Voronoi geometry information

	void voronoi_area_diameter(
		std::vector<double> &particle_area,
		std::vector<double> &particle_diameter_square,
		std::vector<std::set<size_t> > &particle_incident,
		const std::vector<point_2d_t> &particle_pseudorapidity_azimuth
		)
	{
		// Make the Voronoi diagram

		voronoi_diagram_t diagram;

		// Reverse Voronoi face lookup
		std::map<voronoi_diagram_t::Face_handle, size_t> face_index;

	static const size_t ncms_hcal_edge_pseudorapidity = 82 + 1;
	static const double cms_hcal_edge_pseudorapidity[
		ncms_hcal_edge_pseudorapidity] = {
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

	const std::vector<double> cms_hcal_edge_pseudorapidity_v(
		cms_hcal_edge_pseudorapidity,
		cms_hcal_edge_pseudorapidity + ncms_hcal_edge_pseudorapidity);

		for (std::vector<point_2d_t>::const_iterator iterator =
				 particle_pseudorapidity_azimuth.begin();
			 iterator != particle_pseudorapidity_azimuth.end();
			 iterator++) {
			// Make two additional replicas with azimuth +/- 2 pi (and
			// use only the middle) to mimick the azimuthal cyclicity
			for (int j = -1; j <= 1; j++) {
				const double parity = j == 0 ? 1 : -1;
				const double origin = j == 0 ? 0 :
					j == -1 ? 2 * cms_hcal_edge_pseudorapidity_v.front() :
					2 * cms_hcal_edge_pseudorapidity_v.back();
				for (int k = -1; k <= 1; k++) {
					const point_2d_t p(origin + parity * iterator->x(),
									   k * (2 * M_PI) + iterator->y());
					const voronoi_diagram_t::Face_handle
						handle = diagram.insert(p);

					face_index[handle] = iterator -
						particle_pseudorapidity_azimuth.begin();
				}
			}
		}

		particle_area.clear();
		particle_diameter_square.clear();
		particle_incident = std::vector<std::set<size_t> >(
			particle_pseudorapidity_azimuth.size(),
			std::set<size_t>());

		// Extract the Voronoi cells as polygon and calculate the
		// area associated with individual particles

		for (std::vector<point_2d_t>::const_iterator iterator =
				 particle_pseudorapidity_azimuth.begin();
			 iterator != particle_pseudorapidity_azimuth.end();
			 iterator++) {
			const voronoi_diagram_t::Locate_result result =
				diagram.locate(*iterator);
			const voronoi_diagram_t::Face_handle *face =
				boost::get<voronoi_diagram_t::Face_handle>(&result);
			double polygon_area;
			double polygon_diameter_square;

			if (face != NULL) {
				voronoi_diagram_t::Ccb_halfedge_circulator
					circulator_start = (*face)->outer_ccb();
				bool unbounded = false;
				polygon_t polygon;

				voronoi_diagram_t::Ccb_halfedge_circulator
					circulator = circulator_start;

				// Circle around the edges and extract the polygon
				// vertices
				do {
					if (circulator->has_target()) {
						polygon.push_back(
							circulator->target()->point());
						particle_incident[face_index[*face]].insert(
							face_index[circulator->twin()->face()]);
					}
					else {
						unbounded = true;
						break;
					}
				}
				while (++circulator != circulator_start);
				if (unbounded) {
					polygon_area = INFINITY;
					polygon_diameter_square = INFINITY;
				}
				else {
					polygon_area = polygon.area();
					polygon_diameter_square = 0;

					for (polygon_t::Vertex_iterator
							 iterator_vertex_outer =
							 polygon.vertices_begin();
						 iterator_vertex_outer !=
							 polygon.vertices_end();
						 iterator_vertex_outer++) {
						for (polygon_t::Vertex_iterator
								 iterator_vertex_inner =
								 iterator_vertex_outer + 1;
							 iterator_vertex_inner !=
								 polygon.vertices_end();
							 iterator_vertex_inner++) {
							polygon_diameter_square = std::max(
								polygon_diameter_square,
								squared_distance(
									*iterator_vertex_inner,
									*iterator_vertex_outer));
						}
					}

				}
			}
			else {
				polygon_area = NAN;
				polygon_diameter_square = NAN;
			}
			particle_area.push_back(fabs(polygon_area));
			particle_diameter_square.push_back(
				polygon_diameter_square);
		}
	}

fastjet::PseudoJet pseudo_jet(
	const double perp, const double pseudorapidity,
	const double azimuth, const double mass = 0)
{
	fastjet::PseudoJet r;

	r.reset_PtYPhiM(perp, pseudorapidity, azimuth, mass);

	return r;
}

void antikt_cluster(
	float *jet_perp, float *jet_pseudorapidity,
	float *jet_azimuth, int &njet,
	const std::vector<fastjet::PseudoJet> particle_positive,
	const std::vector<fastjet::PseudoJet> particle,
	float *pf_particle,
	fastjet::JetDefinition jet_definition)
{
	const fastjet::ClusterSequence
		cluster_sequence(particle_positive, jet_definition);
	const std::vector<fastjet::PseudoJet> jet =
		fastjet::sorted_by_pt(cluster_sequence.inclusive_jets(0));

	njet = jet.size();
	for (int i = 0; i < njet; i++) {
		std::vector<fastjet::PseudoJet> constituent =
			cluster_sequence.constituents(jet[i]);
		fastjet::PseudoJet jet_resum(0, 0, 0, 0);
		// double perp_resum = 0;

		for (std::vector<fastjet::PseudoJet>::const_iterator
				 iterator_constituent = constituent.begin();
			 iterator_constituent != constituent.end();
			 iterator_constituent++) {
			jet_resum += particle[iterator_constituent->user_index()];
			// perp_resum += pf_particle[iterator_constituent->user_index()];
		}

		jet_perp[i] = jet_resum.perp();
		jet_pseudorapidity[i] = jet_resum.pseudorapidity();
		jet_azimuth[i] = jet_resum.phi_std();
	}
}

}

void voronoi_recalc(const char *filename, const int data,
					const int calorimetric, const double antikt_distance)
{
	static const size_t nfourier = 5;

	const char *root_tree_name = calorimetric ?
		"rechitanalyzer/tower" : "pfcandAnalyzer/pfTree";

	TFile *root_file = TFile::Open(filename);

	if (root_file == NULL) {
		return;
	}

	TTree *root_tree = dynamic_cast<TTree *>(gDirectory->Get(root_tree_name));

	Int_t nPFpart;
	Int_t pfId[NOBJECT_MAX];
	Float_t pfPt[NOBJECT_MAX];
	Float_t pfVsPtInitial[NOBJECT_MAX];
	Float_t pfVsPtInitialRecalc[NOBJECT_MAX];
	Float_t pfEta[NOBJECT_MAX];
	Float_t pfPhi[NOBJECT_MAX];
	Float_t pfArea[NOBJECT_MAX];

	if (calorimetric) {
		root_tree->SetBranchAddress("n", &nPFpart);
		root_tree->SetBranchAddress("et", pfPt);
		root_tree->SetBranchAddress("eta", pfEta);
		root_tree->SetBranchAddress("phi", pfPhi);
		root_tree->SetBranchAddress("vsPtInitial", pfVsPtInitial);
		root_tree->SetBranchAddress("vsArea", pfArea);
	}
	else {
		root_tree->SetBranchAddress("nPFpart", &nPFpart);
		root_tree->SetBranchAddress("pfId", pfId);
		root_tree->SetBranchAddress("pfPt", pfPt);
		root_tree->SetBranchAddress("pfVsPtInitial", pfVsPtInitial);
		root_tree->SetBranchAddress("pfEta", pfEta);
		root_tree->SetBranchAddress("pfPhi", pfPhi);
		root_tree->SetBranchAddress("pfArea", pfArea);
	}

	fprintf(stderr, "%s:%d:\n", __FILE__, __LINE__);

	TTree *hiTree = dynamic_cast<TTree *>(gDirectory->Get("hiEvtAnalyzer/HiTree"));

	Int_t run;
	Int_t lumi;
	Float_t vx;
	Float_t vy;
	Float_t vz;
	Int_t hiBin;
	Float_t hiHF;
	Float_t hiHFplus;
	Float_t hiHFminus;
	Float_t hiZDC;
	Float_t hiZDCplus;
	Float_t hiZDCminus;
	Float_t hiHFhit;
	Float_t hiHFhitPlus;
	Float_t hiHFhitMinus;
	Float_t hiET;
	Float_t hiEE;
	Float_t hiEB;
	Float_t hiEEplus;
	Float_t hiEEminus;
	Int_t hiNpix;
	Int_t hiNpixelTracks;
	Int_t hiNtracks;
	Int_t hiNtracksPtCut;
	Int_t hiNtracksEtaCut;
	Int_t hiNtracksEtaPtCut;
	Int_t hiNevtPlane;
	Float_t hiEvtPlanes[38];

	if (hiTree != NULL) {
	// Set branch addresses.
	hiTree->SetBranchAddress("run",&run);
	hiTree->SetBranchAddress("lumi",&lumi);
	hiTree->SetBranchAddress("vx",&vx);
	hiTree->SetBranchAddress("vy",&vy);
	hiTree->SetBranchAddress("vz",&vz);
	hiTree->SetBranchAddress("hiBin",&hiBin);
	hiTree->SetBranchAddress("hiHF",&hiHF);
	hiTree->SetBranchAddress("hiHFplus",&hiHFplus);
	hiTree->SetBranchAddress("hiHFminus",&hiHFminus);
	hiTree->SetBranchAddress("hiZDC",&hiZDC);
	hiTree->SetBranchAddress("hiZDCplus",&hiZDCplus);
	hiTree->SetBranchAddress("hiZDCminus",&hiZDCminus);
	hiTree->SetBranchAddress("hiHFhit",&hiHFhit);
	hiTree->SetBranchAddress("hiHFhitPlus",&hiHFhitPlus);
	hiTree->SetBranchAddress("hiHFhitMinus",&hiHFhitMinus);
	hiTree->SetBranchAddress("hiET",&hiET);
	hiTree->SetBranchAddress("hiEE",&hiEE);
	hiTree->SetBranchAddress("hiEB",&hiEB);
	hiTree->SetBranchAddress("hiEEplus",&hiEEplus);
	hiTree->SetBranchAddress("hiEEminus",&hiEEminus);
	hiTree->SetBranchAddress("hiNpix",&hiNpix);
	hiTree->SetBranchAddress("hiNpixelTracks",&hiNpixelTracks);
	hiTree->SetBranchAddress("hiNtracks",&hiNtracks);
	hiTree->SetBranchAddress("hiNtracksPtCut",&hiNtracksPtCut);
	hiTree->SetBranchAddress("hiNtracksEtaCut",&hiNtracksEtaCut);
	hiTree->SetBranchAddress("hiNtracksEtaPtCut",&hiNtracksEtaPtCut);
	hiTree->SetBranchAddress("hiNevtPlane",&hiNevtPlane);
	hiTree->SetBranchAddress("hiEvtPlanes",hiEvtPlanes);
	}

	char tree_name[BUFSIZ];

	snprintf(tree_name, BUFSIZ, "akVs%d%sJetAnalyzer/t", static_cast<int>(rint(antikt_distance * 10.0)), calorimetric ? "Calo" : "PF");

	TTree *t = dynamic_cast<TTree *>(gDirectory->Get(tree_name));

	Int_t evt;
	// Float_t b;
	Int_t nref;
	Float_t rawpt[NOBJECT_MAX];
	Float_t jtpt[NOBJECT_MAX];
	Float_t jteta[NOBJECT_MAX];
	Float_t jty[NOBJECT_MAX];
	Float_t jtphi[NOBJECT_MAX];
	Float_t jtpu[NOBJECT_MAX];
	Int_t beamId1;
	Int_t beamId2;
	Float_t pthat;
	Float_t refpt[NOBJECT_MAX];
	Float_t refeta[NOBJECT_MAX];
	Float_t refy[NOBJECT_MAX];
	Float_t refphi[NOBJECT_MAX];
	Float_t refdphijt[NOBJECT_MAX];
	Float_t refdrjt[NOBJECT_MAX];
	Float_t refparton_pt[NOBJECT_MAX];
	Int_t refparton_flavor[NOBJECT_MAX];
	Int_t refparton_flavorForB[NOBJECT_MAX];
	Int_t subid[NOBJECT_MAX];

	Int_t ngen;
	Int_t genmatchindex[NOBJECT_MAX];
	Float_t genpt[NOBJECT_MAX];
	Float_t geneta[NOBJECT_MAX];
	Float_t geny[NOBJECT_MAX];
	Float_t genphi[NOBJECT_MAX];
	Float_t gendphijt[NOBJECT_MAX];
	Float_t gendrjt[NOBJECT_MAX];
	Int_t gensubid[NOBJECT_MAX];

	// Set branch addresses.
	t->SetBranchAddress("evt", &evt);
	// t->SetBranchAddress("b", &b);
	t->SetBranchAddress("nref", &nref);
	t->SetBranchAddress("rawpt", rawpt);
	t->SetBranchAddress("jtpt", jtpt);
	t->SetBranchAddress("jteta", jteta);
	t->SetBranchAddress("jty", jty);
	t->SetBranchAddress("jtphi", jtphi);
	t->SetBranchAddress("jtpu", jtpu);
	t->SetBranchAddress("beamId1", &beamId1);
	t->SetBranchAddress("beamId2", &beamId2);
	t->SetBranchAddress("pthat", &pthat);
	t->SetBranchAddress("refpt", refpt);
	t->SetBranchAddress("refeta", refeta);
	t->SetBranchAddress("refy", refy);
	t->SetBranchAddress("refphi", refphi);
	t->SetBranchAddress("refdphijt", refdphijt);
	t->SetBranchAddress("refdrjt", refdrjt);
	t->SetBranchAddress("refparton_pt", refparton_pt);
	t->SetBranchAddress("refparton_flavor", refparton_flavor);
	t->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);
	t->SetBranchAddress("subid", subid);

	t->SetBranchAddress("ngen",&ngen);
	t->SetBranchAddress("genmatchindex",genmatchindex);
	t->SetBranchAddress("genpt",genpt);
	t->SetBranchAddress("geneta",geneta);
	t->SetBranchAddress("geny",geny);
	t->SetBranchAddress("genphi",genphi);
	t->SetBranchAddress("gendphijt",gendphijt);
	t->SetBranchAddress("gendrjt",gendrjt);
	t->SetBranchAddress("gensubid",gensubid);

	fprintf(stderr, "%s:%d: data = %d, calorimetric = %d\n", __FILE__, __LINE__, data, calorimetric);

	std::ifstream in_stream(
		(data ? (calorimetric ?
				 "ue_calibrations_calo_data.txt" :
				 "ue_calibrations_pf_data.txt") :
		 		(calorimetric ?
				 "ue_calibrations_calo_mc.txt" :
				 "ue_calibrations_pf_mc.txt"))
		);
	std::string line;
	size_t index = 0;
	const size_t nline_predictor = 3 * 15 * (1 + (5 - 1) * 2) * 82;

	fprintf(stderr, "%s:%d:\n", __FILE__, __LINE__);

	TRandom3 rand;

	while (std::getline(in_stream, line)) {
		if (line.empty() || line[0] == '#') {
			continue;
		}

		std::istringstream line_stream(line);
		double val;
		int bin0, bin1, bin2, bin3, bin4;

		if (index < nline_predictor) {
			line_stream >> bin0 >> bin1 >> bin2 >> bin3 >> bin4 >> val;
			ue_predictor_pf[bin0][bin1][bin2][bin3][bin4] = val;
		}
		else if (index < nline_predictor + sizeof(ue_interpolation_pf0) / sizeof(double)) {
			line_stream >> bin0 >> bin1 >> val;
			ue_interpolation_pf0[bin0][bin1] = val;
		}
		else if (index < nline_predictor + sizeof(ue_interpolation_pf0) / sizeof(double) + sizeof(ue_interpolation_pf1) / sizeof(double)) {
            line_stream >> bin0 >> bin1 >> val;
            ue_interpolation_pf1[bin0][bin1] = val;
		}
		else if (index < nline_predictor + sizeof(ue_interpolation_pf0) / sizeof(double) + sizeof(ue_interpolation_pf1) / sizeof(double) + sizeof(ue_interpolation_pf2) / sizeof(double)) {
			line_stream >> bin0 >> bin1 >> val;
			ue_interpolation_pf2[bin0][bin1] = val;
		}
		index++;
	}

	fprintf(stderr, "%s:%d:\n", __FILE__, __LINE__);

	static const size_t nreduced_id = 3;

	static const size_t nedge_pseudorapidity = 15 + 1;
	static const double edge_pseudorapidity[nedge_pseudorapidity] = {
		-5.191, -2.650, -2.043, -1.740, -1.479, -1.131, -0.783, -0.522,
		0.522, 0.783, 1.131, 1.479, 1.740, 2.043, 2.650, 5.191
	};

	const std::vector<double> edge_pseudorapidity_v(
		edge_pseudorapidity,
		edge_pseudorapidity + nedge_pseudorapidity);

	static const size_t ncms_hcal_edge_pseudorapidity = 82 + 1;
	static const double cms_hcal_edge_pseudorapidity[
		ncms_hcal_edge_pseudorapidity] = {
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

	const std::vector<double> cms_hcal_edge_pseudorapidity_v(
		cms_hcal_edge_pseudorapidity,
		cms_hcal_edge_pseudorapidity + ncms_hcal_edge_pseudorapidity);

	static const size_t ncms_ecal_edge_pseudorapidity = 344 + 1;
	double cms_ecal_edge_pseudorapidity[
		ncms_ecal_edge_pseudorapidity];

	for (size_t i = 0; i < ncms_ecal_edge_pseudorapidity; i++) {
		cms_ecal_edge_pseudorapidity[i] =
			i * (2 * 2.9928 /
				 (ncms_ecal_edge_pseudorapidity - 1)) -
			2.9928;
	};

	const std::vector<double> cms_ecal_edge_pseudorapidity_v(
		cms_ecal_edge_pseudorapidity,
		cms_ecal_edge_pseudorapidity + ncms_ecal_edge_pseudorapidity);

	size_t nentries = root_tree->GetEntries();
	//nentries = 50;

	fastjet::JetDefinition jet_definition(fastjet::antikt_algorithm,
										  antikt_distance);

	fprintf(stderr, "%s:%d:\n", __FILE__, __LINE__);

	for (size_t i = 0; i < nentries; i++) {
		if (root_tree != NULL) {
			root_tree->GetEntry(i);
		}
		if (hiTree != NULL) {
			hiTree->GetEntry(i);
		}
		if (t != NULL) {
			t->GetEntry(i);
		}

		fprintf(stderr, "%s:%d: %d\n", __FILE__, __LINE__, hiBin);

		if ((hiBin * 0.5 >= 10 && hiBin * 0.5 < 30 && i % 3 != 0) ||
			(hiBin * 0.5 >= 50 && i % 5 != 0)) {
			continue;
		}

			std::vector<point_2d_t> particle_pseudorapidity_azimuth;

			for (Int_t k = 0; k < nPFpart; k++) {
				particle_pseudorapidity_azimuth.push_back(point_2d_t(
						pfEta[k], pfPhi[k]));
			}

			std::vector<double> particle_area;
			std::vector<double> particle_diameter_square;
			std::vector<std::set<size_t> > particle_incident;

			voronoi_area_diameter(particle_area, particle_diameter_square,
								  particle_incident,
								  particle_pseudorapidity_azimuth);

			for (Int_t k = 0; k < nPFpart; k++) {
				pfArea[k] = particle_area[k];
			}

		// Event collective Fourier components, per particle flow ID
		// group. Note that since by heavy ion convention, dN/dphi =
		// v0/(2 pi) + v0 v1/pi cos(phi - Psi_RP) + v0 v2/pi cos(2(phi
		// - Psi_RP)) + ..., and orthonormal relation for the Fourier
		// basis f0 = v0, f1 = v0 v1, ..., if f is the Fourier and v0
		// the phi-averaged 1/(2 pi) dpT/dy

		double perp_fourier[nedge_pseudorapidity - 1][nreduced_id][nfourier][2];

		for (size_t k = 1; k < nedge_pseudorapidity; k++) {
			for (size_t l = 0; l < nreduced_id; l++) {
				for (size_t m = 0; m < nfourier; m++) {
					for (size_t re_or_im = 0; re_or_im < 2;
						 re_or_im++) {
						perp_fourier[k - 1][l][m][re_or_im] = 0;
					}
				}
			}
		}
		memset(perp_fourier, 0,
			   (nedge_pseudorapidity - 1) * nreduced_id * nfourier * 2 *
			   sizeof(double));

		for (Int_t j = 0; j < nPFpart; j++) {
			size_t reduced_id = pf_id_reduce(pfId[j]);

			for (size_t k = 1; k < nedge_pseudorapidity; k++) {
				if (pfEta[j] >= edge_pseudorapidity[k - 1] &&
					pfEta[j] < edge_pseudorapidity[k]) {
					for (size_t l = 0; l < nfourier; l++) {
						perp_fourier[k - 1][reduced_id][l][0] +=
							pfPt[j] * cos(l * pfPhi[j]);
						perp_fourier[k - 1][reduced_id][l][1] +=
							pfPt[j] * sin(l * pfPhi[j]);
					}
				}
			}
		}

		// Event selection

		static const size_t nfeature = 2 * nfourier - 1;
		double feature[nfeature];

		// Scale factor to get 95% of the coefficient below 1.0

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

		feature[0] = 0;
		for (size_t l = 0; l < nreduced_id; l++) {
		feature[0] += scale[0] *
			(perp_fourier[0                       ][l][0][0] +
			 perp_fourier[nedge_pseudorapidity - 2][l][0][0]);
		}
		for (size_t k = 1; k < nfourier; k++) {
			feature[2 * k - 1] = 0;
			for (size_t l = 0; l < nreduced_id; l++) {
			feature[2 * k - 1] += scale[k] *
				(perp_fourier[0                       ][l][k][0] +
				 perp_fourier[nedge_pseudorapidity - 2][l][k][0]);
			}
			feature[2 * k] = 0;
			for (size_t l = 0; l < nreduced_id; l++) {
			feature[2 * k] += scale[k] *
				(perp_fourier[0                       ][l][k][1] +
				 perp_fourier[nedge_pseudorapidity - 2][l][k][1]);
			}
		}

#if 0
		const double event_plane_2 = atan2(feature[4], feature[3]);
		const double v2 =
			sqrt(feature[3] * feature[3] +
				 feature[4] * feature[4]) / feature[0];
#endif

		for (Int_t k = 0; k < nPFpart; k++) {
			int predictor_index = -1;
			int interpolation_index = -1;
			double density = 0;

			if (pfEta[k] >= edge_pseudorapidity[0] &&
				pfEta[k] < edge_pseudorapidity[nedge_pseudorapidity - 1]) {
				std::vector<double>::const_iterator p = std::lower_bound(edge_pseudorapidity_v.begin(), edge_pseudorapidity_v.end(), pfEta[k]);

				predictor_index = (p - edge_pseudorapidity_v.begin()) - 1;
			}

			for (size_t j = 0; j < nreduced_id; j++) {
				if (j == 2) {
					// HCAL
					if (pfEta[k] >=
						cms_hcal_edge_pseudorapidity[0] &&
						pfEta[k] <
						cms_hcal_edge_pseudorapidity[ncms_hcal_edge_pseudorapidity - 1]) {
						std::vector<double>::const_iterator p = std::lower_bound(cms_hcal_edge_pseudorapidity_v.begin(), cms_hcal_edge_pseudorapidity_v.end(), pfEta[k]);

						interpolation_index = (p - cms_hcal_edge_pseudorapidity_v.begin()) - 1;
					}
				}
				else {
					// Tracks or ECAL clusters
					if (pfEta[k] >=
						cms_ecal_edge_pseudorapidity[0] &&
						pfEta[k] <
						cms_ecal_edge_pseudorapidity[ncms_ecal_edge_pseudorapidity - 1]) {
						std::vector<double>::const_iterator p = std::lower_bound(cms_ecal_edge_pseudorapidity_v.begin(), cms_ecal_edge_pseudorapidity_v.end(), pfEta[k]);

						interpolation_index = (p - cms_ecal_edge_pseudorapidity_v.begin()) - 1;
					}
				}

				if (predictor_index >= 0 && interpolation_index >= 0) {
					// Calculate the aggregated prediction and
					// interpolation for the pseudorapidity segment

					const double azimuth = pfPhi[k];
					const double (*p)[2][82] =
						ue_predictor_pf[j][predictor_index];
					double pred = 0;

					for (size_t l = 0; l < nfourier; l++) {
						const size_t norder = l == 0 ? 9 : 1;

						for (size_t m = 0; m < 2; m++) {
							float u = p[l][m][0];

							for (size_t n = 0; n < 2 * nfourier - 1; n++) {
								if ((l == 0 && n == 0) || (l == 2 && (n == 3 || n == 4))) {
									u += p[l][m][9 * n + 1] * feature[n];
									for (size_t o = 2; o < norder + 1; o++) {
										u += p[l][m][9 * n + o] * hermite_h_normalized(
											2 * o - 1, feature[n]) *
											exp(-feature[n] * feature[n]);
									}
								}
							}
#if 0
							// This looks at a specific flow component and see how the polynomial is evaluated
							if (j == 0 && predictor_index == 3 && l == 0 && m == 0) {
								//fprintf(stderr, "%s:%d: %f %f\n", __FILE__, __LINE__, perp_fourier[0][2][2][0], perp_fourier[nedge_pseudorapidity - 2][2][2][1]);
								fprintf(stderr, "%s:%d: << %f %f %f %f %f %f %f\n", __FILE__, __LINE__, feature[0], feature[1], feature[2], feature[3], feature[4], u, perp_fourier[predictor_index][j][l][m]);
							}
#endif

							pred += u * (l == 0 ? 1.0 : 2.0) *
								(m == 0 ? cos(l * azimuth) :
								 sin(l * azimuth));
						}
					}

					double interp;

					if (j == 0) {
						interp =
							ue_interpolation_pf0[predictor_index][
								interpolation_index];
					}
					else if (j == 1) {
						interp =
							ue_interpolation_pf1[predictor_index][
								interpolation_index];
					}
					else if (j == 2) {
						interp =
							ue_interpolation_pf2[predictor_index][
								interpolation_index];
					}

					// Interpolate down to the finely binned
					// pseudorapidity

					density += pred /
						(2.0 * M_PI *
						 (edge_pseudorapidity[predictor_index + 1] -
						  edge_pseudorapidity[predictor_index])) *
						interp;
				}
			}

							// Prints the subtracted density * area
			//fprintf(stderr, "%s:%d: %.8e %.8e %.8e %.8e %.8e\n", __FILE__, __LINE__, hiBin * 0.5, pfEta[k], pfPhi[k], pfPt[k], density * pfArea[k]);
			pfVsPtInitialRecalc[k] = pfPt[k] - (TMath::Finite(pfArea[k]) ? density * pfArea[k] : 0);
		}

/////////////////////////////////////////////////////////////////////
		for (Int_t j = 0; j < 16; j++) {
			float cone_pseudorapidity = rand.Uniform(-5.191, 5.191);
			float cone_azimuth = rand.Uniform(-M_PI, M_PI);
			float perp_sum_unsubtracted = 0;
			float perp_sum = 0;
			float perp_sum_recalc = 0;

			for (Int_t k = 0; k < nPFpart; k++) {
				if (std::pow(pfEta[k] - cone_pseudorapidity, 2) +
					std::pow(angular_range_reduce(pfPhi[k] - cone_azimuth), 2) < antikt_distance * antikt_distance) {
					perp_sum_unsubtracted += pfPt[k];
					perp_sum += pfVsPtInitial[k];
					perp_sum_recalc += pfVsPtInitialRecalc[k];
				}
			}

			fprintf(stderr, "%s:%d: %f %f %f %f %f\n", __FILE__, __LINE__, perp_sum_unsubtracted, perp_sum, perp_sum_recalc, cone_pseudorapidity, cone_azimuth);
		}

		const float infinitesimalPt = 0.005;

		std::vector<fastjet::PseudoJet> pf_particle_pseudojet;
		std::vector<fastjet::PseudoJet> pf_particle_pseudojet_positive;

		for (Int_t k = 0; k < nPFpart; k++) {
			fastjet::PseudoJet p =
				pseudo_jet(pfVsPtInitialRecalc[k], pfEta[k], pfPhi[k]);

			p.set_user_index(k);
			pf_particle_pseudojet.push_back(p);

			fastjet::PseudoJet pp =
				pseudo_jet(std::max(infinitesimalPt,
									pfVsPtInitialRecalc[k]),
						   pfEta[k], pfPhi[k]);

			pp.set_user_index(k);
			pf_particle_pseudojet_positive.push_back(pp);
		}

		int nref_old = nref;

		for (Int_t j = 0; j < nref; j++) {
			if (refpt[j] > 0) {
				fprintf(stderr, "%s:%d: %d %f %f %f %f\n", __FILE__, __LINE__, j, refpt[j], rawpt[j], jteta[j], jtphi[j]);
			}
		}

		antikt_cluster(rawpt, jteta, jtphi, nref,
					   pf_particle_pseudojet_positive,
					   pf_particle_pseudojet, pfVsPtInitialRecalc,
					   jet_definition);

		std::vector<bool> taken(nref, false);

		for (Int_t j = 0; j < nref; j++) {
			int index_max = -1;
			float perp_max = FLT_MIN;

			for (Int_t k = 0; k < ngen; k++) {
				if (std::pow(jteta[j] - geneta[k], 2) +
					std::pow(angular_range_reduce(
						jtphi[j] - genphi[k]), 2) <
					antikt_distance * antikt_distance) {
					if (genpt[k] > perp_max && gensubid[k] == 0 && !taken[k]) {
						index_max = k;
						perp_max = genpt[k];
					}
				}
			}

			if (index_max != -1) {
				refpt[j] = genpt[index_max];
				refeta[j] = geneta[index_max];
				refphi[j] = genphi[index_max];
				refparton_flavor[j] = -999;
				subid[j] = gensubid[index_max];

				taken[index_max] = true;
			}
			else {
				refpt[j] = -999;
				refeta[j] = -999;
				refphi[j] = -999;
				refparton_flavor[j] = -999;
				subid[j] = -999;
			}
		}
		fprintf(stderr, "%s:%d: %d %d %d\n", __FILE__, __LINE__, hiBin, nref, nref_old);
		for (Int_t j = 0; j < nref; j++) {
			if (refpt[j] > 0) {
				fprintf(stderr, "%s:%d: %d %f %f %f %f\n", __FILE__, __LINE__, j, refpt[j], rawpt[j], jteta[j], jtphi[j]);
			}
		}
/////////////////////////////////////////////////////////////////////

	}

	fprintf(stderr, "%s:%d:\n", __FILE__, __LINE__);

	root_file->Close();
}

int main(int argc, char *argv[])
{
	int data = 0;
	int calorimetric = 0;
	double antikt_distance = 0.4;

	if (argc > 2) {
		sscanf(argv[2], "%d", &data);
	}
	if (argc > 3) {
		sscanf(argv[3], "%d", &calorimetric);
	}
	if (argc > 4) {
		sscanf(argv[4], "%lf", &antikt_distance);
	}
	voronoi_recalc(argv[1], data, calorimetric, antikt_distance);

	return 0;
}

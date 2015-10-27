void copytree(char *filename_old = "/mnt/hadoop/cms/store/user/yetkin/MC_Production/HydjetDrum03/HydjetDrum03_HiForest_v05_merged_test02.root", char *filename_new = "hydjet_small.root", const unsigned int nevent = 256, int track_tree = 0)
{
	static const size_t BUFSIZ = 4096;
	TFile *root_file_old = new TFile(filename_old);
	TTree *pf_tree_old = (TTree *)root_file_old->Get("pfcandAnalyzer/pfTree");
	TTree *he_tree_old = (TTree *)root_file_old->Get("hiEvtAnalyzer/HiTree");
	TTree *hg_tree_old = (TTree *)root_file_old->Get("HiGenParticleAna/hi");
	TTree *tr_tree_old = NULL;

	if (track_tree != 0) {
		tr_tree_old = (TTree *)root_file_old->Get("anaTrack/trackTree");
	}

	TTree *t_old[6];

	for (unsigned int i = 0; i < 6; i++) {
		char buf[BUFSIZ];

		snprintf(buf, BUFSIZ, "akPu%uPFJetAnalyzer/t", i + 1);

		t_old[i] = (TTree *)root_file_old->Get(buf);
	}

	Long64_t nentries = pf_tree_old->GetEntries();

	fprintf(stderr, "%s:%d: nentries = %lu\n", __FILE__, __LINE__, nentries);

	TFile *root_file_new = new TFile(filename_new, "recreate");
	root_file_new->mkdir("pfcandAnalyzer");
	root_file_new->cd("pfcandAnalyzer");
	TTree *pf_tree_new = pf_tree_old->CloneTree(0);
	root_file_new->cd();
	root_file_new->mkdir("hiEvtAnalyzer");
	root_file_new->cd("hiEvtAnalyzer");
	TTree *he_tree_new = he_tree_old->CloneTree(0);
	TTree *hg_tree_new;
	TTree *tr_tree_new;

	if (hg_tree_old != NULL) {
		root_file_new->cd();
		root_file_new->mkdir("HiGenParticleAna");
		root_file_new->cd("HiGenParticleAna");
		hg_tree_new = hg_tree_old->CloneTree(0);
	}
	if (tr_tree_old != NULL) {
		root_file_new->cd();
		root_file_new->mkdir("anaTrack");
		root_file_new->cd("anaTrack");
		tr_tree_new = tr_tree_old->CloneTree(0);
	}

	TTree *t_new[6];

	for (unsigned int i = 0; i < 6; i++) {
		char buf[BUFSIZ];

		snprintf(buf, BUFSIZ, "akPu%uPFJetAnalyzer", i + 1);

		root_file_new->cd();
		root_file_new->mkdir(buf);
		root_file_new->cd(buf);
		t_new[i] = t_old[i]->CloneTree(0);
	}

	Float_t jtpt[65536];

	t_old[2]->SetBranchAddress("jtpt", jtpt);

	for (unsigned int i = 0; i < nevent; i++) {
		if (i % 100 == 0) {
			fprintf(stderr, "%s:%d: %u\n", __FILE__, __LINE__, i);
		}
#if 1
		t_old[2]->GetEntry(i);
		if (!(jtpt[0] >= 200)) {
			continue;
		}
#endif

		pf_tree_old->GetEntry(i);
		pf_tree_new->Fill();
		he_tree_old->GetEntry(i);
		he_tree_new->Fill();
		if (hg_tree_old != NULL) {
			hg_tree_old->GetEntry(i);
			hg_tree_new->Fill();
		}
		if (tr_tree_old != NULL) {
			tr_tree_old->GetEntry(i);
			tr_tree_new->Fill();
		}
		for (int j = 0; j < 6; j++) {
			t_old[j]->GetEntry(i);
			t_new[j]->Fill();
		}
	}
	pf_tree_new->AutoSave();
	he_tree_new->AutoSave();
	if (hg_tree_old != NULL) {
		hg_tree_new->AutoSave();
	}
	if (tr_tree_old != NULL) {
		tr_tree_new->AutoSave();
	}
	for (unsigned int i = 0; i < 6; i++) {
		t_new[i]->AutoSave();
	}

	delete root_file_old;
	delete root_file_new;

	gSystem->Exit(0);
}


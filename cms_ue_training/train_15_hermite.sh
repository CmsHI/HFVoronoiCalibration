#./cms_ue_train_fine_11 /mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/mc/PbPb_MinBias_Hydjet_STARTHI53_LV1_track8_Jet5_FOREST_6Feb2014_mergedv1.root > fine-mc-calo-11 2> log-fine-mc-calo-11 &
#./cms_ue_train_coarse_11 /mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/mc/PbPb_MinBias_Hydjet_STARTHI53_LV1_track8_Jet5_FOREST_6Feb2014_mergedv1.root > coarse-mc-calo-11 2> log-coarse-mc-calo-11 &
#nohup nice -19 ./cms_ue_train_fine_11 /mnt/hadoop/cms/store/user/belt/Validation53X/Track8_Jet11/*.root > fine-data-pf-11-t8j11 2> log-fine-data-pf-11-t8j11 &
#nohup nice -19 ./cms_ue_train_coarse_15_hermite /mnt/hadoop/cms/store/user/belt/Validation53X/Track8_Jet11/*.root > coarse-data-pf-15h-t8j11 2> log-coarse-data-pf-15h-t8j11 &

p="/mnt/hadoop/cms/store/user/dgulhan/HIMinBiasUPC-HIRun2011-14Mar2014-v2_tag_HI_MatchEqR_DatabaseJEC"
f="$p/HiForest_[0-9]_*.root $p/HiForest_[1-3][0-9]_*.root"

nohup nice -19 ./cms_ue_train_coarse_15_hermite $f > coarse-data-pf-15h-new 2> log-coarse-data-pf-15h-new &
nohup nice -19 ./cms_ue_train_coarse_15_hermite -c $f > coarse-data-calo-15h-new 2> log-coarse-data-calo-15h-new &

#exit 0

p="/mnt/hadoop/cms/store/user/ginnocen/Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV/HiMinBias_Forest_26June2014/d9ab4aca1923b3220eacf8ee0d550950"
f="$p/HiForest_[0-9]_*.root $p/HiForest_1[0-8]_*.root $p/HiForest_[2-9][0-9]_*.root $p/HiForest_[1-3][0-9][0-9]_*.root"

nohup nice -19 ./cms_ue_train_coarse_15_hermite $f > coarse-mc-pf-15h-new 2> log-coarse-mc-pf-15h-new &
nohup nice -19 ./cms_ue_train_coarse_15_hermite -c $f > coarse-mc-calo-15h-new 2> log-coarse-mc-calo-15h-new &


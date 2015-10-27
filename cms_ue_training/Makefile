#LAPACK =		/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/lapack/3.3.1
LAPACK  =               /cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/lapack/3.3.1-cms

CXX =			g++
CXXFLAGS =		-g -O2 -Wall -Werror $(shell root-config --cflags) -I.
LDFLAGS =		-L$(shell root-config --libdir) -L$(LAPACK)/lib -L$(LAPACK)/../../gcc/4.6.2/lib64
LIBS =			$(shell root-config --noldflags --libs) -llapack -lblas -lquadmath -lgfortran

all:			cms_ue_train_coarse cms_ue_train_fine \
			cms_ue_train_coarse_11 cms_ue_train_fine_11 \
			cms_ue_train_coarse_15_hermite

cms_ue_train_coarse:	cms_ue_train_coarse.o
			$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

cms_ue_train_fine:	cms_ue_train_fine.o
			$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

cms_ue_train_coarse_11:	cms_ue_train_coarse_11.o
			$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

cms_ue_train_fine_11:	cms_ue_train_fine_11.o
			$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

cms_ue_train_coarse_15_hermite:	cms_ue_train_coarse_15_hermite.o
			$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
			rm -f *~ *.o cms_ue_train_coarse cms_ue_train_fine

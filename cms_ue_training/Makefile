#LAPACK =		/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/lapack/3.3.1
LAPACK  =               /cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/lapack/3.3.1-cms
FASTJET =		/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/fastjet/3.1.0-odfocd
CGAL	=		/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/cgal/4.2-kpegke
BOOST	=		/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/boost/1.57.0-kpegke

CXX =			g++
CXXFLAGS =		-g -O2 -Wall -Werror -frounding-math $(shell root-config --cflags) -I$(CGAL)/include $(shell $(FASTJET)/bin/fastjet-config --cxxflags) -I.
LDFLAGS =		-L$(shell root-config --libdir) -L$(LAPACK)/lib -L$(LAPACK)/../../gcc/4.6.2/lib64 -L$(CGAL)/lib -L$(BOOST)/lib
LIBS =			$(shell root-config --noldflags --libs) -llapack -lblas -lquadmath -lgfortran

all:			cms_ue_train_coarse cms_ue_train_fine \
			cms_ue_train_coarse_11 cms_ue_train_fine_11 \
			cms_ue_train_coarse_15_hermite \
			cms_ue_recalc \
			cms_ue_train_coarse_15_hermite_datavsMCRun2

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

cms_ue_train_coarse_15_hermite_datavsMCRun2: cms_ue_train_coarse_15_hermite_datavsMCRun2.o
			$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

cms_ue_recalc:	cms_ue_recalc.o
			$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS) $(shell $(FASTJET)/bin/fastjet-config --libs) -lCGAL -lgmp -lboost_thread

clean:
			rm -f *~ *.o cms_ue_train_coarse cms_ue_train_fine

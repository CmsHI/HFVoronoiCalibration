# HFVoronoiCalibration
Code for underlying event calibration for HF-Voronoi jet algorithm.
For instructions, see also: https://twiki.cern.ch/twiki/bin/view/CMS/VoronoiUETraining
```
git clone git@githubSPAMNOT.com:CmsHI/HFVoronoiCalibration
```
## Building the code
The code requires a location to LAPACK, which is usually installed as an external with CMSSW. Normally you can just locate the LAPACK path by browsing through `/cvmfs/cms.cern.ch/<your architecture>`, or you can extract it from SCRAM:
```
gzip -dc .SCRAM/*/ToolCache.db.gz | grep LAPACK_BASE
```
Insert this directory behind `LAPACK =` in the Makefile.
For instance, at MIT, this directory works in the Makefile (set this up once):
```
/cvmfs/cms.cern.ch/slc5_amd64_gcc462/external/lapack/3.3.1
```
at CERN:
```
/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/lapack/3.3.1-cms
```

Build the code with
```
make
```

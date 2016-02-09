#ifndef ForestPFsVector_h
#define ForestPFsVector_h

#include <iostream>
#include <vector>
#include "TBranch.h"

class ForestPFsVector {
public :
   ForestPFsVector(){};
   ~ForestPFsVector(){};

   // Declaration of leaf types
   Int_t                      nPFpart;
   std::vector<int>           *pfId = 0;
   std::vector<float>         *pfPt = 0;
   std::vector<float>         *pfVsPt = 0;
   std::vector<float>         *pfEta = 0;
   std::vector<float>         *pfPhi = 0;
   Float_t                    vn[5][15];
   Float_t                    psin[5][15];
   Float_t                    sumpt[15];

   // List of branches
   TBranch        *b_nPFpart; //!
   TBranch        *b_pfId;    //!
   TBranch        *b_pfPt;    //!
   TBranch        *b_pfVsPt;  //!
   TBranch        *b_pfEta;   //!
   TBranch        *b_pfPhi;   //!
   TBranch        *b_vn;      //!
   TBranch        *b_psin;    //!
   TBranch        *b_sumpt;   //!
};
#endif

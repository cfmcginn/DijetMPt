//=============================================                                 
// Author: Chris McGinn                                                         
//                                                                              
// DiJet Analysis Skim Class (MC)                                                    
//                                                                              
//=============================================  
#ifndef cfmDiJetAnaSkim_h
#define cfmDiJetAnaSkim_h

#include "/net/hisrv0001/home/cfmcginn/DijetMPt/CMSSW_5_3_20/src/DijetInitialSkim/cfmDiJetIniSkim.h"
#include "swapPtCorr.h"
#include "resTrkPtCorr.h"

#include <string>

TTree* trackTreeAna_p = 0;
TTree* jetTreeAna_p = 0;
TTree* genTreeAna_p = 0;

const Int_t nJtAlg = 31;
const Int_t nJtMax = 4;
const Int_t nSumAlg = 31;
const Int_t nPtBins = 7;
const Int_t nRBins = 10;

const Int_t resStartPos = 15;
const Int_t resEndPos = 20;

const Int_t swapStartPos = 21;
const Int_t swapEndPos = 26;

const Float_t ptBins[nPtBins] = {0.00, 0.50, 1.00, 2.00, 4.00, 8.00, 10000.00};

//Track Tree Variables

//Tracks proj. onto Alg (enum ordered above, w/ corrected in back 5), All, Cone, and NotCone

//2nd Dim 0 == 0_1, 1 == 1_2 ... 5 == F

Float_t rAlgJtMult_[2*nSumAlg][2];

Float_t rAlgImbProjA_[2*nSumAlg][nPtBins];
Float_t rAlgImbProjA13_[2*nSumAlg][nPtBins];
Float_t rAlgImbProjAC_[2*nSumAlg][nPtBins];
Float_t rAlgImbProjANC_[2*nSumAlg][nPtBins];

Float_t rAlgImbProjAC0_[2*nSumAlg][nPtBins];
Float_t rAlgImbProjAC1_[2*nSumAlg][nPtBins];
Float_t rAlgImbProjAC2_[2*nSumAlg][nPtBins];
Float_t rAlgImbProjAC3_[2*nSumAlg][nPtBins];

//DelR ProjA

Float_t rAlgImbProjAR_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAR13_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAEta_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAPhi_[2*nSumAlg][nPtBins][nRBins];

Float_t rAlgImbProjAR_Cut_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAR_CutEta_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAR_CutPhi_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAEta_Cut_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAPhi_Cut_[2*nSumAlg][nPtBins][nRBins];

Float_t rAlgImbProjAR_FOR_[2*nSumAlg][nPtBins][nRBins];

//DelR Mult

Float_t rAlgMultAR_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAEta_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAPhi_[2*nSumAlg][nPtBins][nRBins];

Float_t rAlgMultAR_Cut_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAR_CutEta_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAR_CutPhi_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAEta_Cut_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAPhi_Cut_[2*nSumAlg][nPtBins][nRBins];

Float_t rAlgImbRawAR_[2*nSumAlg][nPtBins][nRBins];


//Jet Tree Variables

Int_t run_;
Int_t evt_;
Int_t lumi_;
Int_t hiBin_;

Float_t pthat_;

Float_t hiEvtPlane_;
Float_t psin_;

//Set Bool

//Event Set Bool array, [0] == PuPF, [1] == PuCalo, .etc according to enum

Bool_t eventSet_[nJtAlg];
Float_t centWeight_[nJtAlg];
Float_t centWeight_Merge_[nJtAlg];

Bool_t isQuarkJet_[nJtAlg];
Bool_t isGluonJet_[nJtAlg];

Float_t pthatWeight_;
Float_t leadEtaWeight_;
Float_t subleadEtaWeight_;

Float_t swap12Weight_[nJtAlg];
Float_t swap23Weight_[nJtAlg];

//Jet Set, Array by algorithm, according to enum above, 2nd, 0 = lead, 1 =sublead etc.

Int_t AlgJtN_[nJtAlg];
Float_t AlgJtPt_[nJtAlg][nJtMax];
Float_t AlgJtPhi_[nJtAlg][nJtMax];
Float_t AlgJtEta_[nJtAlg][nJtMax];
Float_t AlgJtTrkMax_[nJtAlg][nJtMax];
Float_t AlgJtRawPt_[nJtAlg][nJtMax];

Float_t AlgJtAvePhi12_[nJtAlg];
Float_t AlgJtAvePhi13_[nJtAlg];
Float_t AlgJtDelPhi12_[nJtAlg];
Float_t AlgJtDelPhi13_[nJtAlg];
Float_t AlgJtDelPhi23_[nJtAlg];

Float_t AlgJtAsymm12_[nJtAlg];
Float_t AlgJtAsymm13_[nJtAlg];
Float_t AlgJtAsymm23_[nJtAlg];
Float_t AlgJtAsymm123_[nJtAlg];

Float_t AlgRefPt_[nJtAlg][nJtMax];
Float_t AlgRefPhi_[nJtAlg][nJtMax];
Float_t AlgRefEta_[nJtAlg][nJtMax];

Float_t AlgRefAvePhi12_[nJtAlg];
Float_t AlgRefAvePhi13_[nJtAlg];
Float_t AlgRefDelPhi12_[nJtAlg];
Float_t AlgRefDelPhi13_[nJtAlg];
Float_t AlgRefDelPhi23_[nJtAlg];
Float_t AlgRefAsymm12_[nJtAlg];
Float_t AlgRefAsymm13_[nJtAlg];
Float_t AlgRefAsymm23_[nJtAlg];
Float_t AlgRefAsymm123_[nJtAlg];

Float_t eventPt_[2];

//Gen Tree Variables

//Gen. proj. onto Jets, ordered by algorithm according to enum, PuPF == [0], PuCalo == [1], etc.

Float_t gAlgJtMult_[nSumAlg][2];

Float_t gAlgImbProjA_[nSumAlg][nPtBins];
Float_t gAlgImbProjAC_[nSumAlg][nPtBins];
Float_t gAlgImbProjANC_[nSumAlg][nPtBins];

//truth delRs
//ProjA delR

Float_t gAlgImbProjAR_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAEta_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAPhi_[nSumAlg][nPtBins][nRBins];

Float_t gAlgImbProjAR_Cut_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAR_CutEta_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAR_CutPhi_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAEta_Cut_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAPhi_Cut_[nSumAlg][nPtBins][nRBins];

Float_t gAlgImbProjAR_FOR_[nSumAlg][nPtBins][nRBins];

//delR Mult

Float_t gAlgMultAR_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAEta_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAPhi_[nSumAlg][nPtBins][nRBins];

Float_t gAlgMultAR_Cut_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAR_CutEta_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAR_CutPhi_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAEta_Cut_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAPhi_Cut_[nSumAlg][nPtBins][nRBins];

Float_t gAlgImbRawAR_[nSumAlg][nPtBins][nRBins];


//Cut var
const Float_t rBounds[nRBins] = {.20, .40, .60, .80, 1.00, 1.20, 1.40, 1.60, 1.80, 100000};

//const Float_t leadJtPtCut[nJtAlg] = {120., 120., 120., 117.448, 120., 123.66, 127.89, 120., 120., 120.};
//const Float_t subLeadJtPtCut[nJtAlg] = {50., 50., 50., 48.9367, 50., 51.5252, 53.2873, 50., 50., 50.};

const Float_t leadJtPtCut = 120.;
const Float_t subLeadJtPtCut = 50.;

const Float_t jtDelPhiCut = 0;
const Float_t jtEtaCut = 2.0; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed                                                 
const std::string algType[nJtAlg] = {"Pu3Calo", "Pu4Calo", "Pu5Calo", "Pu3PF", "Vs2Calo", "Vs3Calo", "Vs4Calo", "Vs5Calo", "Vs3PF", "Pu3CaloFrag", "Vs2CaloFrag", "Vs3CaloFrag", "Vs4CaloFrag", "Vs5CaloFrag", "Vs3PFFrag", "Pu3CaloRes", "Vs2CaloRes", "Vs3CaloRes", "Vs4CaloRes", "Vs5CaloRes", "Vs3PFRes", "Pu3CaloSwap", "Vs2CaloSwap", "Vs3CaloSwap", "Vs4CaloSwap", "Vs5CaloSwap", "Vs3PFSwap", "T2", "T3", "T4", "T5"};


void SetAnaBranches(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  //Track Tree Branches

  std::cout << "Branches Set" << std::endl;

  //Tracks proj. onto Alg, ordered according to enum above, All, Cone, and NotCone

  trackTreeAna_p->Branch("rAlgJtMult", rAlgJtMult_, Form("rAlgJtMult[%d][2]/F", 2*nSumAlg));

  trackTreeAna_p->Branch("rAlgImbProjA", rAlgImbProjA_, Form("rAlgImbProjA[%d][%d]/F", 2*nSumAlg, nPtBins));
  trackTreeAna_p->Branch("rAlgImbProjA13", rAlgImbProjA13_, Form("rAlgImbProjA13[%d][%d]/F", 2*nSumAlg, nPtBins));
  trackTreeAna_p->Branch("rAlgImbProjAC", rAlgImbProjAC_, Form("rAlgImbProjAC[%d][%d]/F", 2*nSumAlg, nPtBins));
  trackTreeAna_p->Branch("rAlgImbProjANC", rAlgImbProjANC_, Form("rAlgImbProjANC[%d][%d]/F", 2*nSumAlg, nPtBins));
  
  trackTreeAna_p->Branch("rAlgImbProjAC0", rAlgImbProjAC0_, Form("rAlgImbProjAC0[%d][%d]/F", 2*nSumAlg, nPtBins));
  trackTreeAna_p->Branch("rAlgImbProjAC1", rAlgImbProjAC1_, Form("rAlgImbProjAC1[%d][%d]/F", 2*nSumAlg, nPtBins));
  trackTreeAna_p->Branch("rAlgImbProjAC2", rAlgImbProjAC2_, Form("rAlgImbProjAC2[%d][%d]/F", 2*nSumAlg, nPtBins));
  trackTreeAna_p->Branch("rAlgImbProjAC3", rAlgImbProjAC3_, Form("rAlgImbProjAC3[%d][%d]/F", 2*nSumAlg, nPtBins));

  //ProjA DelRs
  
  trackTreeAna_p->Branch("rAlgImbProjAR", rAlgImbProjAR_, Form("rAlgImbProjAR[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAR13", rAlgImbProjAR13_, Form("rAlgImbProjAR13[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAEta", rAlgImbProjAEta_, Form("rAlgImbProjAEta[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAPhi", rAlgImbProjAPhi_, Form("rAlgImbProjAPhi[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  trackTreeAna_p->Branch("rAlgImbProjAR_Cut", rAlgImbProjAR_Cut_, Form("rAlgImbProjAR_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAR_CutEta", rAlgImbProjAR_CutEta_, Form("rAlgImbProjAR_CutEta[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAR_CutPhi", rAlgImbProjAR_CutPhi_, Form("rAlgImbProjAR_CutPhi[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAEta_Cut", rAlgImbProjAEta_Cut_, Form("rAlgImbProjAEta_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAPhi_Cut", rAlgImbProjAPhi_Cut_, Form("rAlgImbProjAPhi_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  trackTreeAna_p->Branch("rAlgImbProjAR_FOR", rAlgImbProjAR_FOR_, Form("rAlgImbProjAR_FOR[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  //Mult DelRs

  trackTreeAna_p->Branch("rAlgMultAR", rAlgMultAR_, Form("rAlgMultAR[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAEta", rAlgMultAEta_, Form("rAlgMultAEta[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAPhi", rAlgMultAPhi_, Form("rAlgMultAPhi[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  trackTreeAna_p->Branch("rAlgMultAR_Cut", rAlgMultAR_Cut_, Form("rAlgMultAR_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAR_CutEta", rAlgMultAR_CutEta_, Form("rAlgMultAR_CutEta[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAR_CutPhi", rAlgMultAR_CutPhi_, Form("rAlgMultAR_CutPhi[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAEta_Cut", rAlgMultAEta_Cut_, Form("rAlgMultAEta_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAPhi_Cut", rAlgMultAPhi_Cut_, Form("rAlgMultAPhi_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  trackTreeAna_p->Branch("rAlgImbRawAR", rAlgImbRawAR_, Form("rAlgImbRawAR[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  //Jet Tree Branches

  jetTreeAna_p->Branch("run", &run_, "run/I");
  jetTreeAna_p->Branch("evt", &evt_, "evt/I");
  jetTreeAna_p->Branch("lumi", &lumi_, "lumi/I");

  if(hi){
    jetTreeAna_p->Branch("hiBin", &hiBin_, "hiBin/I");
    jetTreeAna_p->Branch("hiEvtPlane", &hiEvtPlane_, "hiEvtPlane/F");
    jetTreeAna_p->Branch("psin", &psin_, "psin/F");
  }

  jetTreeAna_p->Branch("eventSet", eventSet_, Form("eventSet[%d]/O", nJtAlg));

  if(sType == kHIMC){
    jetTreeAna_p->Branch("centWeight", centWeight_, Form("centWeight[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("centWeight_Merge", centWeight_Merge_, Form("centWeight_Merge[%d]/F", nJtAlg));
  }

  jetTreeAna_p->Branch("AlgJtN", AlgJtN_, Form("AlgJtN[%d]/I", nJtAlg));
  jetTreeAna_p->Branch("AlgJtPt", AlgJtPt_, Form("AlgJtPt[%d][%d]/F", nJtAlg, nJtMax));
  jetTreeAna_p->Branch("AlgJtPhi", AlgJtPhi_, Form("AlgJtPhi[%d][%d]/F", nJtAlg, nJtMax));
  jetTreeAna_p->Branch("AlgJtEta", AlgJtEta_, Form("AlgJtEta[%d][%d]/F", nJtAlg, nJtMax));
  jetTreeAna_p->Branch("AlgJtTrkMax", AlgJtTrkMax_, Form("AlgJtTrkMax[%d][%d]/F", nJtAlg, nJtMax));
  jetTreeAna_p->Branch("AlgJtRawPt", AlgJtRawPt_, Form("AlgJtRawPt[%d][%d]/F", nJtAlg, nJtMax));

  jetTreeAna_p->Branch("AlgJtAvePhi12", AlgJtAvePhi12_, Form("AlgJtAvePhi12[%d]/F", nJtAlg));
  jetTreeAna_p->Branch("AlgJtAvePhi13", AlgJtAvePhi13_, Form("AlgJtAvePhi13[%d]/F", nJtAlg));
  jetTreeAna_p->Branch("AlgJtDelPhi12", AlgJtDelPhi12_, Form("AlgJtDelPhi12[%d]/F", nJtAlg));
  jetTreeAna_p->Branch("AlgJtDelPhi13", AlgJtDelPhi13_, Form("AlgJtDelPhi13[%d]/F", nJtAlg));
  jetTreeAna_p->Branch("AlgJtDelPhi23", AlgJtDelPhi23_, Form("AlgJtDelPhi23[%d]/F", nJtAlg));

  jetTreeAna_p->Branch("AlgJtAsymm12", AlgJtAsymm12_, Form("AlgJtAsymm12[%d]/F", nJtAlg));
  jetTreeAna_p->Branch("AlgJtAsymm13", AlgJtAsymm13_, Form("AlgJtAsymm13[%d]/F", nJtAlg));
  jetTreeAna_p->Branch("AlgJtAsymm23", AlgJtAsymm23_, Form("AlgJtAsymm23[%d]/F", nJtAlg));
  jetTreeAna_p->Branch("AlgJtAsymm123", AlgJtAsymm123_, Form("AlgJtAsymm123[%d]/F", nJtAlg));

  if(montecarlo){
    jetTreeAna_p->Branch("isQuarkJet", isQuarkJet_, Form("isQuarkJet[%d]/O", nJtAlg));
    jetTreeAna_p->Branch("isGluonJet", isGluonJet_, Form("isGluonJet[%d]/O", nJtAlg));
  }

  jetTreeAna_p->Branch("pthatWeight", &pthatWeight_, "pthatWeight/F");
  jetTreeAna_p->Branch("leadEtaWeight", &leadEtaWeight_, "leadEtaWeight/F");
  jetTreeAna_p->Branch("subleadEtaWeight", &subleadEtaWeight_, "subleadEtaWeight/F");

  jetTreeAna_p->Branch("swap12Weight", swap12Weight_, "swap12Weight/F");
  jetTreeAna_p->Branch("swap23Weight", swap23Weight_, "swap23Weight/F");

  if(justJt)
    jetTreeAna_p->Branch("eventPt", eventPt_, "eventPt[2]/F");

  if(montecarlo){
    jetTreeAna_p->Branch("pthat", &pthat_, "pthat/F");

    //refpt for jets immediately above
    jetTreeAna_p->Branch("AlgRefPt", AlgRefPt_, Form("AlgRefPt[%d][%d]/F", nJtAlg, nJtMax));
    jetTreeAna_p->Branch("AlgRefPhi", AlgRefPhi_, Form("AlgRefPhi[%d][%d]/F", nJtAlg, nJtMax));
    jetTreeAna_p->Branch("AlgRefEta", AlgRefEta_, Form("AlgRefEta[%d][%d]/F", nJtAlg, nJtMax));

    jetTreeAna_p->Branch("AlgRefAvePhi12", AlgRefAvePhi12_, Form("AlgRefAvePhi12[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("AlgRefAvePhi13", AlgRefAvePhi13_, Form("AlgRefAvePhi13[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("AlgRefDelPhi12", AlgRefDelPhi12_, Form("AlgRefDelPhi12[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("AlgRefDelPhi13", AlgRefDelPhi13_, Form("AlgRefDelPhi13[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("AlgRefDelPhi23", AlgRefDelPhi23_, Form("AlgRefDelPhi23[%d]/F", nJtAlg));

    jetTreeAna_p->Branch("AlgRefAsymm12", AlgRefAsymm12_, Form("AlgRefAsymm12[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("AlgRefAsymm13", AlgRefAsymm13_, Form("AlgRefAsymm13[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("AlgRefAsymm23", AlgRefAsymm23_, Form("AlgRefAsymm23[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("AlgRefAsymm123", AlgRefAsymm123_, Form("AlgRefAsymm123[%d]/F", nJtAlg));


    //Gen. proj. onto jetAlg, array ordered according to enum

    if(!justJt){
      genTreeAna_p->Branch("gAlgJtMult", gAlgJtMult_, Form("gAlgJtMult[%d][2]/F", nSumAlg));
      genTreeAna_p->Branch("gAlgImbProjA", gAlgImbProjA_, Form("gAlgImbProjA[%d][%d]/F", nSumAlg, nPtBins));
      genTreeAna_p->Branch("gAlgImbProjAC", gAlgImbProjAC_, Form("gAlgImbProjAC[%d][%d]/F", nSumAlg, nPtBins));
      genTreeAna_p->Branch("gAlgImbProjANC", gAlgImbProjANC_, Form("gAlgImbProjANC[%d][%d]/F", nSumAlg, nPtBins));
      
      //Truth Del Rs, Proj
      //Proj A's
      
      genTreeAna_p->Branch("gAlgImbProjAR", gAlgImbProjAR_, Form("gAlgImbProjAR[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAEta", gAlgImbProjAEta_, Form("gAlgImbProjAEta[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAPhi", gAlgImbProjAPhi_, Form("gAlgImbProjAPhi[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));

      genTreeAna_p->Branch("gAlgImbProjAR_Cut", gAlgImbProjAR_Cut_, Form("gAlgImbProjAR_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAR_CutEta", gAlgImbProjAR_CutEta_, Form("gAlgImbProjAR_CutEta[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAR_CutPhi", gAlgImbProjAR_CutPhi_, Form("gAlgImbProjAR_CutPhi[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAEta_Cut", gAlgImbProjAEta_Cut_, Form("gAlgImbProjAEta_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAPhi_Cut", gAlgImbProjAPhi_Cut_, Form("gAlgImbProjAPhi_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));

      genTreeAna_p->Branch("gAlgImbProjAR_FOR", gAlgImbProjAR_FOR_, Form("gAlgImbProjAR_FOR[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));

      //Mult DelRs

      genTreeAna_p->Branch("gAlgMultAR", gAlgMultAR_, Form("gAlgMultAR[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAEta", gAlgMultAEta_, Form("gAlgMultAEta[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAPhi", gAlgMultAPhi_, Form("gAlgMultAPhi[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));

      genTreeAna_p->Branch("gAlgMultAR_Cut", gAlgMultAR_Cut_, Form("gAlgMultAR_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAR_CutEta", gAlgMultAR_CutEta_, Form("gAlgMultAR_CutEta[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAR_CutPhi", gAlgMultAR_CutPhi_, Form("gAlgMultAR_CutPhi[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAEta_Cut", gAlgMultAEta_Cut_, Form("gAlgMultAEta_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAPhi_Cut", gAlgMultAPhi_Cut_, Form("gAlgMultAPhi_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));

      genTreeAna_p->Branch("gAlgImbRawAR", gAlgImbRawAR_, Form("gAlgImbRawAR[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
    }
  }
}


void GetAnaBranches(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);
  
  std::cout << "Get Branches" << std::endl;
  
  //Track Tree Branches
  
  trackTreeAna_p->SetBranchAddress("rAlgJtMult", rAlgJtMult_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjA", rAlgImbProjA_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjA13", rAlgImbProjA13_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAC", rAlgImbProjAC_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjANC", rAlgImbProjANC_);

  trackTreeAna_p->SetBranchAddress("rAlgImbProjAR", rAlgImbProjAR_);    
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAR13", rAlgImbProjAR13_);    
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAEta", rAlgImbProjAEta_);    
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAPhi", rAlgImbProjAPhi_);    

  trackTreeAna_p->SetBranchAddress("rAlgImbProjAR_Cut", rAlgImbProjAR_Cut_);    
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAR_CutEta", rAlgImbProjAR_CutEta_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAR_CutPhi", rAlgImbProjAR_CutPhi_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAEta_Cut", rAlgImbProjAEta_Cut_);    
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAPhi_Cut", rAlgImbProjAPhi_Cut_);    

  trackTreeAna_p->SetBranchAddress("rAlgImbProjAR_FOR", rAlgImbProjAR_FOR_);    

  //Mult DelRs

  trackTreeAna_p->SetBranchAddress("rAlgMultAR", rAlgMultAR_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAEta", rAlgMultAEta_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAPhi", rAlgMultAPhi_);    

  trackTreeAna_p->SetBranchAddress("rAlgMultAR_Cut", rAlgMultAR_Cut_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAR_CutEta", rAlgMultAR_CutEta_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAR_CutPhi", rAlgMultAR_CutPhi_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAEta_Cut", rAlgMultAEta_Cut_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAPhi_Cut", rAlgMultAPhi_Cut_);    

  trackTreeAna_p->SetBranchAddress("rAlgImbRawAR", rAlgImbRawAR_);    
  
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAC0", rAlgImbProjAC0_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAC1", rAlgImbProjAC1_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAC2", rAlgImbProjAC2_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAC3", rAlgImbProjAC3_);
  
  //Jet Tree Branches
  
  jetTreeAna_p->SetBranchAddress("run", &run_);
  jetTreeAna_p->SetBranchAddress("evt", &evt_);
  jetTreeAna_p->SetBranchAddress("lumi", &lumi_);
  
  if(montecarlo)
    jetTreeAna_p->SetBranchAddress("pthat", &pthat_);
  
  if(hi){
    jetTreeAna_p->SetBranchAddress("hiBin", &hiBin_);
    jetTreeAna_p->SetBranchAddress("hiEvtPlane", &hiEvtPlane_);
    jetTreeAna_p->SetBranchAddress("psin", &psin_);
  }
  
  jetTreeAna_p->SetBranchAddress("eventSet", eventSet_);
  
  if(sType == kHIMC){
    jetTreeAna_p->SetBranchAddress("centWeight", centWeight_);
    jetTreeAna_p->SetBranchAddress("centWeight_Merge", centWeight_Merge_);
  }
  
  if(montecarlo){
    jetTreeAna_p->SetBranchAddress("isQuarkJet", isQuarkJet_);
    jetTreeAna_p->SetBranchAddress("isGluonJet", isGluonJet_);
  }     

  jetTreeAna_p->SetBranchAddress("pthatWeight", &pthatWeight_);
  jetTreeAna_p->SetBranchAddress("leadEtaWeight", &leadEtaWeight_);
  jetTreeAna_p->SetBranchAddress("subleadEtaWeight", &subleadEtaWeight_);

  jetTreeAna_p->SetBranchAddress("swap12Weight", swap12Weight_);
  jetTreeAna_p->SetBranchAddress("swap23Weight", swap23Weight_);

  jetTreeAna_p->SetBranchAddress("AlgJtN", AlgJtN_);
  jetTreeAna_p->SetBranchAddress("AlgJtPt", AlgJtPt_);
  jetTreeAna_p->SetBranchAddress("AlgJtPhi", AlgJtPhi_);
  jetTreeAna_p->SetBranchAddress("AlgJtEta", AlgJtEta_);
  jetTreeAna_p->SetBranchAddress("AlgJtTrkMax", AlgJtTrkMax_);
  jetTreeAna_p->SetBranchAddress("AlgJtRawPt", AlgJtRawPt_);
  
  jetTreeAna_p->SetBranchAddress("AlgJtAvePhi12", AlgJtAvePhi12_);
  jetTreeAna_p->SetBranchAddress("AlgJtAvePhi13", AlgJtAvePhi13_);
  jetTreeAna_p->SetBranchAddress("AlgJtDelPhi12", AlgJtDelPhi12_);
  jetTreeAna_p->SetBranchAddress("AlgJtDelPhi13", AlgJtDelPhi13_);
  jetTreeAna_p->SetBranchAddress("AlgJtDelPhi23", AlgJtDelPhi23_);

  jetTreeAna_p->SetBranchAddress("AlgJtAsymm12", AlgJtAsymm12_);
  jetTreeAna_p->SetBranchAddress("AlgJtAsymm13", AlgJtAsymm13_);
  jetTreeAna_p->SetBranchAddress("AlgJtAsymm23", AlgJtAsymm23_);
  jetTreeAna_p->SetBranchAddress("AlgJtAsymm123", AlgJtAsymm123_);  

  if(justJt)
    jetTreeAna_p->SetBranchAddress("eventPt", eventPt_);
  
  if(montecarlo){
    jetTreeAna_p->SetBranchAddress("AlgRefPt", AlgRefPt_);
    jetTreeAna_p->SetBranchAddress("AlgRefPhi", AlgRefPhi_);
    jetTreeAna_p->SetBranchAddress("AlgRefEta", AlgRefEta_);
    
    jetTreeAna_p->SetBranchAddress("AlgRefAvePhi12", AlgRefAvePhi12_);
    jetTreeAna_p->SetBranchAddress("AlgRefAvePhi13", AlgRefAvePhi13_);
    jetTreeAna_p->SetBranchAddress("AlgRefDelPhi12", AlgRefDelPhi12_);
    jetTreeAna_p->SetBranchAddress("AlgRefDelPhi13", AlgRefDelPhi13_);
    jetTreeAna_p->SetBranchAddress("AlgRefDelPhi23", AlgRefDelPhi23_);
    jetTreeAna_p->SetBranchAddress("AlgRefAsymm12", AlgRefAsymm12_);
    jetTreeAna_p->SetBranchAddress("AlgRefAsymm13", AlgRefAsymm13_);
    jetTreeAna_p->SetBranchAddress("AlgRefAsymm23", AlgRefAsymm23_);
    jetTreeAna_p->SetBranchAddress("AlgRefAsymm123", AlgRefAsymm123_);
    
    //Gen Tree Variables
    
    if(!justJt){
      genTreeAna_p->SetBranchAddress("gAlgJtMult", gAlgJtMult_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjA", gAlgImbProjA_);     
      genTreeAna_p->SetBranchAddress("gAlgImbProjAC", gAlgImbProjAC_);     
      genTreeAna_p->SetBranchAddress("gAlgImbProjANC", gAlgImbProjANC_);     
      
      genTreeAna_p->SetBranchAddress("gAlgImbProjAR", gAlgImbProjAR_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAEta", gAlgImbProjAEta_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAPhi", gAlgImbProjAPhi_);

      genTreeAna_p->SetBranchAddress("gAlgImbProjAR_Cut", gAlgImbProjAR_Cut_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAR_CutEta", gAlgImbProjAR_CutEta_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAR_CutPhi", gAlgImbProjAR_CutPhi_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAEta_Cut", gAlgImbProjAEta_Cut_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAPhi_Cut", gAlgImbProjAPhi_Cut_);

      genTreeAna_p->SetBranchAddress("gAlgImbProjAR_FOR", gAlgImbProjAR_FOR_);

      //Mult DelRs

      genTreeAna_p->SetBranchAddress("gAlgMultAR", gAlgMultAR_);
      genTreeAna_p->SetBranchAddress("gAlgMultAEta", gAlgMultAEta_);
      genTreeAna_p->SetBranchAddress("gAlgMultAPhi", gAlgMultAPhi_);

      genTreeAna_p->SetBranchAddress("gAlgMultAR_Cut", gAlgMultAR_Cut_);
      genTreeAna_p->SetBranchAddress("gAlgMultAR_CutEta", gAlgMultAR_CutEta_);
      genTreeAna_p->SetBranchAddress("gAlgMultAR_CutPhi", gAlgMultAR_CutPhi_);
      genTreeAna_p->SetBranchAddress("gAlgMultAEta_Cut", gAlgMultAEta_Cut_);
      genTreeAna_p->SetBranchAddress("gAlgMultAPhi_Cut", gAlgMultAPhi_Cut_);

      genTreeAna_p->SetBranchAddress("gAlgImbRawAR", gAlgImbRawAR_);
    }
  }

  return;
}


void InitDiJetAnaSkim(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  std::cout << "Init DiJet AnaSkim" << std::endl;

  trackTreeAna_p = new TTree("trackTreeAna", "trackTreeAna");
  jetTreeAna_p = new TTree("jetTreeAna", "jetTreeAna");
  if(isMonteCarlo(sType)) genTreeAna_p = new TTree("genTreeAna", "genTreeAna");

  SetAnaBranches(sType, justJt);
}


void CleanupDiJetAnaSkim()
{
  if(trackTreeAna_p != 0) delete trackTreeAna_p;
  if(jetTreeAna_p != 0) delete jetTreeAna_p;
  if(genTreeAna_p != 0) delete genTreeAna_p;
}


void GetDiJetAnaSkim(TFile* anaFile_p, sampleType sType = kHIDATA, Bool_t justJt = false){
  std::cout << "Get DiJet AnaSkim" << std::endl;

  if(!justJt) trackTreeAna_p = (TTree*)anaFile_p->Get("trackTreeAna");
  jetTreeAna_p = (TTree*)anaFile_p->Get("jetTreeAna");
  if(isMonteCarlo(sType) && !justJt) genTreeAna_p = (TTree*)anaFile_p->Get("genTreeAna");

  GetAnaBranches(sType, justJt);
}


Float_t getAbsDphi(Float_t phi1, Float_t phi2, Int_t tag, Int_t tag2 = -1, Int_t tag3 = -1)
{
  return TMath::Abs(getDPHI(phi1, phi2, tag, tag2, tag3));
}


Float_t getFlippedPhi(Float_t inPhi)
{
  Float_t outPhi;

  if(TMath::Abs(inPhi) > TMath::Pi()){
    std::cout << "getFlippedPhi: inPhi is outside accepted range, return -10" << std::endl;
    return -10;
  }
  else if(inPhi > 0)
    outPhi = inPhi - TMath::Pi();
  else
    outPhi = inPhi + TMath::Pi();

  return outPhi;
}


Bool_t sameSign(Float_t num1, Float_t num2){
  if((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0)) return true;

  return false;
}


Float_t getAvePhi(Float_t inLeadPhi, Float_t inSubLeadPhi)
{
  Float_t flipPhi = getFlippedPhi(inSubLeadPhi);
  Float_t avePhi;

  if(sameSign(inLeadPhi, flipPhi) || (TMath::Abs(inLeadPhi) < TMath::Pi()/2 && TMath::Abs(flipPhi) < TMath::Pi()/2))
    avePhi = (flipPhi + inLeadPhi)/2;
  else if(TMath::Abs(inLeadPhi) > TMath::Pi()/2 && TMath::Abs(flipPhi) > TMath::Pi()/2){
    avePhi = (flipPhi + inLeadPhi)/2;
    if(avePhi > 0)
      avePhi = TMath::Pi() - avePhi;
    else
      avePhi = -TMath::Pi() - avePhi;
  }
  else{
    avePhi = 0.;
  }

  return avePhi;
}


Int_t getPtRange(Float_t cutPt)
{
  for(Int_t iter = 0; iter < nPtBins - 1; iter++){
    if(cutPt < ptBins[iter+1]) return iter;
  }

  std::cout << "pt out of range. return -1" << std::endl;
  return -1;
}



void InitJetVar(sampleType sType = kHIDATA)
{
  Bool_t montecarlo = isMonteCarlo(sType);

  for(Int_t initIter = 0; initIter < nJtAlg; initIter++){
    eventSet_[initIter] = false;

    if(sType == kHIMC){
      centWeight_[initIter] = 1;    
      centWeight_Merge_[initIter] = 1;    
    }

    AlgJtN_[initIter] = -999;

    for(Int_t initIter2 = 0; initIter2 < nJtMax; initIter2++){
      AlgJtPt_[initIter][initIter2] = -999;
      AlgJtPhi_[initIter][initIter2] = -999;
      AlgJtEta_[initIter][initIter2] = -999;
      AlgJtTrkMax_[initIter][initIter2] = -999;
      AlgJtRawPt_[initIter][initIter2] = -999;
    }

    AlgJtAvePhi12_[initIter] = -999;
    AlgJtAvePhi13_[initIter] = -999;
    AlgJtDelPhi12_[initIter] = -999;
    AlgJtDelPhi13_[initIter] = -999;
    AlgJtDelPhi23_[initIter] = -999;

    AlgJtAsymm12_[initIter] = -999;
    AlgJtAsymm13_[initIter] = -999;
    AlgJtAsymm23_[initIter] = -999;
    AlgJtAsymm123_[initIter] = -999;

    if(montecarlo){
      isQuarkJet_[initIter] = false;
      isGluonJet_[initIter] = false;

      pthatWeight_ = -999;
      leadEtaWeight_ = -999;
      subleadEtaWeight_ = -999;

      swap12Weight_[initIter] = -999;
      swap23Weight_[initIter] = -999;

      for(Int_t initIter2 = 0; initIter2 < nJtMax; initIter2++){
	AlgRefPt_[initIter][initIter2] = -999;
	AlgRefPhi_[initIter][initIter2] = -999;
	AlgRefEta_[initIter][initIter2] = -999;
      }

      AlgRefAvePhi12_[initIter] = -999;
      AlgRefAvePhi13_[initIter] = -999;
      AlgRefDelPhi12_[initIter] = -999;
      AlgRefDelPhi13_[initIter] = -999;
      AlgRefDelPhi23_[initIter] = -999;

      AlgRefAsymm12_[initIter] = -999;
      AlgRefAsymm13_[initIter] = -999;
      AlgRefAsymm23_[initIter] = -999;
      AlgRefAsymm123_[initIter] = -999;
    }
  }
}

void InitProjPerp(sampleType sType = kHIDATA)
{
  //Tracks proj. onto Alg, ordered according to enum above, corr in the back 5, All, Cone, and NotCone

  for(Int_t initIter = 0; initIter < 2*nSumAlg; initIter++){
    for(Int_t initIter2 = 0; initIter2 < nPtBins; initIter2++){
      if(initIter2 < 2) rAlgJtMult_[initIter][initIter2] = 0;

      rAlgImbProjA_[initIter][initIter2] = 0;
      rAlgImbProjA13_[initIter][initIter2] = 0;
      rAlgImbProjAC_[initIter][initIter2] = 0;
      rAlgImbProjANC_[initIter][initIter2] = 0;

      rAlgImbProjAC0_[initIter][initIter2] = 0;
      rAlgImbProjAC1_[initIter][initIter2] = 0;
      rAlgImbProjAC2_[initIter][initIter2] = 0;
      rAlgImbProjAC3_[initIter][initIter2] = 0;

      for(Int_t initIter3 = 0; initIter3 < nRBins; initIter3++){
	rAlgImbProjAR_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAR13_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAEta_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAPhi_[initIter][initIter2][initIter3] = 0;

	rAlgImbProjAR_Cut_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAR_CutEta_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAR_CutPhi_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAEta_Cut_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAPhi_Cut_[initIter][initIter2][initIter3] = 0;

	rAlgImbProjAR_FOR_[initIter][initIter2][initIter3] = 0;

	//Init Mults

	rAlgMultAR_[initIter][initIter2][initIter3] = 0;
	rAlgMultAEta_[initIter][initIter2][initIter3] = 0;
	rAlgMultAPhi_[initIter][initIter2][initIter3] = 0;

	rAlgMultAR_Cut_[initIter][initIter2][initIter3] = 0;
	rAlgMultAR_CutEta_[initIter][initIter2][initIter3] = 0;
	rAlgMultAR_CutPhi_[initIter][initIter2][initIter3] = 0;
	rAlgMultAEta_Cut_[initIter][initIter2][initIter3] = 0;
	rAlgMultAPhi_Cut_[initIter][initIter2][initIter3] = 0;

	rAlgImbRawAR_[initIter][initIter2][initIter3] = 0;
      }
    }
  }

  if(isMonteCarlo(sType)){
    //Gen. proj. onto Truth
    for(Int_t initIter = 0; initIter < nSumAlg; initIter++){
      for(Int_t initIter2 = 0; initIter2 < nPtBins; initIter2++){
	if(initIter2 < 2) gAlgJtMult_[initIter][initIter2] = 0;

	gAlgImbProjA_[initIter][initIter2] = 0;
	gAlgImbProjAC_[initIter][initIter2] = 0;
	gAlgImbProjANC_[initIter][initIter2] = 0;

	for(Int_t initIter3 = 0; initIter3 < nRBins; initIter3++){
	  gAlgImbProjAR_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAEta_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAPhi_[initIter][initIter2][initIter3] = 0;

	  gAlgImbProjAR_Cut_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAR_CutEta_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAR_CutPhi_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAEta_Cut_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAPhi_Cut_[initIter][initIter2][initIter3] = 0;

	  gAlgImbProjAR_FOR_[initIter][initIter2][initIter3] = 0;

	  //Init Mults

	  gAlgMultAR_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAEta_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAPhi_[initIter][initIter2][initIter3] = 0;

	  gAlgMultAR_Cut_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAR_CutEta_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAR_CutPhi_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAEta_Cut_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAPhi_Cut_[initIter][initIter2][initIter3] = 0;

	  gAlgImbRawAR_[initIter][initIter2][initIter3] = 0;
	}
      }
    }    
  }
}


Float_t getSwapComb(Float_t imbProjAR, Float_t imbProjAR13, Float_t swap12Weight, Float_t swap23Weight, Float_t jtDelPhi13, Float_t jtAsymm12, Float_t jtAsymm23)
{
  Float_t combImbProjAR = imbProjAR;
  if(jtAsymm12 > 0.05 && (swap23Weight == 1 || jtDelPhi13 < 5*TMath::Pi()/6)) combImbProjAR /= (2*swap12Weight - 1);
  else if(swap23Weight > 0 && jtAsymm23 > 0.05 && jtDelPhi13 > 5*TMath::Pi()/6) combImbProjAR = swap23Weight*combImbProjAR + (1-swap23Weight)*imbProjAR13;

  return combImbProjAR;
}

void getJtVar(Int_t nJt, Float_t jtPt[], Float_t jtPhi[], Float_t jtEta[], Float_t trkMax[], Float_t rawPt[], Float_t refPt[], Float_t refPhi[], Float_t refEta[], Int_t refPart[], Int_t algNum, Bool_t montecarlo = false, Bool_t truth = false)
{
  AlgJtN_[algNum] = nJt;
  if(nJt == 0 || nJt == 1 || jtPt[0] < leadJtPtCut || jtPt[1] < subLeadJtPtCut) return;

  eventSet_[algNum] = true;

  Int_t iterMax = nJt;
  if(nJt > nJtMax) iterMax = nJtMax;

  for(Int_t leadIter = 0; leadIter < iterMax; leadIter++){
    AlgJtPt_[algNum][leadIter] = jtPt[leadIter];
    AlgJtPhi_[algNum][leadIter] = jtPhi[leadIter];
    AlgJtEta_[algNum][leadIter] = jtEta[leadIter];

    if(!truth){
      AlgJtTrkMax_[algNum][leadIter] = trkMax[leadIter];
      AlgJtRawPt_[algNum][leadIter] = rawPt[leadIter];

      if(montecarlo){
        AlgRefPt_[algNum][leadIter] = refPt[leadIter];
        AlgRefPhi_[algNum][leadIter] = refPhi[leadIter];
        AlgRefEta_[algNum][leadIter] = refEta[leadIter];
      }
    }
  }

  if(montecarlo){
    if(TMath::Abs(refPart[0]) < 9) isQuarkJet_[algNum] = true;
    else if(TMath::Abs(refPart[0]) == 21) isGluonJet_[algNum] = true;
  }

  AlgJtAvePhi12_[algNum] = getAvePhi(AlgJtPhi_[algNum][0], AlgJtPhi_[algNum][1]);
  AlgJtDelPhi12_[algNum] = getAbsDphi(AlgJtPhi_[algNum][0], AlgJtPhi_[algNum][1], 5);
  AlgJtAsymm12_[algNum] = (AlgJtPt_[algNum][0] - AlgJtPt_[algNum][1])/(AlgJtPt_[algNum][0] + AlgJtPt_[algNum][1]);
  
  if(AlgJtPt_[algNum][2] > 30 && TMath::Abs(AlgJtEta_[algNum][2]) < 2.0){
    AlgJtDelPhi13_[algNum] = getAbsDphi(AlgJtPhi_[algNum][0], AlgJtPhi_[algNum][2], 6);
    AlgJtAvePhi13_[algNum] = getAvePhi(AlgJtPhi_[algNum][0], AlgJtPhi_[algNum][2]);
    AlgJtDelPhi23_[algNum] = getAbsDphi(AlgJtPhi_[algNum][1], AlgJtPhi_[algNum][2], 7);

    AlgJtAsymm13_[algNum] = (AlgJtPt_[algNum][0] - AlgJtPt_[algNum][2])/(AlgJtPt_[algNum][0] + AlgJtPt_[algNum][2]);
    AlgJtAsymm23_[algNum] = (AlgJtPt_[algNum][1] - AlgJtPt_[algNum][2])/(AlgJtPt_[algNum][1] + AlgJtPt_[algNum][2]);
    Float_t temp23PtX = AlgJtPt_[algNum][1]*cos(AlgJtPhi_[algNum][1]) + AlgJtPt_[algNum][2]*cos(AlgJtPhi_[algNum][2]);
    Float_t temp23PtY = AlgJtPt_[algNum][1]*sin(AlgJtPhi_[algNum][1]) + AlgJtPt_[algNum][2]*sin(AlgJtPhi_[algNum][2]);
    Float_t temp23Pt = TMath::Sqrt(temp23PtX*temp23PtX + temp23PtY*temp23PtY);;
    AlgJtAsymm123_[algNum] = (AlgJtPt_[algNum][0] - temp23Pt)/(AlgJtPt_[algNum][0] + temp23Pt);
  }
  else AlgJtAsymm123_[algNum] = AlgJtAsymm12_[algNum];
    
  if(montecarlo && !truth){
    if(AlgRefPhi_[algNum][0] > -9 && AlgRefPhi_[algNum][1] > -9){
      AlgRefAvePhi12_[algNum] = getAvePhi(AlgRefPhi_[algNum][0], AlgRefPhi_[algNum][1]);
      AlgRefDelPhi12_[algNum] = getAbsDphi(AlgRefPhi_[algNum][0], AlgRefPhi_[algNum][1], 8);
      AlgRefAsymm12_[algNum] = (AlgRefPt_[algNum][0] - AlgRefPt_[algNum][1])/(AlgRefPt_[algNum][0] + AlgRefPt_[algNum][1]);
      if(AlgRefPhi_[algNum][2] > -9 && nJt > 2){
	AlgRefDelPhi13_[algNum] = getAbsDphi(AlgRefPhi_[algNum][0], AlgRefPhi_[algNum][2], 9);
	AlgRefAvePhi13_[algNum] = getAvePhi(AlgRefPhi_[algNum][0], AlgRefPhi_[algNum][2]);
	AlgRefDelPhi23_[algNum] = getAbsDphi(AlgRefPhi_[algNum][1], AlgRefPhi_[algNum][2], 10);
    
	AlgRefAsymm13_[algNum] = (AlgRefPt_[algNum][0] - AlgRefPt_[algNum][2])/(AlgRefPt_[algNum][0] + AlgRefPt_[algNum][2]);
	AlgRefAsymm23_[algNum] = (AlgRefPt_[algNum][1] - AlgRefPt_[algNum][2])/(AlgRefPt_[algNum][1] + AlgRefPt_[algNum][2]);


	Float_t temp23PtX = AlgRefPt_[algNum][1]*cos(AlgRefPhi_[algNum][1]) + AlgRefPt_[algNum][2]*cos(AlgRefPhi_[algNum][2]);
	Float_t temp23PtY = AlgRefPt_[algNum][1]*sin(AlgRefPhi_[algNum][1]) + AlgRefPt_[algNum][2]*sin(AlgRefPhi_[algNum][2]);
	Float_t temp23Pt = TMath::Sqrt(temp23PtX*temp23PtX + temp23PtY*temp23PtY);;
	AlgRefAsymm123_[algNum] = (AlgRefPt_[algNum][0] - temp23Pt)/(AlgRefPt_[algNum][0] + temp23Pt);
      }
      else AlgRefAsymm123_[algNum] = AlgRefAsymm12_[algNum];

    }
    else{
      AlgRefAvePhi12_[algNum] = -999;
      AlgRefDelPhi12_[algNum] = -999;
      AlgRefDelPhi13_[algNum] = -999;
      AlgRefAvePhi13_[algNum] = -999;
      AlgRefDelPhi23_[algNum] = -999;
      AlgRefAsymm12_[algNum] = -999;
      AlgRefAsymm13_[algNum] = -999;
      AlgRefAsymm23_[algNum] = -999;
      AlgRefAsymm123_[algNum] = -999;
    }
  }

  return;
}

void GetTrkProjPerp(Int_t jtAlg, Int_t jtAlgCorr, Float_t rawPt, Float_t corrPt, Float_t phi, Float_t eta, Bool_t doCorrPt2, sampleType sType)
{
  Int_t multVal;

  if(jtAlg == jtAlgCorr) multVal = 1;
  else multVal = corrPt/rawPt;

  Int_t signVal = 1;
  if(cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2)) < 0) signVal = -1;

  if(getAbsDphi(AlgJtAvePhi12_[jtAlg], phi, 2) < TMath::Pi()/2) rAlgJtMult_[jtAlgCorr][0] += multVal;
  else rAlgJtMult_[jtAlgCorr][1] += multVal;

  Int_t ptIter = getPtRange(rawPt);
  Float_t tempLeadDelR = getDR(eta, phi, AlgJtEta_[jtAlg][0], AlgJtPhi_[jtAlg][0], 2);
  Float_t tempSubLeadDelR = getDR(eta, phi, AlgJtEta_[jtAlg][1], AlgJtPhi_[jtAlg][1], 2);
  Float_t tempThirdDelR = 99.0;
  if(AlgJtPt_[jtAlg][2] > 30) tempThirdDelR = getDR(eta, phi, AlgJtEta_[jtAlg][2], AlgJtPhi_[jtAlg][2]);
  Float_t tempLeadDelEta = TMath::Abs(eta - AlgJtEta_[jtAlg][0]);
  Float_t tempSubLeadDelEta = TMath::Abs(eta - AlgJtEta_[jtAlg][1]);
  Float_t tempLeadDelPhi = TMath::Abs(getDPHI(phi, AlgJtPhi_[jtAlg][0], 2));
  Float_t tempSubLeadDelPhi = TMath::Abs(getDPHI(phi, AlgJtPhi_[jtAlg][1], 2));

  Float_t corr2 = 1.0;
  if(doCorrPt2){
    Int_t pos = -1;
    if(jtAlg == 15 || jtAlg == 21) pos = 0;
    else if(jtAlg == 16 || jtAlg == 22) pos = 1;
    else if(jtAlg == 17 || jtAlg == 23) pos = 2;
    else if(jtAlg == 18 || jtAlg == 24) pos = 3;
    else if(jtAlg == 19 || jtAlg == 25) pos = 4;

    if(sType == kPPMC || sType == kPPDATA) pos-=1;

    corr2 = GetRESTrkPtCorr(pos, tempLeadDelR, tempSubLeadDelR, AlgJtAsymm12_[jtAlg], rawPt, sType);
  }

  if(ptIter != 0) rAlgImbProjA_[jtAlgCorr][nPtBins - 1] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2))*corr2;
  rAlgImbProjA_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2))*corr2;

  if(AlgJtPt_[jtAlg][2] > 30){
    if(ptIter != 0) rAlgImbProjA13_[jtAlgCorr][nPtBins - 1] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi13_[jtAlg], 2));
    rAlgImbProjA13_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi13_[jtAlg], 2));
  }


  if(tempLeadDelR > 0 && tempSubLeadDelR > 0){
    if(tempLeadDelR < 0.8 || tempSubLeadDelR < 0.8){
      if(ptIter != 0) rAlgImbProjAC_[jtAlgCorr][nPtBins - 1] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
      rAlgImbProjAC_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
    }
    else{
      if(ptIter != 0) rAlgImbProjANC_[jtAlgCorr][nPtBins - 1] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
      rAlgImbProjANC_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
    }

    if(tempLeadDelR < 0.5 || tempSubLeadDelR < 0.5){
      if(ptIter != 0) rAlgImbProjAC0_[jtAlgCorr][nPtBins - 1] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
      rAlgImbProjAC0_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
    }
    else if(tempLeadDelR < 1.0 || tempSubLeadDelR < 1.0){
      if(ptIter != 0) rAlgImbProjAC1_[jtAlgCorr][nPtBins - 1] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
      rAlgImbProjAC1_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
    }
    else if(tempLeadDelR < 1.5 || tempSubLeadDelR < 1.5){
      if(ptIter != 0) rAlgImbProjAC2_[jtAlgCorr][nPtBins - 1] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
      rAlgImbProjAC2_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
    }
    else{
      if(ptIter != 0) rAlgImbProjAC3_[jtAlgCorr][nPtBins - 1] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
      rAlgImbProjAC3_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){

	if(ptIter != 0) rAlgImbProjAR_[jtAlgCorr][nPtBins - 1][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2))*corr2;	
	rAlgImbProjAR_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2))*corr2;

	if(ptIter != 0) rAlgMultAR_[jtAlgCorr][nPtBins - 1][rIter] -= signVal*multVal;	
	rAlgMultAR_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	if(ptIter != 0) rAlgImbRawAR_[jtAlgCorr][nPtBins - 1][rIter] += TMath::Abs(corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2)));	
	rAlgImbRawAR_[jtAlgCorr][ptIter][rIter] += TMath::Abs(corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2)));
	break;
      }
    }

    if(AlgJtPt_[jtAlg][2] > 30){
      for(Int_t rIter = 0; rIter < 10; rIter++){
	if(tempLeadDelR < rBounds[rIter] || tempThirdDelR < rBounds[rIter]){
	  if(ptIter != 0) rAlgImbProjAR13_[jtAlgCorr][nPtBins - 1][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi13_[jtAlg], 2));
	  rAlgImbProjAR13_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi13_[jtAlg], 2));

	  break;
	}
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	if(ptIter != 0) rAlgImbProjAR_Cut_[jtAlgCorr][nPtBins - 1][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));	
	rAlgImbProjAR_Cut_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));

	if(ptIter != 0) rAlgMultAR_Cut_[jtAlgCorr][nPtBins - 1][rIter] -= signVal*multVal;	
	rAlgMultAR_Cut_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	if(ptIter != 0) rAlgImbProjAR_CutEta_[jtAlgCorr][nPtBins - 1][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));	
	rAlgImbProjAR_CutEta_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));

	if(ptIter != 0) rAlgMultAR_CutEta_[jtAlgCorr][nPtBins - 1][rIter] -= signVal*multVal;	
	rAlgMultAR_CutEta_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	if(ptIter != 0) rAlgImbProjAR_CutPhi_[jtAlgCorr][nPtBins - 1][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));	
	rAlgImbProjAR_CutPhi_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));

	if(ptIter != 0) rAlgMultAR_CutPhi_[jtAlgCorr][nPtBins - 1][rIter] -= signVal*multVal;	
	rAlgMultAR_CutPhi_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi < rBounds[rIter] || tempSubLeadDelPhi < rBounds[rIter]){
	if(ptIter != 0) rAlgImbProjAPhi_[jtAlgCorr][nPtBins - 1][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));	
	rAlgImbProjAPhi_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));

	if(ptIter != 0) rAlgMultAPhi_[jtAlgCorr][nPtBins - 1][rIter] -= signVal*multVal;
	rAlgMultAPhi_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;
      if(tempLeadDelPhi < rBounds[rIter] || tempSubLeadDelPhi < rBounds[rIter]){
	if(ptIter != 0) rAlgImbProjAPhi_Cut_[jtAlgCorr][nPtBins - 1][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));	
	rAlgImbProjAPhi_Cut_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));

	if(ptIter != 0) rAlgMultAPhi_Cut_[jtAlgCorr][nPtBins - 1][rIter] -= signVal*multVal;	
	rAlgMultAPhi_Cut_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta < rBounds[rIter] || tempSubLeadDelEta < rBounds[rIter]){
	if(ptIter != 0) rAlgImbProjAEta_[jtAlgCorr][nPtBins - 1][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));	
	rAlgImbProjAEta_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));

	if(ptIter != 0) rAlgMultAEta_[jtAlgCorr][nPtBins - 1][rIter] -= signVal*multVal;	
	rAlgMultAEta_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;
      if(tempLeadDelEta < rBounds[rIter] || tempSubLeadDelEta < rBounds[rIter]){
	if(ptIter != 0) rAlgImbProjAEta_Cut_[jtAlgCorr][nPtBins - 1][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));	
	rAlgImbProjAEta_Cut_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));

	if(ptIter != 0) rAlgMultAEta_Cut_[jtAlgCorr][nPtBins - 1][rIter] -= signVal*multVal;	
	rAlgMultAEta_Cut_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    if((AlgJtEta_[jtAlg][0] < 0.6 && AlgJtEta_[jtAlg][1] < 0.6)){
      if(eta > TMath::Min(AlgJtEta_[jtAlg][0], AlgJtEta_[jtAlg][1]) - 0.3){
	for(Int_t rIter = 0; rIter < 10; rIter++){
	  if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	    if(ptIter != 0) rAlgImbProjAR_FOR_[jtAlgCorr][nPtBins - 1][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
	    rAlgImbProjAR_FOR_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
	    break;
	  }
	}
      }
    }
    else if((AlgJtEta_[jtAlg][0] > -0.6 && AlgJtEta_[jtAlg][1] > -0.6)){
      if(eta < TMath::Max(AlgJtEta_[jtAlg][0], AlgJtEta_[jtAlg][1]) + 0.3){
	for(Int_t rIter = 0; rIter < 10; rIter++){
	  if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	    if(ptIter != 0) rAlgImbProjAR_FOR_[jtAlgCorr][nPtBins - 1][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
	    rAlgImbProjAR_FOR_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 2));
	    break;
	  }
	}
      }
    }
  
  }

  return;
}


void GetMixProjPerp(Float_t evPt[2], Int_t nConst, Float_t constPt[], Float_t constPhi[], Float_t constEta[], Float_t constCorr[])
{
  for(Int_t constIter = 0; constIter < nConst; constIter++){
    evPt[0] -= constPt[constIter]*constCorr[constIter]*cos(constPhi[constIter]);
    evPt[1] -= constPt[constIter]*constCorr[constIter]*sin(constPhi[constIter]);
    GetTrkProjPerp(1, 1, constPt[constIter], constPt[constIter]*constCorr[constIter], constPhi[constIter], constEta[constIter], true, kHIDATA);
    GetTrkProjPerp(1, 4, constPt[constIter], constPt[constIter]*constCorr[constIter], constPhi[constIter], constEta[constIter], true, kHIDATA);
  }

  return;
}


void GetGenProjPerp(Int_t jtAlg, Float_t pt, Float_t phi, Float_t eta)
{
  if(TMath::IsNaN(pt)) return;

  Int_t signVal = 1;
  if(cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3)) < 0) signVal = -1;

  if(getAbsDphi(AlgJtAvePhi12_[jtAlg], phi, 3) < TMath::Pi()/2) gAlgJtMult_[jtAlg][0] += 1;
  else gAlgJtMult_[jtAlg][1] += 1;

  Int_t ptIter = getPtRange(pt);
  Float_t tempLeadDelR = getDR(eta, phi, AlgJtEta_[jtAlg][0], AlgJtPhi_[jtAlg][0], 3);
  Float_t tempSubLeadDelR = getDR(eta, phi, AlgJtEta_[jtAlg][1], AlgJtPhi_[jtAlg][1], 3);
  Float_t tempLeadDelEta = TMath::Abs(eta - AlgJtEta_[jtAlg][0]);
  Float_t tempSubLeadDelEta = TMath::Abs(eta - AlgJtEta_[jtAlg][1]);
  Float_t tempLeadDelPhi = TMath::Abs(getDPHI(phi, AlgJtPhi_[jtAlg][0], 3));
  Float_t tempSubLeadDelPhi = TMath::Abs(getDPHI(phi, AlgJtPhi_[jtAlg][1], 3));

  if(ptIter != 0) gAlgImbProjA_[jtAlg][nPtBins - 1] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
  gAlgImbProjA_[jtAlg][ptIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));

  if(tempLeadDelR > 0 && tempSubLeadDelR > 0){
    if(tempLeadDelR < 0.8 || tempSubLeadDelR < 0.8){
      if(ptIter != 0) gAlgImbProjAC_[jtAlg][nPtBins - 1] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
      gAlgImbProjAC_[jtAlg][ptIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
    }
    else{
      if(ptIter != 0) gAlgImbProjANC_[jtAlg][nPtBins - 1] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
      gAlgImbProjANC_[jtAlg][ptIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	if(ptIter != 0) gAlgImbProjAR_[jtAlg][nPtBins - 1][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	gAlgImbProjAR_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));

	if(ptIter != 0) gAlgMultAR_[jtAlg][nPtBins - 1][rIter] -= signVal;
	gAlgMultAR_[jtAlg][ptIter][rIter] -= signVal;

	if(ptIter != 0) gAlgImbRawAR_[jtAlg][nPtBins - 1][rIter] += TMath::Abs(pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3)));
	gAlgImbRawAR_[jtAlg][ptIter][rIter] += TMath::Abs(pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3)));
	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	if(ptIter != 0) gAlgImbProjAR_Cut_[jtAlg][nPtBins - 1][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	gAlgImbProjAR_Cut_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));

	if(ptIter != 0) gAlgMultAR_Cut_[jtAlg][nPtBins - 1][rIter] -= signVal;
	gAlgMultAR_Cut_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	if(ptIter != 0) gAlgImbProjAR_CutEta_[jtAlg][nPtBins - 1][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	gAlgImbProjAR_CutEta_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));

	if(ptIter != 0) gAlgMultAR_CutEta_[jtAlg][nPtBins - 1][rIter] -= signVal;
	gAlgMultAR_CutEta_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	if(ptIter != 0) gAlgImbProjAR_CutPhi_[jtAlg][nPtBins - 1][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	gAlgImbProjAR_CutPhi_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));

	if(ptIter != 0) gAlgMultAR_CutPhi_[jtAlg][nPtBins - 1][rIter] -= signVal;
	gAlgMultAR_CutPhi_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }


    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta < rBounds[rIter] || tempSubLeadDelEta < rBounds[rIter]){
	if(ptIter != 0) gAlgImbProjAEta_[jtAlg][nPtBins - 1][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	gAlgImbProjAEta_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));

	if(ptIter != 0) gAlgMultAEta_[jtAlg][nPtBins - 1][rIter] -= signVal;
	gAlgMultAEta_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }


    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;

      if(tempLeadDelEta < rBounds[rIter] || tempSubLeadDelEta < rBounds[rIter]){
	if(ptIter != 0) gAlgImbProjAEta_Cut_[jtAlg][nPtBins - 1][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	gAlgImbProjAEta_Cut_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));

	if(ptIter != 0) gAlgMultAEta_Cut_[jtAlg][nPtBins - 1][rIter] -= signVal;
	gAlgMultAEta_Cut_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi < rBounds[rIter] || tempSubLeadDelPhi < rBounds[rIter]){
	if(ptIter != 0) gAlgImbProjAPhi_[jtAlg][nPtBins - 1][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	gAlgImbProjAPhi_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));

	if(ptIter != 0) gAlgMultAPhi_[jtAlg][nPtBins - 1][rIter] -= signVal;
	gAlgMultAPhi_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }


    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;

      if(tempLeadDelPhi < rBounds[rIter] || tempSubLeadDelPhi < rBounds[rIter]){
	if(ptIter != 0) gAlgImbProjAPhi_Cut_[jtAlg][nPtBins - 1][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	gAlgImbProjAPhi_Cut_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));

	if(ptIter != 0) gAlgMultAPhi_Cut_[jtAlg][nPtBins - 1][rIter] -= signVal;
	gAlgMultAPhi_Cut_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }


    if((AlgJtEta_[jtAlg][0] < 0.6 && AlgJtEta_[jtAlg][1] < 0.6)){
      if(eta > TMath::Min(AlgJtEta_[jtAlg][0], AlgJtEta_[jtAlg][1]) - 0.3){
	for(Int_t rIter = 0; rIter < 10; rIter++){
	  if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	    if(ptIter != 0) gAlgImbProjAR_FOR_[jtAlg][nPtBins - 1][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	    gAlgImbProjAR_FOR_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	    break;
	  }
	}
      }
    }
    else if((AlgJtEta_[jtAlg][0] > -0.6 && AlgJtEta_[jtAlg][1] > -0.6)){
      if(eta < TMath::Max(AlgJtEta_[jtAlg][0], AlgJtEta_[jtAlg][1]) + 0.3){
	for(Int_t rIter = 0; rIter < 10; rIter++){
	  if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	    if(ptIter != 0) gAlgImbProjAR_FOR_[jtAlg][nPtBins - 1][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	    gAlgImbProjAR_FOR_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi12_[jtAlg], 3));
	    break;
	  }
	}
      }
    }


  }

  return;
}


#endif

//=============================================
// Author: Chris McGinn / Doga Gulhan
// 
// DiJet Analysis Skim Class (MC)
//
//=============================================

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "cfmDiJetAnaSkim.h"
#include "stdlib.h"
#include <fstream>
#include "TComplex.h"

const Int_t pthatCuts_PYTH_HITrk[10] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 10000000};
const Float_t pthatWeights_PYTH_HITrk[9] = {.551019, .034814, .00242254, .000304825, .0000426788, .00000492814, .000000879673, .00000017353, .0000000292439};

const Int_t pthatCuts_PYTH_PPTrk[7] = {15, 30, 50, 80, 120, 170, 1000000};
const Float_t pthatWeights_PYTH_PPTrk[6] = {.161482, .00749461, .000752396, .0000837038, .0000101988, .00000175206};

const Int_t pthatCuts_PYTH_HYD[11] = {15, 30, 50, 80, 100, 120, 170, 220, 280, 370, 10000000};
const Float_t pthatWeights_PYTH_HYD[10] = {.611066, .0374106, .00232016, .00014917, .0000822379, .0000142819, .00000296162, .00000102099, .000000522123, .000000232907};

const Int_t jtAlgRBin[nJtAlg] = {1, 2, 3, 1, 0, 1, 2, 3, 1, 1, 0, 1, 2, 3, 1, 1, 0, 1, 2, 3, 1, 1, 0, 1, 2, 3, 1, 0, 1, 2, 3};

int makeDiJetAnaSkim(std::string fList = "", sampleType sType = kHIDATA, Int_t num = 0, Bool_t justJt = false, Bool_t isHITrk = false)
{
  //Define MC or Data
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  std::cout << sType << std::endl;
  std::cout << montecarlo << std::endl;

  std::string buffer;
  std::vector<std::string> listOfFiles;
  int nLines = 0;
  ifstream inFile(fList.data());

  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
    return 1;
  }
  else{
    while(true){
      inFile >> buffer;
      if(inFile.eof()) break;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  std::cout << "FileList Loaded" << std::endl;

  std::cout << "FileJob: " << listOfFiles[num] << std::endl;

  TFile* iniSkim_p = new TFile(listOfFiles[num].data(), "READ");

  GetDiJetIniSkim(iniSkim_p, sType, justJt);

  std::cout << "IniSkim Loaded" << std::endl;

  //Setup correction tables

  InitFactCorrFiles(sType, isHITrk);
  InitFactCorrHists(sType);
  InitSWAPCorrFiles(sType);
  InitSWAPCorrFits(sType);
  InitRESTrkCorrFiles(sType);
  InitRESTrkCorrHists(sType);
  
  TFile *histWeightFile_p = new TFile("histWeightFile.root", "READ");
  TH1F *hist_DataOverMC_p[nSumAlg - 10];

  
  if(sType == kHIMC){
    for(Int_t algIter = 0; algIter < nSumAlg-10; algIter++){
      hist_DataOverMC_p[algIter] = (TH1F*)histWeightFile_p->Get(Form("ak%s_dataHiBin_h", algType[algIter].c_str()));
      std::cout << Form("ak%s_dataHiBin_h", algType[algIter].c_str()) << std::endl;
    }
  }
    
  std::string outName = listOfFiles[num];
  const std::string cutString = "/";
  const std::string iniString = "Ini";
  std::size_t strIndex = 0;

  std::cout << "Cull string" << std::endl;

  while(true){
    strIndex = outName.find(cutString);

    if(strIndex == std::string::npos) break;

    outName.replace(0, strIndex + 1, "");
  }

  std::cout << "Replace string" << std::endl;

  strIndex = outName.find(iniString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, iniString.length(), "Ana"); 
  }

  std::cout << "Output name: " << outName.c_str() << std::endl;

  TFile *outFile = new TFile(outName.c_str(), "RECREATE");

  InitDiJetAnaSkim(sType, justJt);

  Long64_t nentries = jetTreeIni_p->GetEntries();

  std::cout << nentries << std::endl;

  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    jetTreeIni_p->GetEntry(jentry);

    if(!justJt){
      trackTreeIni_p->GetEntry(jentry);
      
      if(montecarlo)
	genTreeIni_p->GetEntry(jentry);
    }      

    if(jentry%1000 == 0) std::cout << jentry << std::endl;

    InitJetVar(sType);

    getJtVar(nPu3Calo_, Pu3CaloPt_, Pu3CaloPhi_, Pu3CaloEta_, Pu3CaloTrkMax_, Pu3CaloRawPt_, Pu3CaloRefPt_, Pu3CaloRefPhi_, Pu3CaloRefEta_, Pu3CaloRefPart_, 0, montecarlo, false);
    getJtVar(nPu4Calo_, Pu4CaloPt_, Pu4CaloPhi_, Pu4CaloEta_, Pu4CaloTrkMax_, Pu4CaloRawPt_, Pu4CaloRefPt_, Pu4CaloRefPhi_, Pu4CaloRefEta_, Pu4CaloRefPart_, 1, montecarlo, false);
    getJtVar(nPu5Calo_, Pu5CaloPt_, Pu5CaloPhi_, Pu5CaloEta_, Pu5CaloTrkMax_, Pu5CaloRawPt_, Pu5CaloRefPt_, Pu5CaloRefPhi_, Pu5CaloRefEta_, Pu5CaloRefPart_, 2, montecarlo, false);
    getJtVar(nPu3PF_, Pu3PFPt_, Pu3PFPhi_, Pu3PFEta_, Pu3PFTrkMax_, Pu3PFRawPt_, Pu3PFRefPt_, Pu3PFRefPhi_, Pu3PFRefEta_, Pu3PFRefPart_, 3, montecarlo, false);
    getJtVar(nVs2Calo_, Vs2CaloPt_, Vs2CaloPhi_, Vs2CaloEta_, Vs2CaloTrkMax_, Vs2CaloRawPt_, Vs2CaloRefPt_, Vs2CaloRefPhi_, Vs2CaloRefEta_, Vs2CaloRefPart_, 4, montecarlo, false);
    getJtVar(nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_, Vs3CaloTrkMax_, Vs3CaloRawPt_, Vs3CaloRefPt_, Vs3CaloRefPhi_, Vs3CaloRefEta_, Vs3CaloRefPart_, 5, montecarlo, false);
    getJtVar(nVs4Calo_, Vs4CaloPt_, Vs4CaloPhi_, Vs4CaloEta_, Vs4CaloTrkMax_, Vs4CaloRawPt_, Vs4CaloRefPt_, Vs4CaloRefPhi_, Vs4CaloRefEta_, Vs4CaloRefPart_, 6, montecarlo, false);
    getJtVar(nVs5Calo_, Vs5CaloPt_, Vs5CaloPhi_, Vs5CaloEta_, Vs5CaloTrkMax_, Vs5CaloRawPt_, Vs5CaloRefPt_, Vs5CaloRefPhi_, Vs5CaloRefEta_, Vs5CaloRefPart_, 7, montecarlo, false);
    getJtVar(nVs3PF_, Vs3PFPt_, Vs3PFPhi_, Vs3PFEta_, Vs3PFTrkMax_, Vs3PFRawPt_, Vs3PFRefPt_, Vs3PFRefPhi_, Vs3PFRefEta_, Vs3PFRefPart_, 8, montecarlo, false);


    getJtVar(nPu3CaloFrag_, Pu3CaloFragPt_, Pu3CaloFragPhi_, Pu3CaloFragEta_, Pu3CaloFragTrkMax_, Pu3CaloFragRawPt_, Pu3CaloFragRefPt_, Pu3CaloFragRefPhi_, Pu3CaloFragRefEta_, Pu3CaloFragRefPart_, 9, montecarlo, false);
    getJtVar(nVs2CaloFrag_, Vs2CaloFragPt_, Vs2CaloFragPhi_, Vs2CaloFragEta_, Vs2CaloFragTrkMax_, Vs2CaloFragRawPt_, Vs2CaloFragRefPt_, Vs2CaloFragRefPhi_, Vs2CaloFragRefEta_, Vs2CaloFragRefPart_, 10, montecarlo, false);
    getJtVar(nVs3CaloFrag_, Vs3CaloFragPt_, Vs3CaloFragPhi_, Vs3CaloFragEta_, Vs3CaloFragTrkMax_, Vs3CaloFragRawPt_, Vs3CaloFragRefPt_, Vs3CaloFragRefPhi_, Vs3CaloFragRefEta_, Vs3CaloFragRefPart_, 11, montecarlo, false);
    getJtVar(nVs4CaloFrag_, Vs4CaloFragPt_, Vs4CaloFragPhi_, Vs4CaloFragEta_, Vs4CaloFragTrkMax_, Vs4CaloFragRawPt_, Vs4CaloFragRefPt_, Vs4CaloFragRefPhi_, Vs4CaloFragRefEta_, Vs4CaloFragRefPart_, 12, montecarlo, false);
    getJtVar(nVs5CaloFrag_, Vs5CaloFragPt_, Vs5CaloFragPhi_, Vs5CaloFragEta_, Vs5CaloFragTrkMax_, Vs5CaloFragRawPt_, Vs5CaloFragRefPt_, Vs5CaloFragRefPhi_, Vs5CaloFragRefEta_, Vs5CaloFragRefPart_, 13, montecarlo, false);
    getJtVar(nVs3PFFrag_, Vs3PFFragPt_, Vs3PFFragPhi_, Vs3PFFragEta_, Vs3PFFragTrkMax_, Vs3PFFragRawPt_, Vs3PFFragRefPt_, Vs3PFFragRefPhi_, Vs3PFFragRefEta_, Vs3PFFragRefPart_, 14, montecarlo, false);


    getJtVar(nPu3CaloRes_, Pu3CaloResPt_, Pu3CaloResPhi_, Pu3CaloResEta_, Pu3CaloResTrkMax_, Pu3CaloResRawPt_, Pu3CaloResRefPt_, Pu3CaloResRefPhi_, Pu3CaloResRefEta_, Pu3CaloResRefPart_, 15, montecarlo, false);
    getJtVar(nVs2CaloRes_, Vs2CaloResPt_, Vs2CaloResPhi_, Vs2CaloResEta_, Vs2CaloResTrkMax_, Vs2CaloResRawPt_, Vs2CaloResRefPt_, Vs2CaloResRefPhi_, Vs2CaloResRefEta_, Vs2CaloResRefPart_, 16, montecarlo, false);
    getJtVar(nVs3CaloRes_, Vs3CaloResPt_, Vs3CaloResPhi_, Vs3CaloResEta_, Vs3CaloResTrkMax_, Vs3CaloResRawPt_, Vs3CaloResRefPt_, Vs3CaloResRefPhi_, Vs3CaloResRefEta_, Vs3CaloResRefPart_, 17, montecarlo, false);
    getJtVar(nVs4CaloRes_, Vs4CaloResPt_, Vs4CaloResPhi_, Vs4CaloResEta_, Vs4CaloResTrkMax_, Vs4CaloResRawPt_, Vs4CaloResRefPt_, Vs4CaloResRefPhi_, Vs4CaloResRefEta_, Vs4CaloResRefPart_, 18, montecarlo, false);
    getJtVar(nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloResTrkMax_, Vs5CaloResRawPt_, Vs5CaloResRefPt_, Vs5CaloResRefPhi_, Vs5CaloResRefEta_, Vs5CaloResRefPart_, 19, montecarlo, false);
    getJtVar(nVs3PFRes_, Vs3PFResPt_, Vs3PFResPhi_, Vs3PFResEta_, Vs3PFResTrkMax_, Vs3PFResRawPt_, Vs3PFResRefPt_, Vs3PFResRefPhi_, Vs3PFResRefEta_, Vs3PFResRefPart_, 20, montecarlo, false);

    getJtVar(nPu3CaloRes_, Pu3CaloResPt_, Pu3CaloResPhi_, Pu3CaloResEta_, Pu3CaloResTrkMax_, Pu3CaloResRawPt_, Pu3CaloResRefPt_, Pu3CaloResRefPhi_, Pu3CaloResRefEta_, Pu3CaloResRefPart_, 21, montecarlo, false);
    getJtVar(nVs2CaloRes_, Vs2CaloResPt_, Vs2CaloResPhi_, Vs2CaloResEta_, Vs2CaloResTrkMax_, Vs2CaloResRawPt_, Vs2CaloResRefPt_, Vs2CaloResRefPhi_, Vs2CaloResRefEta_, Vs2CaloResRefPart_, 22, montecarlo, false);
    getJtVar(nVs3CaloRes_, Vs3CaloResPt_, Vs3CaloResPhi_, Vs3CaloResEta_, Vs3CaloResTrkMax_, Vs3CaloResRawPt_, Vs3CaloResRefPt_, Vs3CaloResRefPhi_, Vs3CaloResRefEta_, Vs3CaloResRefPart_, 23, montecarlo, false);
    getJtVar(nVs4CaloRes_, Vs4CaloResPt_, Vs4CaloResPhi_, Vs4CaloResEta_, Vs4CaloResTrkMax_, Vs4CaloResRawPt_, Vs4CaloResRefPt_, Vs4CaloResRefPhi_, Vs4CaloResRefEta_, Vs4CaloResRefPart_, 24, montecarlo, false);
    getJtVar(nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloResTrkMax_, Vs5CaloResRawPt_, Vs5CaloResRefPt_, Vs5CaloResRefPhi_, Vs5CaloResRefEta_, Vs5CaloResRefPart_, 25, montecarlo, false);
    getJtVar(nVs3PFRes_, Vs3PFResPt_, Vs3PFResPhi_, Vs3PFResEta_, Vs3PFResTrkMax_, Vs3PFResRawPt_, Vs3PFResRefPt_, Vs3PFResRefPhi_, Vs3PFResRefEta_, Vs3PFResRefPart_, 26, montecarlo, false);

    Float_t dummyArray2[nT2_];
    Float_t dummyArray3[nT3_];
    Float_t dummyArray4[nT4_];
    Float_t dummyArray5[nT5_];
    getJtVar(nT2_, T2Pt_, T2Phi_, T2Eta_, dummyArray2, dummyArray2, dummyArray2, dummyArray2, dummyArray2, T2Part_, 27, montecarlo, true);
    getJtVar(nT3_, T3Pt_, T3Phi_, T3Eta_, dummyArray3, dummyArray3, dummyArray3, dummyArray3, dummyArray3, T3Part_, 28, montecarlo, true);
    getJtVar(nT4_, T4Pt_, T4Phi_, T4Eta_, dummyArray4, dummyArray4, dummyArray4, dummyArray4, dummyArray4, T4Part_, 29, montecarlo, true);
    getJtVar(nT5_, T5Pt_, T5Phi_, T5Eta_, dummyArray5, dummyArray5, dummyArray5, dummyArray5, dummyArray5, T5Part_, 30, montecarlo, true);

    Bool_t isEventPass = false;
    for(Int_t evtIter = 0; evtIter < nJtAlg; evtIter++){
      if(eventSet_[evtIter]){
	isEventPass = true;
	break;
      }
    }

    if(!isEventPass){
      std::cout << "No event pass after IniSkim; Potential bug" << std::endl;
      continue;
    }

    run_ = runIni_;
    evt_ = evtIni_;
    lumi_ = lumiIni_;

    if(montecarlo)
      pthat_ = pthatIni_;
    
    if(hi){
      hiBin_ = hiBinIni_;
      hiEvtPlane_ = hiEvtPlaneIni_;
      psin_ = psinIni_;
    }
    
    if(montecarlo){
      if(isHITrk){
	for(Int_t hatIter = 0; hatIter < 9; hatIter++){
	  if(pthat_ >= pthatCuts_PYTH_HITrk[hatIter] && pthat_ < pthatCuts_PYTH_HITrk[hatIter + 1]){
	    pthatWeight_ = pthatWeights_PYTH_HITrk[hatIter];
	    break;
	  }
	}
      }
      else if(hi){
	for(Int_t hatIter = 0; hatIter < 10; hatIter++){
          if(pthat_ >= pthatCuts_PYTH_HYD[hatIter] && pthat_ < pthatCuts_PYTH_HYD[hatIter + 1]){
            pthatWeight_ = pthatWeights_PYTH_HYD[hatIter];
            break;
          }
        }
      }
      else{
        for(Int_t hatIter = 0; hatIter < 6; hatIter++){
          if(pthat_ >= pthatCuts_PYTH_PPTrk[hatIter] && pthat_ < pthatCuts_PYTH_PPTrk[hatIter + 1]){
            pthatWeight_ = pthatWeights_PYTH_PPTrk[hatIter];
            break;
          }
        }
      }
    }
    
    for(Int_t algIter = 0; algIter < nJtAlg; algIter++){
      swap12Weight_[algIter] = GetJtSWAP12Prob(sType, jtAlgRBin[algIter], hiBin_, AlgJtAsymm12_[algIter]);
      if(AlgJtPt_[algIter][2] > 50 && TMath::Abs(AlgJtEta_[algIter][2]) < 2.0) swap23Weight_[algIter] = GetJtSWAP23Prob(sType, jtAlgRBin[algIter], hiBin_, AlgJtAsymm23_[algIter]);
    }

    if(sType == kHIMC){
      for(Int_t algIter = 0; algIter < nSumAlg - 10; algIter++){
	centWeight_[algIter] = hist_DataOverMC_p[algIter]->GetBinContent(hist_DataOverMC_p[algIter]->FindBin(hiBin_));
      }
      centWeight_[nSumAlg - 10] = centWeight_[nSumAlg - 16];
      centWeight_[nSumAlg - 9] = centWeight_[nSumAlg - 15];
      centWeight_[nSumAlg - 8] = centWeight_[nSumAlg - 14];
      centWeight_[nSumAlg - 7] = centWeight_[nSumAlg - 13];
      centWeight_[nSumAlg - 6] = centWeight_[nSumAlg - 12];
      centWeight_[nSumAlg - 5] = centWeight_[nSumAlg - 11];


      centWeight_[nSumAlg - 4] = centWeight_[nSumAlg - 9];
      centWeight_[nSumAlg - 3] = centWeight_[nSumAlg - 8];
      centWeight_[nSumAlg - 2] = centWeight_[nSumAlg - 7];
      centWeight_[nSumAlg - 1] = centWeight_[nSumAlg - 6];
    }
    

   //Iterate over tracks

    InitProjPerp(sType);

    //Switch below to iterated OR EDIT HERE

    isEventPass = false;
    for(Int_t evtIter = 0; evtIter < nSumAlg; evtIter++){
      if(eventSet_[evtIter]){
        isEventPass = true;
        break;
      }
    }


    if(isEventPass && !justJt){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
        
	//Grab proj. Pt Spectra For Tracks in each Event Subset
	
	for(Int_t jtIter = 0; jtIter < nSumAlg; jtIter++){
	  Bool_t doCorrPt2 = false;
	  if(hi && jtIter >= swapStartPos && jtIter <= swapEndPos - 1) doCorrPt2 = true;
	  if(!hi && jtIter >= swapStartPos +1 && jtIter <= swapEndPos - 1) doCorrPt2 = true;

	  if(eventSet_[jtIter]) GetTrkProjPerp(jtIter, jtIter, trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], doCorrPt2, sType);
	}
      }

      if(hi) InitPosArrPbPb(hiBin_);

      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	Int_t ptPos = getPtBin(trkPt_[trkEntry], sType);

	Float_t tempRMin[nSumAlg];
	Float_t tempFact[nSumAlg];
	Float_t tempCorr[nSumAlg];

	for(Int_t tempIter = 0; tempIter < nSumAlg; tempIter++){
	  tempRMin[tempIter] = 199.;
	  tempFact[tempIter] = 0.;
	  tempCorr[tempIter] = 0.;
	}
	
	if(eventSet_[Pu3Calo]) tempRMin[Pu3Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Pu4Calo]) tempRMin[Pu4Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Pu5Calo]) tempRMin[Pu5Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Pu3PF]) tempRMin[Pu3PF] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs2Calo]) tempRMin[Vs2Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs3Calo]) tempRMin[Vs3Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs4Calo]) tempRMin[Vs4Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs5Calo]) tempRMin[Vs5Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs3PF]) tempRMin[Vs3PF] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);

	if(eventSet_[Pu3CaloFrag]) tempRMin[Pu3CaloFrag] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs2CaloFrag]) tempRMin[Vs2CaloFrag] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs3CaloFrag]) tempRMin[Vs3CaloFrag] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs4CaloFrag]) tempRMin[Vs4CaloFrag] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs5CaloFrag]) tempRMin[Vs5CaloFrag] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs3PFFrag]) tempRMin[Vs3PFFrag] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);

	if(eventSet_[Pu3CaloRes]) tempRMin[Pu3CaloRes] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs2CaloRes]) tempRMin[Vs2CaloRes] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs3CaloRes]) tempRMin[Vs3CaloRes] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs4CaloRes]) tempRMin[Vs4CaloRes] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs5CaloRes]) tempRMin[Vs5CaloRes] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs3PFRes]) tempRMin[Vs3PFRes] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	
	if(eventSet_[Pu3CaloSwap]) tempRMin[Pu3CaloSwap] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs2CaloSwap]) tempRMin[Vs2CaloSwap] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs3CaloSwap]) tempRMin[Vs3CaloSwap] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs4CaloSwap]) tempRMin[Vs4CaloSwap] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs5CaloSwap]) tempRMin[Vs5CaloSwap] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[Vs3PFSwap]) tempRMin[Vs3PFSwap] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);

	if(eventSet_[T2]) tempRMin[T2] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[T3]) tempRMin[T3] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[T4]) tempRMin[T4] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);
	if(eventSet_[T5]) tempRMin[T5] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs5CaloRes_, Vs5CaloResPt_, Vs5CaloResPhi_, Vs5CaloResEta_, Vs5CaloRefPt_, montecarlo);


	
	if(eventSet_[Pu3Calo]){
	  tempFact[Pu3Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Pu3Calo], sType);
	  tempCorr[Pu3Calo] = trkPt_[trkEntry]*tempFact[Pu3Calo];
	}
	if(eventSet_[Pu4Calo]){
	  tempFact[Pu4Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Pu4Calo], sType);
	  tempCorr[Pu4Calo] = trkPt_[trkEntry]*tempFact[Pu4Calo];
	}
	if(eventSet_[Pu5Calo]){
	  tempFact[Pu5Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Pu5Calo], sType);
	  tempCorr[Pu5Calo] = trkPt_[trkEntry]*tempFact[Pu5Calo];
	}
	if(eventSet_[Pu3PF]){
	  tempFact[Pu3PF] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Pu3PF], sType);
	  tempCorr[Pu3PF] = trkPt_[trkEntry]*tempFact[Pu3PF];
	}
	if(eventSet_[Vs2Calo]){
	  tempFact[Vs2Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs2Calo], sType);
	  tempCorr[Vs2Calo] = trkPt_[trkEntry]*tempFact[Vs2Calo];
	}
	if(eventSet_[Vs3Calo]){
	  tempFact[Vs3Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs3Calo], sType);
	  tempCorr[Vs3Calo] = trkPt_[trkEntry]*tempFact[Vs3Calo];
	}
	if(eventSet_[Vs4Calo]){
	  tempFact[Vs4Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs4Calo], sType);
	  tempCorr[Vs4Calo] = trkPt_[trkEntry]*tempFact[Vs4Calo];
	}
	if(eventSet_[Vs5Calo]){
	  tempFact[Vs5Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs5Calo], sType);
	  tempCorr[Vs5Calo] = trkPt_[trkEntry]*tempFact[Vs5Calo];
	}
	if(eventSet_[Vs3PF]){
	  tempFact[Vs3PF] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs3PF], sType);
	  tempCorr[Vs3PF] = trkPt_[trkEntry]*tempFact[Vs3PF];
	}



	if(eventSet_[Pu3CaloFrag]){
	  tempFact[Pu3CaloFrag] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Pu3CaloFrag], sType);
	  tempCorr[Pu3CaloFrag] = trkPt_[trkEntry]*tempFact[Pu3CaloFrag];
	}
	if(eventSet_[Vs2CaloFrag]){
	  tempFact[Vs2CaloFrag] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs2CaloFrag], sType);
	  tempCorr[Vs2CaloFrag] = trkPt_[trkEntry]*tempFact[Vs2CaloFrag];
	}
	if(eventSet_[Vs3CaloFrag]){
	  tempFact[Vs3CaloFrag] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs3CaloFrag], sType);
	  tempCorr[Vs3CaloFrag] = trkPt_[trkEntry]*tempFact[Vs3CaloFrag];
	}
	if(eventSet_[Vs4CaloFrag]){
	  tempFact[Vs4CaloFrag] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs4CaloFrag], sType);
	  tempCorr[Vs4CaloFrag] = trkPt_[trkEntry]*tempFact[Vs4CaloFrag];
	}
	if(eventSet_[Vs5CaloFrag]){
	  tempFact[Vs5CaloFrag] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs5CaloFrag], sType);
	  tempCorr[Vs5CaloFrag] = trkPt_[trkEntry]*tempFact[Vs5CaloFrag];
	}
	if(eventSet_[Vs3PFFrag]){
	  tempFact[Vs3PFFrag] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs3PFFrag], sType);
	  tempCorr[Vs3PFFrag] = trkPt_[trkEntry]*tempFact[Vs3PFFrag];
	}


	if(eventSet_[Pu3CaloRes]){
	  tempFact[Pu3CaloRes] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Pu3CaloRes], sType);
	  tempCorr[Pu3CaloRes] = trkPt_[trkEntry]*tempFact[Pu3CaloRes];
	}
	if(eventSet_[Vs2CaloRes]){
	  tempFact[Vs2CaloRes] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs2CaloRes], sType);
	  tempCorr[Vs2CaloRes] = trkPt_[trkEntry]*tempFact[Vs2CaloRes];
	}
	if(eventSet_[Vs3CaloRes]){
	  tempFact[Vs3CaloRes] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs3CaloRes], sType);
	  tempCorr[Vs3CaloRes] = trkPt_[trkEntry]*tempFact[Vs3CaloRes];
	}
	if(eventSet_[Vs4CaloRes]){
	  tempFact[Vs4CaloRes] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs4CaloRes], sType);
	  tempCorr[Vs4CaloRes] = trkPt_[trkEntry]*tempFact[Vs4CaloRes];
	}
	if(eventSet_[Vs5CaloRes]){
	  tempFact[Vs5CaloRes] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs5CaloRes], sType);
	  tempCorr[Vs5CaloRes] = trkPt_[trkEntry]*tempFact[Vs5CaloRes];
	}
	if(eventSet_[Vs3PFRes]){
	  tempFact[Vs3PFRes] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs3PFRes], sType);
	  tempCorr[Vs3PFRes] = trkPt_[trkEntry]*tempFact[Vs3PFRes];
	}
	

	if(eventSet_[Pu3CaloSwap]){
	  tempFact[Pu3CaloSwap] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Pu3CaloSwap], sType);
	  tempCorr[Pu3CaloSwap] = trkPt_[trkEntry]*tempFact[Pu3CaloSwap];
	}
	if(eventSet_[Vs2CaloSwap]){
	  tempFact[Vs2CaloSwap] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs2CaloSwap], sType);
	  tempCorr[Vs2CaloSwap] = trkPt_[trkEntry]*tempFact[Vs2CaloSwap];
	}
	
	if(eventSet_[Vs3CaloSwap]){
	  tempFact[Vs3CaloSwap] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs3CaloSwap], sType);
	  tempCorr[Vs3CaloSwap] = trkPt_[trkEntry]*tempFact[Vs3CaloSwap];
	}
	
	if(eventSet_[Vs4CaloSwap]){
	  tempFact[Vs4CaloSwap] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs4CaloSwap], sType);
	  tempCorr[Vs4CaloSwap] = trkPt_[trkEntry]*tempFact[Vs4CaloSwap];
	}
	if(eventSet_[Vs5CaloSwap]){
	  tempFact[Vs5CaloSwap] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs5CaloSwap], sType);
	  tempCorr[Vs5CaloSwap] = trkPt_[trkEntry]*tempFact[Vs5CaloSwap];
	}
	if(eventSet_[Vs3PFSwap]){
	  tempFact[Vs3PFSwap] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs3PFSwap], sType);
	  tempCorr[Vs3PFSwap] = trkPt_[trkEntry]*tempFact[Vs3PFSwap];
	}

	if(montecarlo){
	  if(eventSet_[T2]){
	    tempFact[T2] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T2], sType);
	    tempCorr[T2] = trkPt_[trkEntry]*tempFact[T2];
	  }

	  if(eventSet_[T3]){
	    tempFact[T3] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T3], sType);
	    tempCorr[T3] = trkPt_[trkEntry]*tempFact[T3];
	  }

	  if(eventSet_[T4]){
	    tempFact[T4] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T4], sType);
	    tempCorr[T4] = trkPt_[trkEntry]*tempFact[T4];
	  }

	  if(eventSet_[T5]){
	    tempFact[T5] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T5], sType);
	    tempCorr[T5] = trkPt_[trkEntry]*tempFact[T5];
	  }
	}
	
	for(Int_t jtIter = 0; jtIter < nSumAlg; jtIter++){
	  Bool_t doCorrPt2 = false;
          if(hi){
	    if(jtIter >= resStartPos && jtIter <= resEndPos - 1) doCorrPt2 = true;
	    if(jtIter >= swapStartPos && jtIter <= swapEndPos - 1) doCorrPt2 = true;
	  }
	  if(!hi){
	    if(jtIter >= resStartPos + 1 && jtIter <= resEndPos - 1) doCorrPt2 = true;
	    if(jtIter >= swapStartPos + 1 && jtIter <= swapEndPos - 1) doCorrPt2 = true;
	  }

	  if(eventSet_[jtIter]) GetTrkProjPerp(jtIter, jtIter+nSumAlg, trkPt_[trkEntry], tempCorr[jtIter], trkPhi_[trkEntry], trkEta_[trkEntry], doCorrPt2, sType);
	}	
      }

      
      for(Int_t jtIter = swapStartPos+1; jtIter <= swapEndPos-1; jtIter++){
	Int_t pos1 = jtIter + nSumAlg;

	for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
	  rAlgImbProjA_[jtIter][ptIter] = getSwapComb(rAlgImbProjA_[jtIter][ptIter], rAlgImbProjA13_[jtIter][ptIter], swap12Weight_[jtIter], swap23Weight_[jtIter], AlgJtDelPhi13_[jtIter], AlgJtAsymm12_[jtIter], AlgJtAsymm23_[jtIter]);
	  rAlgImbProjA_[pos1][ptIter] = getSwapComb(rAlgImbProjA_[pos1][ptIter], rAlgImbProjA13_[pos1][ptIter], swap12Weight_[jtIter], swap23Weight_[jtIter], AlgJtDelPhi13_[jtIter], AlgJtAsymm12_[jtIter], AlgJtAsymm23_[jtIter]);

	  for(Int_t rIter = 0; rIter < nRBins; rIter++){
	    rAlgImbProjAR_[jtIter][ptIter][rIter] = getSwapComb(rAlgImbProjAR_[jtIter][ptIter][rIter], rAlgImbProjAR13_[jtIter][ptIter][rIter], swap12Weight_[jtIter], swap23Weight_[jtIter], AlgJtDelPhi13_[jtIter], AlgJtAsymm12_[jtIter], AlgJtAsymm23_[jtIter]);
	    rAlgImbProjAR_[pos1][ptIter][rIter] = getSwapComb(rAlgImbProjAR_[pos1][ptIter][rIter], rAlgImbProjAR13_[pos1][ptIter][rIter], swap12Weight_[jtIter], swap23Weight_[jtIter], AlgJtDelPhi13_[jtIter], AlgJtAsymm12_[jtIter], AlgJtAsymm23_[jtIter]);
	  }
	}
      }
      

     if(montecarlo){
	//Iterate over Truth
	for(Int_t genEntry = 0; genEntry < nGen_; genEntry++){
	  for(Int_t jtIter = 0; jtIter < nSumAlg; jtIter++){
	    if(eventSet_[jtIter]) GetGenProjPerp(jtIter, genPt_[genEntry], genPhi_[genEntry], genEta_[genEntry]);
	  }
	}
      }

    }
    jetTreeAna_p->Fill();

    if(!justJt){
      trackTreeAna_p->Fill();
    
      if(montecarlo) genTreeAna_p->Fill();
    }
  }
  
  outFile->cd();
  jetTreeAna_p->Write("", TObject::kOverwrite);

  if(!justJt){
    trackTreeAna_p->Write("", TObject::kOverwrite);

    if(montecarlo) genTreeAna_p->Write("", TObject::kOverwrite);
  }

   
  histWeightFile_p->Close();
  delete histWeightFile_p;
  

  CleanupDiJetAnaSkim();
  outFile->Close();
  delete outFile;

  iniSkim_p->Close();
  delete iniSkim_p;

  printf("Done.\n");
  return(0);
}


int main(int argc, char *argv[])
{
  if(argc != 6){
    std::cout << "Usage: sortForest <inputFile> <sampleType> <#> <justJtBool> <isHITrk>" << std::endl;
    return 1;
  }

  int rStatus = -1;

  rStatus = makeDiJetAnaSkim(argv[1], sampleType(atoi(argv[2])), atoi(argv[3]), Bool_t(atoi(argv[4])), Bool_t(atoi(argv[5])));

  return rStatus;
}

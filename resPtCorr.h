//=============================================
// Author: Chris McGinn
//
// JEC Res. Corr. Class (MC)
//
//=============================================     

#include "TFile.h"
#include "TF1.h"
#include "sType.h"

const Int_t nFileRESPbPb = 6;
const Int_t nFileRESPP = 4;

TFile* RES_File0_PbPb_p[nFileRESPbPb];
TF1* RES_Corr0_010_f[nFileRESPbPb];
TF1* RES_Corr0_1030_f[nFileRESPbPb];
TF1* RES_Corr0_3050_f[nFileRESPbPb];
TF1* RES_Corr0_50100_f[nFileRESPbPb];

TFile* RES_File1_PbPb_p[nFileRESPbPb];
TF1* RES_Corr1_010_f[nFileRESPbPb];
TF1* RES_Corr1_1030_f[nFileRESPbPb];
TF1* RES_Corr1_3050_f[nFileRESPbPb];
TF1* RES_Corr1_50100_f[nFileRESPbPb];

TFile* RES_File0_PP_p[nFileRESPP];
TF1* RES_Corr0_PP_f[nFileRESPP];

TFile* RES_File1_PP_p[nFileRESPP];
TF1* RES_Corr1_PP_f[nFileRESPP];

TFile* RES_File2_PP_p[nFileRESPP];
TF1* RES_Corr2_PP_f[nFileRESPP];

void InitRESCorrFiles(sampleType sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    RES_File0_PbPb_p[0] = new TFile("residualcorr0_akPu3Calo.root");
    RES_File0_PbPb_p[1] = new TFile("residualcorr0_akVs2Calo.root");
    RES_File0_PbPb_p[2] = new TFile("residualcorr0_akVs3Calo.root");
    RES_File0_PbPb_p[3] = new TFile("residualcorr0_akVs4Calo.root");
    RES_File0_PbPb_p[4] = new TFile("residualcorr0_akVs5Calo.root");
    RES_File0_PbPb_p[5] = new TFile("residualcorr0_akVs3PF.root");

    RES_File1_PbPb_p[0] = new TFile("residualcorr1_akPu3Calo.root");
    RES_File1_PbPb_p[1] = new TFile("residualcorr1_akVs2Calo.root");
    RES_File1_PbPb_p[2] = new TFile("residualcorr1_akVs3Calo.root");
    RES_File1_PbPb_p[3] = new TFile("residualcorr1_akVs4Calo.root");
    RES_File1_PbPb_p[4] = new TFile("residualcorr1_akVs5Calo.root");
    RES_File1_PbPb_p[5] = new TFile("residualcorr1_akVs3PF.root");
  }
  else if(sType == kPPDATA || sType == kPPMC){
    RES_File0_PP_p[0] = new TFile("residualcorr0_ak2Calo.root");
    RES_File0_PP_p[1] = new TFile("residualcorr0_ak3Calo.root");
    RES_File0_PP_p[2] = new TFile("residualcorr0_ak4Calo.root");
    RES_File0_PP_p[3] = new TFile("residualcorr0_ak5Calo.root");

    RES_File1_PP_p[0] = new TFile("residualcorr1_ak2Calo.root");
    RES_File1_PP_p[1] = new TFile("residualcorr1_ak3Calo.root");
    RES_File1_PP_p[2] = new TFile("residualcorr1_ak4Calo.root");
    RES_File1_PP_p[3] = new TFile("residualcorr1_ak5Calo.root");

    RES_File2_PP_p[0] = new TFile("residualcorr2_ak2Calo.root");
    RES_File2_PP_p[1] = new TFile("residualcorr2_ak3Calo.root");
    RES_File2_PP_p[2] = new TFile("residualcorr2_ak4Calo.root");
    RES_File2_PP_p[3] = new TFile("residualcorr2_ak5Calo.root");
  }  

  return;
}

void InitRESCorrFits(sampleType sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    for(Int_t iter = 0; iter < nFileRESPbPb; iter++){
      RES_Corr0_010_f[iter] = (TF1*)RES_File0_PbPb_p[iter]->Get("fit0");
      RES_Corr0_1030_f[iter] = (TF1*)RES_File0_PbPb_p[iter]->Get("fit1");
      RES_Corr0_3050_f[iter] = (TF1*)RES_File0_PbPb_p[iter]->Get("fit2");
      RES_Corr0_50100_f[iter] = (TF1*)RES_File0_PbPb_p[iter]->Get("fit3");

      RES_Corr1_010_f[iter] = (TF1*)RES_File1_PbPb_p[iter]->Get("fit0");
      RES_Corr1_1030_f[iter] = (TF1*)RES_File1_PbPb_p[iter]->Get("fit1");
      RES_Corr1_3050_f[iter] = (TF1*)RES_File1_PbPb_p[iter]->Get("fit2");
      RES_Corr1_50100_f[iter] = (TF1*)RES_File1_PbPb_p[iter]->Get("fit3");
    }
  }
  else if(sType == kPPDATA || sType == kPPMC){
    for(Int_t iter = 0; iter < nFileRESPP; iter++){
      RES_Corr0_PP_f[iter] = (TF1*)RES_File0_PP_p[iter]->Get("fit0");
      RES_Corr1_PP_f[iter] = (TF1*)RES_File1_PP_p[iter]->Get("fit0");
      RES_Corr2_PP_f[iter] = (TF1*)RES_File2_PP_p[iter]->Get("fit0");
    }
  }

  return;
}


Float_t GetJtRESCorrPt(sampleType sType, Int_t algBin, Int_t hiBin, Float_t jtPt)
{
  Float_t corrPt = 1.0;
  const Int_t hiBinCut[4] = {19, 59, 99, 199};

  Float_t jtPtCheck = jtPt;
  if(jtPt < 25) jtPtCheck = 25;

  if(jtPt > 15 && jtPt < 700){
    if(sType == kHIDATA || sType == kHIMC){

      if(hiBin <= hiBinCut[0]) corrPt = RES_Corr0_010_f[algBin]->Eval(jtPtCheck);
      else if(hiBin <= hiBinCut[1]) corrPt = RES_Corr0_1030_f[algBin]->Eval(jtPtCheck);
      else if(hiBin <= hiBinCut[2]) corrPt = RES_Corr0_3050_f[algBin]->Eval(jtPtCheck);
      else corrPt = RES_Corr0_50100_f[algBin]->Eval(jtPtCheck);

      jtPtCheck = jtPt/corrPt;
      if(jtPtCheck < 25) jtPtCheck = 25;

      if(hiBin <= hiBinCut[0]) corrPt *= RES_Corr1_010_f[algBin]->Eval(jtPtCheck);
      else if(hiBin <= hiBinCut[1]) corrPt *= RES_Corr1_1030_f[algBin]->Eval(jtPtCheck);
      else if(hiBin <= hiBinCut[2]) corrPt *= RES_Corr1_3050_f[algBin]->Eval(jtPtCheck);
      else corrPt *= RES_Corr1_50100_f[algBin]->Eval(jtPtCheck);
    }
    else if(sType == kPPDATA || sType == kPPMC){
      corrPt = RES_Corr0_PP_f[algBin]->Eval(jtPt);

      jtPtCheck = jtPt/corrPt;
      if(jtPtCheck < 25) jtPtCheck = 25;

      corrPt *= RES_Corr1_PP_f[algBin]->Eval(jtPtCheck);

      jtPtCheck = jtPt/corrPt;
      if(jtPtCheck < 25) jtPtCheck = 25;

      corrPt *= RES_Corr2_PP_f[algBin]->Eval(jtPtCheck);
    }
  }

  return 1/corrPt;
}

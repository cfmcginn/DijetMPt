//=============================================
// Author: Chris McGinn
//
// DiJet Plotter                                                              
//
//=============================================

#include "TDatime.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include <iostream>
#include "TAttFill.h"
#include "TColor.h"

#include "TF1.h"

#include <string>

TFile* plotFile_p = 0;

const Int_t nSumAlg = 40;

const char* algType[nSumAlg] = {"Pu3Calo", "Pu4Calo", "Pu5Calo", "Pu3PF", "Vs2Calo", "Vs3Calo", "Vs4Calo", "Vs5Calo", "Vs3PF", "Pu3CaloFrag", "Vs2CaloFrag", "Vs3CaloFrag", "Vs4CaloFrag", "Vs5CaloFrag", "Vs3PFFrag", "Pu3CaloRes", "Vs2CaloRes", "Vs3CaloRes", "Vs4CaloRes", "Vs5CaloRes", "Vs3PFRes", "Pu3CaloSwap", "Vs2CaloSwap", "Vs3CaloSwap", "Vs4CaloSwap", "Vs5CaloSwap", "Vs3PFSwap", "Pu3CaloResSpill", "Vs2CaloResSpill", "Vs3CaloResSpill", "Vs4CaloResSpill", "Vs5CaloResSpill", "Vs2CaloSwapSpill", "Vs3CaloSwapSpill", "Vs4CaloSwapSpill", "Vs5CaloSwapSpill", "T2", "T3", "T4", "T5"};

const char* algTypePP[nSumAlg] = {"Vs3Calo", "Vs4Calo", "Vs5Calo", "Vs3PF", "Vs2Calo", "Vs3Calo", "Vs4Calo", "Vs5Calo", "Vs3PF", "Vs3CaloFrag", "Vs2CaloFrag", "Vs3CaloFrag", "Vs4CaloFrag", "Vs5CaloFrag", "Vs3PFFrag", "Vs3CaloRes", "Vs2CaloRes", "Vs3CaloRes", "Vs4CaloRes", "Vs5CaloRes", "Vs3PFRes", "Vs3CaloSwap", "Vs2CaloSwap", "Vs3CaloSwap", "Vs4CaloSwap", "Vs5CaloSwap", "Vs3PFSwap", "Vs3CaloResSpill", "Vs2CaloResSpill", "Vs3CaloResSpill", "Vs4CaloResSpill", "Vs5CaloResSpill", "Vs2CaloSwapSpill", "Vs3CaloSwapSpill", "Vs4CaloSwapSpill", "Vs5CaloSwapSpill", "T2", "T3", "T4", "T5"};

const char* truthType[nSumAlg] = {"T3", "T4", "T5", "T3", "T2", "T3", "T4", "T5", "T3", "T3", "T2", "T3", "T4", "T5", "T3", "T3", "T2", "T3", "T4", "T5", "T3", "T3", "T2", "T3", "T4", "T5", "T3", "T3", "T2", "T3", "T4", "T5", "T2", "T3", "T4", "T5", "T2", "T3", "T4", "T5"};

const char* radString[nSumAlg] = {"0.3", "0.4", "0.5", "0.3", "0.2", "0.3", "0.4", "0.5", "0.3", "0.3", "0.2", "0.3", "0.4", "0.5", "0.3", "0.3", "0.2", "0.3", "0.4", "0.5", "0.3", "0.3", "0.2", "0.3", "0.4", "0.5", "0.3", "0.3", "0.2", "0.3", "0.4", "0.5", "0.2", "0.3", "0.4", "0.5", "0.2", "0.3", "0.4", "0.5"};

const char* radString2[nSumAlg] = {"3", "4", "5", "3", "2", "3", "4", "5", "3", "3", "2", "3", "4", "5", "3", "3", "2", "3", "4", "5", "3", "3", "2", "3", "4", "5", "3", "3", "2", "3", "4", "5", "2", "3", "4", "5", "2", "3", "4", "5"};

const char* puVsString[nSumAlg] = {"PuCalo", "PuCalo", "PuCalo", "PuPF", "VsCalo", "VsCalo", "VsCalo", "VsCalo", "VsPF", "PuCaloFrag", "VsCaloFrag", "VsCaloFrag", "VsCaloFrag", "VsCaloFrag", "VsPFFrag", "PuCaloRes", "VsCaloRes", "VsCaloRes", "VsCaloRes", "VsCaloRes", "VsPFRes", "PuCaloSwap", "VsCaloSwap", "VsCaloSwap", "VsCaloSwap", "VsCaloSwap", "VsPFSwap", "PuCaloResSpill", "VsCaloResSpill", "VsCaloResSpill", "VsCaloResSpill", "VsCaloResSpill", "VsCaloSwapSpill", "VsCaloSwapSpill", "VsCaloSwapSpill", "VsCaloSwapSpill", "T", "T", "T", "T"};

Bool_t isAll(const char* CNCR);

Bool_t sameSign(Double_t num1, Double_t num2)
{
  if((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0)) return true;

  return false;
}


Double_t quadSum(Double_t one, Double_t two)
{
  return TMath::Sqrt(one*one + two*two);
}


void drawPatch(float x1, float y1, float x2, float y2){
  TLegend *t1=new TLegend(x1,y1,x2,y2);
  t1->SetFillColor(kWhite);
  t1->SetBorderSize(0);
  t1->SetFillStyle(1001);
  t1->Draw("");

  return;
}


void makePatch(const char* CNCR){
  if(!strcmp(CNCR, "R") || !strcmp(CNCR, "R2") || !strcmp(CNCR, "RU") || !strcmp(CNCR, "RD") || !strcmp(CNCR, "Eta") || !strcmp(CNCR, "EtaD") || !strcmp(CNCR, "EtaU") || !strcmp(CNCR, "Phi") || !strcmp(CNCR, "PhiD") || !strcmp(CNCR, "PhiU"))
    drawPatch(.74, .001, .90, .25);
  else if(!strcmp(CNCR, "C") || !strcmp(CNCR, "NC"))
    drawPatch(.74, .001, .90, .25);
  else
    drawPatch(.74, .001, .90, .25);

  return;
}


void drawNum(const char* CNCR){
  TLatex* temp = new TLatex();
  temp->SetNDC();
  temp->SetTextFont(43);
  temp->SetTextSizePixels(27);

  if(isAll(CNCR))
    temp->DrawLatex(.0000001, .945, "1.5");
  else if(!strcmp(CNCR, "C") || !strcmp(CNCR, "NC"))
    temp->DrawLatex(.0001, .945, "0.4");
  else
    temp->DrawLatex(.0001, .94, "0.4");

  return;
}


Bool_t isR(const char* CNCR)
{
  if(!strcmp(CNCR, "R") || !strcmp(CNCR, "RD") || !strcmp(CNCR, "RU") || !strcmp(CNCR, "R2") || !strcmp(CNCR, "R2D") || !strcmp(CNCR, "R2U") || !strcmp(CNCR, "RFOR") || !strcmp(CNCR, "RFORD") || !strcmp(CNCR, "RFORU") || !strcmp(CNCR, "RFORMID") || !strcmp(CNCR, "RFORMIDD") || !strcmp(CNCR, "RFORMIDU") || !strcmp(CNCR, "RFORFOR") || !strcmp(CNCR, "RFORFORD") || !strcmp(CNCR, "RFORFORU") || !strcmp(CNCR, "RMID") || !strcmp(CNCR, "RMIDD") || !strcmp(CNCR, "RMIDU") || !strcmp(CNCR, "RCut") || !strcmp(CNCR, "RCutD") || !strcmp(CNCR, "RCutU") || !strcmp(CNCR, "RCutEta") || !strcmp(CNCR, "RCutEtaD") || !strcmp(CNCR, "RCutEtaU") || !strcmp(CNCR, "RCutPhi") || !strcmp(CNCR, "RCutPhiD") || !strcmp(CNCR, "RCutPhiU")) return true;
  else return false;
}


Bool_t isEta(const char* CNCR)
{
  if(!strcmp(CNCR, "Eta") || !strcmp(CNCR, "EtaU") || !strcmp(CNCR, "EtaD") || !strcmp(CNCR, "EtaCut") || !strcmp(CNCR, "EtaCutU") || !strcmp(CNCR, "EtaCutD")) return true;
  else return false;
}


Bool_t isPhi(const char* CNCR)
{
  if(!strcmp(CNCR, "Phi") || !strcmp(CNCR, "PhiU") || !strcmp(CNCR, "PhiD") || !strcmp(CNCR, "PhiCut") || !strcmp(CNCR, "PhiCutU") || !strcmp(CNCR, "PhiCutD")) return true;
  else return false;
}


Bool_t isInc(const char* CNCR)
{
  if(!strcmp(CNCR, "R") || !strcmp(CNCR, "R_320") || !strcmp(CNCR, "R_340") || !strcmp(CNCR, "RFOR") || !strcmp(CNCR, "RFORMID") || !strcmp(CNCR, "RFORFOR") || !strcmp(CNCR, "RMID") || !strcmp(CNCR, "Eta") || !strcmp(CNCR, "Phi") || !strcmp(CNCR, "RCut") || !strcmp(CNCR, "RCutEta") || !strcmp(CNCR, "RCutPhi") || !strcmp(CNCR, "EtaCut") || !strcmp(CNCR, "PhiCut")) return true;
  else return false;
}

Bool_t isU(const char* CNCR)
{
  if(!strcmp(CNCR, "RU") || !strcmp(CNCR, "RU_320") || !strcmp(CNCR, "RU_340") || !strcmp(CNCR, "R2U") || !strcmp(CNCR, "RFORU") || !strcmp(CNCR, "RFORMIDU") || !strcmp(CNCR, "RFORFORU") || !strcmp(CNCR, "RMIDU") || !strcmp(CNCR, "EtaU") || !strcmp(CNCR, "PhiU") || !strcmp(CNCR, "RCutU") || !strcmp(CNCR, "RCutEtaU") || !strcmp(CNCR, "RCutPhiU") || !strcmp(CNCR, "EtaCutU") || !strcmp(CNCR, "PhiCutU")) return true;
  else return false;
}


Bool_t isD(const char* CNCR)
{
  if(!strcmp(CNCR, "RD") || !strcmp(CNCR, "RD_320") || !strcmp(CNCR, "RD_340") || !strcmp(CNCR, "R2D")|| !strcmp(CNCR, "RFORD") || !strcmp(CNCR, "RFORMIDD") || !strcmp(CNCR, "RFORFORD") || !strcmp(CNCR, "RMIDD") || !strcmp(CNCR, "EtaD") || !strcmp(CNCR, "PhiD") || !strcmp(CNCR, "RCutD") || !strcmp(CNCR, "RCutEtaD") || !strcmp(CNCR, "RCutPhiD") || !strcmp(CNCR, "EtaCutD") || !strcmp(CNCR, "PhiCutD")) return true;
  else return false;
}


Bool_t isAll(const char* CNCR)
{
  if(isR(CNCR) || isPhi(CNCR) || isEta(CNCR)) return true;
  else return false;
}


void drawBin(const char* CNCR){
  TLatex* temp = new TLatex();
  temp->SetNDC();
  temp->SetTextFont(43);
  temp->SetTextSizePixels(28);

  if(isInc(CNCR)) temp->DrawLatex(.30, .16, "A_{J} Inclusive");
  else if(isU(CNCR)) temp->DrawLatex(.30, .16, "A_{J}>0.22");
  else if(isD(CNCR)) temp->DrawLatex(.30, .16, "A_{J}<0.22");
  else if(!strcmp(CNCR, "C")) temp->DrawLatex(.15, .16, "In-cone, #DeltaR < 0.8");
  else if(!strcmp(CNCR, "NC")) temp->DrawLatex(.15, .16, "Out-cone, #DeltaR > 0.8");
  else if(!strcmp(CNCR, "C0")) temp->DrawLatex(.15, .16, "#DeltaR < 0.5");
  else if(!strcmp(CNCR, "C1")) temp->DrawLatex(.15, .16, "0.5 < #DeltaR < 1.0");
  else if(!strcmp(CNCR, "C2")) temp->DrawLatex(.15, .16, "1.0 < #DeltaR < 1.5");
  else if(!strcmp(CNCR, "C3")) temp->DrawLatex(.15, .16, "1.5 < #DeltaR");

  //Quark In-Venn etc. label here
  temp->DrawLatex(.30, .06, "3rd Jet p_{T}>50 GeV/c");

  if(!strcmp(CNCR, "RCut") || !strcmp(CNCR, "RCutU") || !strcmp(CNCR, "RCutD")) temp->DrawLatex(.30, .08, "|#Delta#phi(#Delta#eta)|_{trk,jet}<1.0");
  else if(!strcmp(CNCR, "EtaCut") || !strcmp(CNCR, "EtaCutU") || !strcmp(CNCR, "EtaCutD") || !strcmp(CNCR, "RCutPhi") || !strcmp(CNCR, "RCutPhiU") || !strcmp(CNCR, "RCutPhiD")) temp->DrawLatex(.30, .08, "|#Delta#phi|_{trk,jet}<1.0");
  else if(!strcmp(CNCR, "PhiCut") || !strcmp(CNCR, "PhiCutU") || !strcmp(CNCR, "PhiCutD") || !strcmp(CNCR, "RCutEta") || !strcmp(CNCR, "RCutEtaU") || !strcmp(CNCR, "RCutEtaD")) temp->DrawLatex(.30, .08, "|#Delta#eta|_{trk,jet}<1.0");

  return;
}


void handsomeTH1( TH1 *a=0, Int_t col =1, Float_t size=1, Int_t markerstyle=20)
{
  a->SetMarkerColor(col);
  a->SetMarkerSize(size);
  a->SetMarkerStyle(markerstyle);
  a->SetLineColor(col);
  a->GetYaxis()->SetTitleOffset(1.25);
  a->GetXaxis()->CenterTitle();
  a->GetYaxis()->CenterTitle();

  return;
}

void handsomeTH1N( TH1 *a=0, Int_t col =1)
{
  handsomeTH1(a,col);
  a->Scale(1./a->GetEntries());

  return;
}


void niceTH1(TH1F* uglyTH1, float max , float min, float ndivX, float ndivY)
{
  handsomeTH1(uglyTH1);
  uglyTH1->SetMaximum(max);
  uglyTH1->SetMinimum(min);
  uglyTH1->SetNdivisions(ndivX);
  uglyTH1->SetNdivisions(ndivY, "Y");

  return;
}


void claverCanvasSaving(TCanvas* c, TString s,TString format="gif")
{
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
}


void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns, const Int_t rows, const Float_t leftOffset, const Float_t bottomOffset, const Float_t leftMargin, const Float_t bottomMargin, const Float_t edge, const char* CNCR = ""){
  if(canv==0){
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
  Float_t PadWidth = (1.0-leftOffset)/((1.0/(1.0-leftMargin)) + (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight = (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) + (1.0/(1.0-edge))+(Float_t)rows-2.0);
  Xlow[0] = leftOffset;
  Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

  for(Int_t i=1; i<columns-1; i++){
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
  }
  Int_t ct = 0;
  for(Int_t i=rows-2; i>0; i--){
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }

  TString padName;
  for(Int_t i=0; i<columns; i++){
    for(Int_t j=0; j<rows; j++){
      canv->cd();
      padName = Form("p_%d_%d",i,j);

      if(i==0 && j==1){
        pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j],Xup[i]*0.890,Yup[j]*0.935);
      }
      else if(i==0 && j==0){
        if(strcmp(CNCR, "") ==0)
          pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j]*0.933,Xup[i],Yup[j]);
        else
          pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j]*0.934,Xup[i],Yup[j]*.9985);
      }
      else if(i==1 && j==1){
        if(strcmp(CNCR, "") == 0)
          pad[i][j] = new TPad(padName.Data(),padName.Data(),0.815*Xlow[i],Ylow[j],Xup[i],Yup[j]*.995);
        else
          pad[i][j] = new TPad(padName.Data(),padName.Data(),0.80*Xlow[i],Ylow[j],Xup[i],Yup[j]*.995);
      }
      else if(j == 0){
        pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j],Xup[i],Yup[j]*.998);
      }
      else pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j],Xup[i],Yup[j]*.995);

      if(i==0){
        if(j == 0)
          pad[i][j]->SetLeftMargin(leftMargin*1.2);
       else
          pad[i][j]->SetLeftMargin(leftMargin);
      }
      else if(i==1 && j==1){
        if(strcmp(CNCR, "") == 0)
          pad[i][j]->SetLeftMargin(PadWidth);
        else
          pad[i][j]->SetLeftMargin(PadWidth*.65);
      }
      else pad[i][j]->SetLeftMargin(0);

      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
      else pad[i][j]->SetRightMargin(0);

      if(j==0){
        if(i==0)pad[i][j]->SetTopMargin(edge);
        else pad[i][j]->SetTopMargin(edge);
      }
      else pad[i][j]->SetTopMargin(0);

      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
      else if(i==0 && j==0) pad[i][j]->SetBottomMargin(0.17*PadHeight);
      else pad[i][j]->SetBottomMargin(0);


      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);

    }
  }
  pad[0][0]->cd();

  return;
}


void grabHistForStackPP(TFile* histFile_p, std::string gorr, std::string alg, std::string projMult, std::string CNCR, std::string Corr, std::string fileTag, TH1F* histPP_p[7])
{  
  std::string FPT[7] = {"0_5", "5_1", "1_2", "2_4", "4_8", "8_100", "F"};

  for(Int_t histIter = 0; histIter < 7; histIter++){
    histPP_p[histIter] = (TH1F*)histFile_p->Get(Form("%s%s%sA%s%s%s_PP_%s_h", gorr.c_str(), alg.c_str(), projMult.c_str(), CNCR.c_str(), FPT[histIter].c_str(), Corr.c_str(), fileTag.c_str()));
  }

  return;
}


void grabSpillForStackPP(TFile* histFile_p, std::string ResSwap, std::string algR, TH1F* histPP_p[7])
{
  std::string FPT[7] = {"0_5", "5_1", "1_2", "2_4", "4_8", "8_100", "F"};

  for(Int_t histIter = 1; histIter < 7; histIter++){
    histPP_p[histIter] = (TH1F*)histFile_p->Get(Form("h_spillover_param_%sCalo%sProj%s", algR.c_str(), ResSwap.c_str(), FPT[histIter].c_str()));
  }

  return;
}


void grabHistForStackPbPbCent(TFile* histFile_p, std::string gorr, std::string alg, std::string projMult, std::string CNCR, std::string Corr, const std::string cent, std::string fileTag, TH1F* hist_p[7])
{  
  std::string FPT[7] = {"0_5", "5_1", "1_2", "2_4", "4_8", "8_100", "F"};

  for(Int_t histIter = 0; histIter < 7; histIter++){
    if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")) hist_p[histIter] = (TH1F*)histFile_p->Get(Form("%s%s%sA%s%s%s_%s_%s_h", gorr.c_str(), alg.c_str(), projMult.c_str(), CNCR.c_str(), FPT[histIter].c_str(), Corr.c_str(), cent.c_str(), fileTag.c_str()));
    else hist_p[histIter] = (TH1F*)histFile_p->Get(Form("%s%s%sA%s%s%s_%s_%s_h", gorr.c_str(), alg.c_str(), projMult.c_str(), CNCR.c_str(), FPT[histIter].c_str(), Corr.c_str(), cent.c_str(), fileTag.c_str()));
  }

  return;
}


void grabHistForStackPbPb(TFile* histFile_p, std::string gorr, std::string alg, std::string projMult, std::string CNCR, std::string Corr, std::string fileTag, TH1F* hist1_p[7], TH1F* hist2_p[7], TH1F* hist3_p[7], TH1F* hist4_p[7])
{  
  std::string FPT[7] = {"0_5", "5_1", "1_2", "2_4", "4_8", "8_100", "F"};

  for(Int_t histIter = 0; histIter < 7; histIter++){
    if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")){
      hist1_p[histIter] = (TH1F*)histFile_p->Get(Form("%s%s%sA%s%s%s_50100_%s_h", gorr.c_str(), alg.c_str(), projMult.c_str(), CNCR.c_str(), FPT[histIter].c_str(), Corr.c_str(), fileTag.c_str()));
      hist2_p[histIter] = (TH1F*)histFile_p->Get(Form("%s%s%sA%s%s%s_3050_%s_h", gorr.c_str(), alg.c_str(), projMult.c_str(), CNCR.c_str(), FPT[histIter].c_str(), Corr.c_str(), fileTag.c_str()));
      hist3_p[histIter] = (TH1F*)histFile_p->Get(Form("%s%s%sA%s%s%s_1030_%s_h", gorr.c_str(), alg.c_str(), projMult.c_str(), CNCR.c_str(), FPT[histIter].c_str(), Corr.c_str(), fileTag.c_str()));
      hist4_p[histIter] = (TH1F*)histFile_p->Get(Form("%s%s%sA%s%s%s_010_%s_h", gorr.c_str(), alg.c_str(), projMult.c_str(), CNCR.c_str(), FPT[histIter].c_str(), Corr.c_str(), fileTag.c_str()));
    }
    else{
      hist1_p[histIter] = (TH1F*)histFile_p->Get(Form("%s%s%sA%s%s%s_30100_%s_h", gorr.c_str(), alg.c_str(), projMult.c_str(), CNCR.c_str(), FPT[histIter].c_str(), Corr.c_str(), fileTag.c_str()));
      hist2_p[histIter] = (TH1F*)histFile_p->Get(Form("%s%s%sA%s%s%s_030_%s_h", gorr.c_str(), alg.c_str(), projMult.c_str(), CNCR.c_str(), FPT[histIter].c_str(), Corr.c_str(), fileTag.c_str()));
    }
  }

  return;
}


void grabSpillForStackPbPb(TFile* histFile_p, std::string ResSwap, std::string algR, std::string CNCR, TH1F* hist1_p[7], TH1F* hist2_p[7], TH1F* hist3_p[7], TH1F* hist4_p[7])
{
  std::string FPT[7] = {"0_5", "5_1", "1_2", "2_4", "4_8", "8_100", "F"};

  for(Int_t histIter = 1; histIter < 7; histIter++){
    if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")){
      std::cout << Form("h_spillover_param_%sCalo%sProj%s_50100", algR.c_str(), ResSwap.c_str(), FPT[histIter].c_str()) << std::endl;

      hist1_p[histIter] = (TH1F*)histFile_p->Get(Form("h_spillover_param_%sCalo%sProj%s_50100", algR.c_str(), ResSwap.c_str(), FPT[histIter].c_str()));
      hist2_p[histIter] = (TH1F*)histFile_p->Get(Form("h_spillover_param_%sCalo%sProj%s_3050", algR.c_str(), ResSwap.c_str(), FPT[histIter].c_str()));
      hist3_p[histIter] = (TH1F*)histFile_p->Get(Form("h_spillover_param_%sCalo%sProj%s_1030", algR.c_str(), ResSwap.c_str(), FPT[histIter].c_str()));
      hist4_p[histIter] = (TH1F*)histFile_p->Get(Form("h_spillover_param_%sCalo%sProj%s_010", algR.c_str(), ResSwap.c_str(), FPT[histIter].c_str()));
    }
    else{
      std::cout << Form("h_spillover_param_%sCalo%sProj%s_30100", algR.c_str(), ResSwap.c_str(), FPT[histIter].c_str()) << std::endl;
      hist1_p[histIter] = (TH1F*)histFile_p->Get(Form("h_spillover_param_%sCalo%sProj%s_30100", algR.c_str(), ResSwap.c_str(), FPT[histIter].c_str()));
      hist2_p[histIter] = (TH1F*)histFile_p->Get(Form("h_spillover_param_%sCalo%sProj%s_030", algR.c_str(), ResSwap.c_str(), FPT[histIter].c_str()));
    }
  }
  return;
}


Double_t sumYForPTStack(Double_t dIn = 0, Double_t comp1 = 0, Double_t comp2 = 0, Double_t comp3 = 0, Double_t comp4 = 0, Double_t comp5 = 0)
{
  Double_t dOut = dIn;

  if(sameSign(comp1, dOut))
    dOut += comp1;

  if(sameSign(comp2, dOut))
    dOut += comp2;

  if(sameSign(comp3, dOut))
    dOut += comp3;

  if(sameSign(comp4, dOut))
    dOut += comp4;

  if(sameSign(comp5, dOut))
    dOut += comp5;

  return dOut;
}


void makeHistForPtStack(TH1F* h_p[7], Int_t pos = 4, const char* Tight = "", const char* CNCR = "", const std::string projMult = "", const Bool_t isDoubleSub = false)
{
  Int_t nBins = 4;

  if(!strcmp(CNCR, "R2") || !strcmp(CNCR, "R2D") || !strcmp(CNCR, "R2U")) nBins = 9;
  else if(isAll(CNCR)) nBins = 10;
  else if(!strcmp(CNCR, "R_320") || !strcmp(CNCR, "RD_320") || !strcmp(CNCR, "RU_320")) nBins = 40;
  else if(!strcmp(CNCR, "R_340") || !strcmp(CNCR, "RD_340") || !strcmp(CNCR, "RU_340")) nBins = 40;
  else if(strcmp(Tight, "") != 0) nBins = 8;

  for(Int_t iter = 0; iter < nBins; iter++){
    //    h_p[0]->SetBinContent(iter + 1, sumYForPTStack(h_p[0]->GetBinContent(iter+1), h_p[1]->GetBinContent(iter+1), h_p[2]->GetBinContent(iter+1), h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1), h_p[5]->GetBinContent(iter+1)));

    h_p[1]->SetBinContent(iter + 1, sumYForPTStack(h_p[1]->GetBinContent(iter+1), h_p[2]->GetBinContent(iter+1), h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1), h_p[5]->GetBinContent(iter+1)));
    h_p[2]->SetBinContent(iter + 1, sumYForPTStack(h_p[2]->GetBinContent(iter+1), h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1), h_p[5]->GetBinContent(iter+1)));
    h_p[3]->SetBinContent(iter + 1, sumYForPTStack(h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1), h_p[5]->GetBinContent(iter+1)));
    h_p[4]->SetBinContent(iter + 1, sumYForPTStack(h_p[4]->GetBinContent(iter+1), h_p[5]->GetBinContent(iter+1)));
  }

  const char* xTitle;

  if(isR(CNCR)) xTitle = "#DeltaR";
  else if(!strcmp(CNCR, "R_320") || !strcmp(CNCR, "RD_320") || !strcmp(CNCR, "RU_320")) xTitle = "#DeltaR";
  else if(!strcmp(CNCR, "R_340") || !strcmp(CNCR, "RD_340") || !strcmp(CNCR, "RU_340")) xTitle = "#DeltaR";
  else if(isEta(CNCR)) xTitle = "#Delta#eta";
  else if(isPhi(CNCR)) xTitle = "#Delta#phi";
  else xTitle = "A_{J}";

  h_p[4]->SetXTitle(xTitle);
  h_p[3]->SetXTitle(xTitle);
  h_p[2]->SetXTitle(xTitle);
  h_p[1]->SetXTitle(xTitle);
  h_p[0]->SetXTitle(xTitle);
  h_p[5]->SetXTitle(xTitle);
  h_p[6]->SetXTitle(xTitle);

  if(pos == 1 && !isDoubleSub){
    if(!strcmp(projMult.c_str(), "Proj")){
      h_p[4]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
      h_p[3]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
      h_p[2]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
      h_p[1]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
      h_p[0]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
      h_p[5]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
      h_p[6]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
    }
    else if(!strcmp(projMult.c_str(), "Mult")){
      h_p[4]->SetYTitle("   #Delta_{Mult}");
      h_p[3]->SetYTitle("   #Delta_{Mult}");
      h_p[2]->SetYTitle("   #Delta_{Mult}");
      h_p[1]->SetYTitle("   #Delta_{Mult}");
      h_p[0]->SetYTitle("   #Delta_{Mult}");
      h_p[5]->SetYTitle("   #Delta_{Mult}");
      h_p[6]->SetYTitle("   #Delta_{Mult}");
    }
  }
  else if(pos == 1){
    if(!strcmp(projMult.c_str(), "Proj")){
      h_p[4]->SetYTitle("   #Delta(<#slash{p}_{T}^{||}> (GeV/c))");
      h_p[3]->SetYTitle("   #Delta(<#slash{p}_{T}^{||}> (GeV/c))");
      h_p[2]->SetYTitle("   #Delta(<#slash{p}_{T}^{||}> (GeV/c))");
      h_p[1]->SetYTitle("   #Delta(<#slash{p}_{T}^{||}> (GeV/c))");
      h_p[0]->SetYTitle("   #Delta(<#slash{p}_{T}^{||}> (GeV/c))");
      h_p[5]->SetYTitle("   #Delta(<#slash{p}_{T}^{||}> (GeV/c))");
      h_p[6]->SetYTitle("   #Delta(<#slash{p}_{T}^{||}> (GeV/c))");
    }
    else if(!strcmp(projMult.c_str(), "Mult")){
      h_p[4]->SetYTitle("   #Delta_{Mult}");
      h_p[3]->SetYTitle("   #Delta_{Mult}");
      h_p[2]->SetYTitle("   #Delta_{Mult}");
      h_p[1]->SetYTitle("   #Delta_{Mult}");
      h_p[0]->SetYTitle("   #Delta_{Mult}");
      h_p[5]->SetYTitle("   #Delta_{Mult}");
      h_p[6]->SetYTitle("   #Delta_{Mult}");
    }
  }

  return;
}

void drawHistToPTStack(TH1F* drawHist_p, Int_t color, const char* drawOpt, const std::string projMult, Bool_t isSub = false, Bool_t isDoubleSub = false, const char* CNCR = "")
{
  Float_t lowBound = -10;
  //  lowBound = -7.999;
  if(!strcmp(projMult.c_str(), "Mult")) lowBound = -5.;
  //  if(!strcmp(CNCR, "RD") || !strcmp(CNCR, "R")) lowBound = -5.;

  Float_t hiBound = 7.999;
  if(!strcmp(CNCR, "RD") || !strcmp(CNCR, "R")) hiBound = 3.999;

  if(isSub){
    if(isAll(CNCR) && !isDoubleSub) niceTH1(drawHist_p, hiBound, lowBound, 505, 503);
    else if(isAll(CNCR) && isDoubleSub) niceTH1(drawHist_p, 14.999, -10, 505, 503);
    else if(!strcmp(CNCR, "C0") || !strcmp(CNCR, "C1") || !strcmp(CNCR, "C2") || !strcmp(CNCR, "C3")) niceTH1(drawHist_p, 14.999, -20., 505, 403);
    else if(isDoubleSub) niceTH1(drawHist_p, 29.999, -30., 505, 406);
    else if(!strcmp(CNCR, "R_320") || !strcmp(CNCR, "RD_320") || !strcmp(CNCR, "RU_320")) niceTH1(drawHist_p, 3.999, -4., 505, 406);
    else if(!strcmp(CNCR, "R_340") || !strcmp(CNCR, "RD_340") || !strcmp(CNCR, "RU_340")) niceTH1(drawHist_p, 2.999, -3., 505, 406);
    else niceTH1(drawHist_p, 39.999, -40., 505, 406);
  }
  else{
    if(isAll(CNCR) && isDoubleSub) niceTH1(drawHist_p, 14.999, -10, 505, 503);
    else if(!isAll(CNCR) && isDoubleSub) niceTH1(drawHist_p, 29.999, -30, 505, 406);
    else if(isR(CNCR) && !isDoubleSub){
      if(!strcmp(CNCR, "RU")) niceTH1(drawHist_p, 19.999, -50.0, 505, 503);
      else niceTH1(drawHist_p, 9.999, -30., 505, 503);

      //      if(!strcmp(CNCR, "RU")) niceTH1(drawHist_p, 9.999, -5., 505, 503);
      //      else niceTH1(drawHist_p, 9.999, -5., 505, 503);
    }
    else if (!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")) niceTH1(drawHist_p, 59.999, -60., 505, 406);
    else if (!strcmp(CNCR, "R_320") || !strcmp(CNCR, "RD_320") || !strcmp(CNCR, "RU_320")) niceTH1(drawHist_p, 5.999, -4., 505, 406);
    else if (!strcmp(CNCR, "R_340") || !strcmp(CNCR, "RD_340") || !strcmp(CNCR, "RU_340")) niceTH1(drawHist_p, 2.999, -3., 505, 406);
  }

  drawHist_p->GetYaxis()->SetTitleOffset(3.0);

  drawHist_p->SetFillColor(color);
  drawHist_p->SetMarkerStyle(6);
  drawHist_p->SetMarkerSize(.5);
  drawHist_p->GetXaxis()->SetTitleOffset(1.6);
  drawHist_p->DrawCopy(drawOpt);
  drawHist_p->DrawCopy("E1 SAME");

  return;
}


void drawFullStack(TH1F* h_p[7], Int_t color, Int_t style, const std::string projMult, Int_t pos = 2, Bool_t isSub = false, const char* CNCR = "", Bool_t isDoubleSub = false, Bool_t isRDEP = false)
{
  Float_t lowBound = -10;
  if(!strcmp(projMult.c_str(), "Mult")) lowBound = -5.;

  //  drawHistToPTStack(h_p[0], kMagenta + 4, "E1 HIST", projMult, isSub, isDoubleSub, CNCR);
  drawHistToPTStack(h_p[1], kBlue - 9, "E1 HIST", projMult, isSub, isDoubleSub, CNCR);
  drawHistToPTStack(h_p[2], kYellow - 9, "E1 HIST SAME", projMult, isSub, isDoubleSub, CNCR);
  drawHistToPTStack(h_p[3], kOrange + 1, "E1 HIST SAME", projMult, isSub, isDoubleSub, CNCR);
  //drawHistToPTStack(h_p[3], kGray, "E1 HIST SAME", projMult, isSub, isDoubleSub, CNCR);
  drawHistToPTStack(h_p[4], kGreen + 3, "E1 HIST SAME", projMult, isSub, isDoubleSub, CNCR);
  //drawHistToPTStack(h_p[4], kGray + 1, "E1 HIST SAME", projMult, isSub, isDoubleSub, CNCR);
  drawHistToPTStack(h_p[5], kRed + 1, "E1 HIST SAME", projMult, isSub, isDoubleSub, CNCR);
  //drawHistToPTStack(h_p[5], kGray+2, "E1 HIST SAME", projMult, isSub, isDoubleSub, CNCR);


  if(isSub){
    if(isAll(CNCR))
      niceTH1(h_p[6], 9.999, lowBound, 505, 403);
    else
      niceTH1(h_p[6], 59.999, -60., 505, 406);
  }

  h_p[6]->SetFillColor(color);
  h_p[6]->SetMarkerStyle(style);
  h_p[6]->SetLineWidth(3);

  h_p[6]->DrawCopy("SAME E1");

  h_p[6]->SetLineWidth(1);

  Int_t nBins = 9;
  if(!strcmp(CNCR, "R2") || !strcmp(CNCR, "R2U") || !strcmp(CNCR, "R2D")) nBins = 8;
  if(!strcmp(CNCR, "R_320") || !strcmp(CNCR, "RU_320") || !strcmp(CNCR, "RD_320")) nBins = 19;
  if(!strcmp(CNCR, "R_340") || !strcmp(CNCR, "RU_340") || !strcmp(CNCR, "RD_340")) nBins = 39;

  if(isSub == false){
    if(isAll(CNCR) || isInc(CNCR) || isU(CNCR) || isD(CNCR)){
      for(Int_t hIter = 0; hIter < nBins; hIter++){
	h_p[6]->SetBinContent(hIter+2, h_p[6]->GetBinContent(hIter+1) + h_p[6]->GetBinContent(hIter+2));
      }
      if(pos == 1 || (isRDEP && pos < 5))
	h_p[6]->SetLineStyle(2);

      if(!isDoubleSub && strcmp(projMult.c_str(), "Mult") != 0){
	h_p[6]->DrawCopy("SAME HIST C");
	

	Float_t amp = h_p[6]->GetBinContent(1);
	TF1* r1 = new TF1("r1", "[0]/x", 0.05, 1.95);
	r1->SetParameter(0, amp*0.1);
	TF1* r2 = new TF1("r2", "[0]/(x*x)", 0.05, 1.95);
	r2->SetParameter(0, amp*0.1*0.1);
	TF1* rHalf = new TF1("rHalf", "[0]/TMath::Sqrt(x)", 0.05, 1.95);
	rHalf->SetParameter(0, amp*TMath::Sqrt(0.1));

	TF1* rFit = new TF1("rFit", "[0]/TMath::Power(x+[1], [2])", 0.05, 1.85);
	rFit->SetParameter(0, amp*0.1);
	rFit->SetParameter(1, 0);
	rFit->SetParameter(2, 1);

	h_p[6]->Fit(rFit, "QNMR", "");

	//	r1->DrawCopy("SAME");
	//	r2->DrawCopy("SAME");
	//	rHalf->DrawCopy("SAME");
	//	rFit->DrawCopy("SAME");


	TF1* rFit2 = new TF1("rFit2", "[0]/TMath::Power(x+[1], [2])", 0.05, 1.95);
	rFit2->SetParameter(0, rFit->GetParameter(0));
	rFit2->SetParameter(1, rFit->GetParameter(1));
	rFit2->SetParameter(2, rFit->GetParameter(2));
	//	rFit2->DrawCopy("SAME");
      }

    }
  }

  return;
}


void makeSysError(Float_t sysArr[], TH1F* hist_p, Bool_t bottom, const std::string CNCR)
{
  Float_t length1 = 2.0;
  Float_t length2 = 0.01;
  if(bottom && strcmp(CNCR.c_str(), "") != 0){
    length1 = 0.25;
    length2 = 0.025;
  }
  else if(isR(CNCR.c_str())){
    length1 = 1.0;
  }

  for(Int_t iter = 0; iter < hist_p->GetNbinsX(); iter++){
    Float_t yVal = hist_p->GetBinContent(iter+1);
    Float_t sys = sysArr[iter];
    TLine* l = new TLine(hist_p->GetBinLowEdge(iter+1) + length2, yVal - TMath::Sqrt(sys*sys), hist_p->GetBinLowEdge(iter+2) - length2, yVal - TMath::Sqrt(sys*sys));
    l->SetLineColor(1);
    l->Draw();
    l->DrawLine(hist_p->GetBinLowEdge(iter+1) + length2, yVal - TMath::Sqrt(sys*sys), hist_p->GetBinLowEdge(iter+1) + length2, yVal - TMath::Sqrt(sys*sys) + length1);
    l->DrawLine(hist_p->GetBinLowEdge(iter + 2) - length2, yVal - TMath::Sqrt(sys*sys), hist_p->GetBinLowEdge(iter + 2) - length2, yVal - TMath::Sqrt(sys*sys) + length1);
    l->DrawLine(hist_p->GetBinLowEdge(iter + 1) + length2, yVal + TMath::Sqrt(sys*sys), hist_p->GetBinLowEdge(iter + 2) -length2, yVal + TMath::Sqrt(sys*sys));

    l->DrawLine(hist_p->GetBinLowEdge(iter + 1) + length2, yVal + TMath::Sqrt(sys*sys) - length1, hist_p->GetBinLowEdge(iter + 1) + length2, yVal + TMath::Sqrt(sys*sys));

    l->DrawLine(hist_p->GetBinLowEdge(iter + 2) - length2, yVal + TMath::Sqrt(sys*sys) - length1, hist_p->GetBinLowEdge(iter + 2) - length2, yVal + TMath::Sqrt(sys*sys));
  }

  return;
}


void histToCanvas_TOP(TCanvas* plotCanvas_p, TH1F* histPP_p[7], TH1F* hist1_p[7], TH1F* hist2_p[7], TH1F* hist3_p[7], TH1F* hist4_p[7], const Int_t setNum, const std::string CNCR, const std::string projMult, const std::string Tight, const Bool_t montecarlo, const Bool_t isHighPtTrk, Bool_t isDoubleSub = false, Bool_t forceSys = false)
{
  const char* mcLabel[4] = {"PYTHIA", "PYT.+HYD.", "PYT.+HYD.", "(P+H)-P"};
  const char* dataLabel[4] = {"pp 5.3 pb^{-1}", "PbPb" /* 150 #mub^{-1}"*/, "PbPb", "PbPb-pp"};
  const char* overLabel[4];
  Float_t overCoord[4] = {.84, .76, .90, .82};
  if(strcmp(CNCR.c_str(), "") != 0 && strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0){
    for(Int_t coordIter = 0; coordIter < 4; coordIter++){
      overCoord[coordIter] += .04;
    }
  }
  for(Int_t iter = 0; iter < 4; iter++){
    if(montecarlo && !forceSys)
      overLabel[iter] = mcLabel[iter];
    else
      overLabel[iter] = dataLabel[iter];
  }

  plotCanvas_p->cd(1);

  TH1F* sysClonePP_p = (TH1F*)histPP_p[6]->Clone();

  makeHistForPtStack(histPP_p, 1, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
  histPP_p[0]->GetYaxis()->SetTitleOffset(2.2);
  histPP_p[1]->GetYaxis()->SetTitleOffset(2.2);
  drawFullStack(histPP_p, 0, 25, projMult.c_str(), 1, false, CNCR.c_str(), isDoubleSub);
  makePatch(CNCR.c_str());

  Float_t sysPP[10];
  Float_t sysPP5_1[10];
  Float_t sysPP1_2[10];

  if(!isDoubleSub && strcmp(CNCR.c_str(), "R_320") != 0 && strcmp(CNCR.c_str(), "RD_320") != 0 && strcmp(CNCR.c_str(), "RU_320") != 0 && strcmp(CNCR.c_str(), "R_340") != 0 && strcmp(CNCR.c_str(), "RD_340") != 0 && strcmp(CNCR.c_str(), "RU_340") != 0 && (strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0) && (!montecarlo || forceSys) && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 22 || setNum == 23 || setNum == 24 || setNum == 25 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31 || setNum == 32 || setNum == 33 || setNum == 34 || setNum == 35)){
    if(setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31){
      TFile* getSys_p = new TFile(Form("systDir/total_syst_pp_A%s_CaloResProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_pp_%sCaloResProjF", radString2[setNum]));
      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	sysPP[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;

      if(strcmp(CNCR.c_str(), "") != 0){
	TFile* getSys5_1_p = new TFile(Form("systDir/total_syst_pp_A%s_CaloResProj_5_1.root", CNCR.c_str()), "READ");
	TH1F* getSysHist5_1_p = (TH1F*)getSys5_1_p->Get(Form("syst_pp_%sCaloResProj5_1", radString2[setNum]));
	for(Int_t iter = 0; iter < getSysHist5_1_p->GetNbinsX(); iter++){
	  sysPP5_1[iter] = getSysHist5_1_p->GetBinContent(iter+1);
	}
	getSys5_1_p->Close();
	delete getSys5_1_p;
	
	
	TFile* getSys1_2_p = new TFile(Form("systDir/total_syst_pp_A%s_CaloResProj_1_2.root", CNCR.c_str()), "READ");
	TH1F* getSysHist1_2_p = (TH1F*)getSys1_2_p->Get(Form("syst_pp_%sCaloResProj1_2", radString2[setNum]));
	for(Int_t iter = 0; iter < getSysHist1_2_p->GetNbinsX(); iter++){
	  sysPP1_2[iter] = getSysHist1_2_p->GetBinContent(iter+1);
	}
	getSys1_2_p->Close();
	delete getSys1_2_p;
	
	for(Int_t iter = 0; iter < 10; iter++){
	  sysPP5_1[iter] += sysPP1_2[iter];
	  if(histPP_p[1]->GetBinContent(iter+1) < histPP_p[2]->GetBinContent(iter+1)) histPP_p[1]->SetBinContent(iter+1, histPP_p[2]->GetBinContent(iter+1));
	}

	//	makeSysError(sysPP5_1, histPP_p[1], true, CNCR);
	//	makeSysError(sysPP1_2, histPP_p[2], false, CNCR);
      }
    }
    else{
      TFile* getSys_p = new TFile(Form("systDir/total_syst_pp_A%s_CaloSwapProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_pp_%sCaloSwapProjF", radString2[setNum]));
      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	sysPP[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;
    }
       makeSysError(sysPP, sysClonePP_p, false, CNCR);

    
  }

  TLine* zeroLine_p;
  if(isAll(CNCR.c_str())) zeroLine_p = new TLine(0., 0., 2.0, 0.);
  else zeroLine_p = new TLine(0., 0., 0.5, 0.);

  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(1);
  zeroLine_p->DrawClone();

  std::cout << "A3" << std::endl;

  TLatex* label1_p = new TLatex();
  label1_p->SetNDC();
  label1_p->SetTextFont(43);
  label1_p->SetTextSizePixels(28);

  if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")){
    label1_p->DrawLatex(.28, overCoord[0], Form("%s", overLabel[0]));
    label1_p->DrawLatex(.28, overCoord[1], "CMS Preliminary");
  }
  else{
    label1_p->DrawLatex(.62, overCoord[0], Form("%s", overLabel[0]));
    label1_p->DrawLatex(.40, overCoord[1], "CMS Preliminary");
  }
  plotCanvas_p->cd(2);

  TH1F* sysClone1_p = (TH1F*)hist1_p[6]->Clone();

  makeHistForPtStack(hist1_p, 2, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
  drawFullStack(hist1_p, 0, 28, projMult.c_str(), 2, false, CNCR.c_str(), isDoubleSub);

  if(isAll(CNCR.c_str()) || isInc(CNCR.c_str()) || isU(CNCR.c_str()) || isD(CNCR.c_str())){
    if(!isDoubleSub && strcmp(projMult.c_str(), "Mult") != 0) histPP_p[6]->DrawCopy("SAME HIST C");
  }

  Float_t sys1[10];
  Float_t sys15_1[10];
  Float_t sys11_2[10];

  if(!isDoubleSub && strcmp(CNCR.c_str(), "R_320") != 0 && strcmp(CNCR.c_str(), "RD_320") != 0 && strcmp(CNCR.c_str(), "RU_320") != 0 && strcmp(CNCR.c_str(), "R_340") != 0 && strcmp(CNCR.c_str(), "RD_340") != 0 && strcmp(CNCR.c_str(), "RU_340") != 0 && (strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0) && (!montecarlo || forceSys) && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 22 || setNum == 23 || setNum == 24 || setNum == 25 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31 || setNum == 32 || setNum == 33 || setNum == 34 || setNum == 35)){
    if(setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31){
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloResProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      if(!strcmp(CNCR.c_str(), "")) getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_%sCaloResProjF_50100", radString2[setNum]));
      else getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_%sCaloResProjF_30100", radString2[setNum]));

      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	sys1[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;

      if(strcmp(CNCR.c_str(), "") != 0){
	TFile* getSys5_1_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloResProj_5_1.root", CNCR.c_str()), "READ");
	TH1F* getSysHist5_1_p;
	getSysHist5_1_p = (TH1F*)getSys5_1_p->Get(Form("syst_PbPb_%sCaloResProj5_1_30100", radString2[setNum]));

	for(Int_t iter = 0; iter < getSysHist5_1_p->GetNbinsX(); iter++){
	  sys15_1[iter] = getSysHist5_1_p->GetBinContent(iter+1);
	}
	getSys5_1_p->Close();
	delete getSys5_1_p;

	TFile* getSys1_2_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloResProj_1_2.root", CNCR.c_str()), "READ");
	TH1F* getSysHist1_2_p;
	getSysHist1_2_p = (TH1F*)getSys1_2_p->Get(Form("syst_PbPb_%sCaloResProj1_2_30100", radString2[setNum]));

	for(Int_t iter = 0; iter < getSysHist1_2_p->GetNbinsX(); iter++){
	  sys11_2[iter] = getSysHist1_2_p->GetBinContent(iter+1);
	}
	getSys1_2_p->Close();
	delete getSys1_2_p;

	for(Int_t iter = 0; iter < 10; iter++){
	  sys15_1[iter] += sys11_2[iter];
	  if(hist1_p[1]->GetBinContent(iter+1) < hist1_p[2]->GetBinContent(iter+1)) hist1_p[1]->SetBinContent(iter+1, hist1_p[2]->GetBinContent(iter+1));
	}

	//	makeSysError(sys15_1, hist1_p[1], true, CNCR);
	//	makeSysError(sys11_2, hist1_p[2], false, CNCR);
      }
    }
    else{
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloSwapProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      if(!strcmp(CNCR.c_str(), "")) getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_50100", radString2[setNum]));
      else getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_30100", radString2[setNum]));
    
      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	sys1[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;
    }

        makeSysError(sys1, sysClone1_p, false, CNCR);
  }

  zeroLine_p->DrawClone();


  if(!strcmp("", CNCR.c_str()) || !strcmp("C", CNCR.c_str()) || !strcmp("NC", CNCR.c_str())){
    label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[1]));
    label1_p->DrawLatex(.05, overCoord[1], "50-100%");
  }
  else{
    label1_p->DrawLatex(.50, overCoord[0], Form("%s", overLabel[1]));
    label1_p->DrawLatex(.50, overCoord[1], "30-100%");
  }

  if(isAll(CNCR.c_str()) && !strcmp(projMult.c_str(), "Proj")) label1_p->DrawLatex(.20, .05, "#sqrt{s_{NN}} = 2.76 TeV");
  else label1_p->DrawLatex(.05, .05, "#sqrt{s_{NN}} = 2.76 TeV");

  if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")) plotCanvas_p->cd(6);
  else plotCanvas_p->cd(4);

  label1_p->DrawLatex(.30, .92, "p_{T}^{trk} (|#eta|<2.4)");
  drawBin(CNCR.c_str());

  plotCanvas_p->cd(3);

  TH1F* sysClone2_p = (TH1F*)hist2_p[6]->Clone();

  makeHistForPtStack(hist2_p, 3, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
  drawFullStack(hist2_p, 0, 28, projMult.c_str(), 3, false, CNCR.c_str(), isDoubleSub);

  if(!isDoubleSub && isAll(CNCR.c_str()) || isInc(CNCR.c_str()) || isU(CNCR.c_str()) || isD(CNCR.c_str()) && strcmp(projMult.c_str(), "Mult") != 0) histPP_p[6]->DrawCopy("SAME HIST C");

  Float_t sys2[10];
  Float_t sys25_1[10];
  Float_t sys21_2[10];

  if(!isDoubleSub && strcmp(CNCR.c_str(), "R_320") != 0 && strcmp(CNCR.c_str(), "RD_320") != 0 && strcmp(CNCR.c_str(), "RU_320") != 0 && strcmp(CNCR.c_str(), "R_340") != 0 && strcmp(CNCR.c_str(), "RD_340") != 0 && strcmp(CNCR.c_str(), "RU_340") != 0 && (strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0) && (!montecarlo || forceSys) && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 22 || setNum == 23 || setNum == 24 || setNum == 25 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31 || setNum == 32 || setNum == 33 || setNum == 34 || setNum == 35)){
    if(setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31){
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloResProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      if(!strcmp(CNCR.c_str(), "")) getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_%sCaloResProjF_3050", radString2[setNum]));
      else getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_%sCaloResProjF_030", radString2[setNum]));
      
      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	sys2[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;

      if(strcmp(CNCR.c_str(), "") != 0){
	TFile* getSys5_1_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloResProj_5_1.root", CNCR.c_str()), "READ");
	TH1F* getSysHist5_1_p;
	getSysHist5_1_p = (TH1F*)getSys5_1_p->Get(Form("syst_PbPb_%sCaloResProj5_1_030", radString2[setNum]));
	
	for(Int_t iter = 0; iter < getSysHist5_1_p->GetNbinsX(); iter++){
	  sys25_1[iter] = getSysHist5_1_p->GetBinContent(iter+1);
	}
	getSys5_1_p->Close();
	delete getSys5_1_p;
	
	TFile* getSys1_2_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloResProj_1_2.root", CNCR.c_str()), "READ");
	TH1F* getSysHist1_2_p;
	getSysHist1_2_p = (TH1F*)getSys1_2_p->Get(Form("syst_PbPb_%sCaloResProj1_2_030", radString2[setNum]));
	
	for(Int_t iter = 0; iter < getSysHist1_2_p->GetNbinsX(); iter++){
	  sys21_2[iter] = getSysHist1_2_p->GetBinContent(iter+1);
	}
	getSys1_2_p->Close();
	delete getSys1_2_p;
	
	for(Int_t iter = 0; iter< 10; iter++){
	  sys25_1[iter] += sys21_2[iter];
	  if(hist2_p[1]->GetBinContent(iter+1) < hist2_p[2]->GetBinContent(iter+1)) hist2_p[1]->SetBinContent(iter+1, hist2_p[2]->GetBinContent(iter+1));
	}

	//	makeSysError(sys25_1, hist2_p[1], true, CNCR);
	//	makeSysError(sys21_2, hist2_p[2], false, CNCR);
      }
    }
    else{
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloSwapProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      if(!strcmp(CNCR.c_str(), "")) getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_3050", radString2[setNum]));
      else getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_030", radString2[setNum]));
      
      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	sys2[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;
    }
    
        makeSysError(sys2, sysClone2_p, false, CNCR);
  }

  zeroLine_p->DrawClone();

  if(!strcmp("", CNCR.c_str()) || !strcmp("C", CNCR.c_str()) || !strcmp("NC", CNCR.c_str())){
    label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[2]));
    label1_p->DrawLatex(.05, overCoord[1], "30-50%");
  }
  else{
    label1_p->DrawLatex(.50, overCoord[1], "0-30%");
    label1_p->DrawLatex(.50, overCoord[0], Form("%s", overLabel[2]));
  }

  if(isAll(CNCR.c_str())){
    label1_p->DrawLatex(.10, .15, Form("anti-k_{T} R=%s",radString[setNum]));
    label1_p->DrawLatex(.10, .05, Form("%s", puVsString[setNum]));
  }
  else{
    label1_p->DrawLatex(.05, .15, Form("anti-k_{T}"));
    label1_p->DrawLatex(.05, .05, Form("%s R=%s", puVsString[setNum], radString[setNum]));
  }
  if(isHighPtTrk) label1_p->DrawLatex(.05, .05, "In jet p_{T}^{trk}>12 GeV/c");

  if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")){
    plotCanvas_p->cd(4);

    TH1F* sysClone3_p = (TH1F*)hist3_p[6]->Clone();

    makeHistForPtStack(hist3_p, 4, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
    drawFullStack(hist3_p, 0, 28, projMult.c_str(), 4, false, CNCR.c_str(), isDoubleSub);

    Float_t sys3[10];

  if(!isDoubleSub && (strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0) && (!montecarlo || forceSys) && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 22 || setNum == 23 || setNum == 24 || setNum == 25 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31 || setNum == 32 || setNum == 33 || setNum == 34 || setNum == 35)){
    if(setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31){
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloResProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_%sCaloResProjF_1030", radString2[setNum]));

      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	sys3[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;
    }
    else{
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloSwapProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_1030", radString2[setNum]));
       
      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	sys3[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;
    }

    makeSysError(sys3, sysClone3_p, false, CNCR);
  }
   
 //    if(setNum == 16 && !isDoubleSub && (!montecarlo || forceSys) && !strcmp("", CNCR.c_str()) && !strcmp(Tight.c_str(), "")) makeSysError(sysA1030, hist3_p[6]);

    zeroLine_p->DrawClone();
    label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[2]));
    label1_p->DrawLatex(.05, overCoord[1], "10-30%");

    plotCanvas_p->cd(5);


    TH1F* sysClone4_p = (TH1F*)hist4_p[6]->Clone();


    std::cout << "A9" <<std::endl;

    makeHistForPtStack(hist4_p, 5, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
    drawFullStack(hist4_p, 0, 28, projMult.c_str(), 5, false, CNCR.c_str(), isDoubleSub);

    std::cout << "A10" <<std::endl;

    Float_t sys4[10];

  if(!isDoubleSub && (strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0)  && (!montecarlo || forceSys) && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 22 || setNum == 23 || setNum == 24 || setNum == 25 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31 || setNum == 32 || setNum == 33 || setNum == 34 || setNum == 35)){
    if(setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31){
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloResProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_%sCaloResProjF_010", radString2[setNum]));

      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	sys4[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;
    }
    else{
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_A%s_CaloSwapProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_010", radString2[setNum]));
       
      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	sys4[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;
    }

    makeSysError(sys4, sysClone4_p, false, CNCR);
  }

    //    if(setNum == 16 && !isDoubleSub && (!montecarlo || forceSys) && !strcmp("", CNCR.c_str()) && !strcmp(Tight.c_str(), "")) makeSysError(sysA010, hist4_p[6]);

    zeroLine_p->DrawClone();
    label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[2]));
    label1_p->DrawLatex(.05, overCoord[1], "0-10%");
  }
  std::cout << "EndTop" << std::endl;
  return;
}



void histToCanvasR(TCanvas* plotCanvas_p, const Int_t canvPos,  TH1F* hist1_p[7], TH1F* hist2_p[7], TH1F* hist3_p[7], TH1F* hist4_p[7], const Int_t setNum, const std::string CNCR, const std::string projMult, Bool_t montecarlo, const std::string cent = "", Bool_t isSub = false)
{
  plotCanvas_p->cd(canvPos);

  const std::string Tight = "";
  const Bool_t isDoubleSub = false;

  std::string resSwap = "Res";
  if(setNum == 22 || setNum == 31) resSwap = "Swap";

  std::string ppPbPbFile = "pp";
  if(canvPos == 5) ppPbPbFile = "PbPb";
  else if(canvPos == 9) ppPbPbFile = "PbPb_pp";

  std::string ppPbPbHist = "pp";
  if(canvPos == 5){
    if(setNum == 22 || setNum == 31) ppPbPbHist = "PbPb_pp";
    else ppPbPbHist = "PbPb";
  }
  else if(canvPos == 9) ppPbPbHist = "PbPb_pp";

  std::string centHist = "";
  if(strcmp("", cent.c_str()) != 0) centHist = Form("_%s", cent.c_str());

  Bool_t doSyst = true;
  if(!strcmp(CNCR.c_str(), "") && (!strcmp("030", cent.c_str()) || !strcmp("30100", cent.c_str()))) doSyst = false;
  else if(setNum != 16 && setNum != 22 && setNum != 28 && setNum != 32) doSyst = false;
  else if(montecarlo) doSyst = false;
  else if(!strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")) doSyst = false;
  else if(!strcmp(CNCR.c_str(), "R_320") || !strcmp(CNCR.c_str(), "RD_320") ||  !strcmp(CNCR.c_str(), "RU_320")) doSyst = false;
  else if(!strcmp(CNCR.c_str(), "R_340") || !strcmp(CNCR.c_str(), "RD_340") ||  !strcmp(CNCR.c_str(), "RU_340")) doSyst = false;

  Bool_t isBottom = true;
  if(canvPos == 9) isBottom = true;

  Int_t style = 25;
  if(canvPos == 5) style = 28;
  else if(canvPos == 9) style = 24;

  TH1F* sysClone1_p = (TH1F*)hist1_p[6]->Clone();

  makeHistForPtStack(hist1_p, 1, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
  drawFullStack(hist1_p, 0, style, projMult.c_str(), canvPos, isSub, CNCR.c_str(), isDoubleSub, true);

  TLine* zeroLine_p;
  if(isAll(CNCR.c_str())) zeroLine_p = new TLine(0., 0., 2.0, 0.);
  else zeroLine_p = new TLine(0., 0., 0.5, 0.);

  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(1);
  zeroLine_p->DrawClone();


  Float_t sys1[10];
  Float_t sys15_1[10];
  Float_t sys11_2[10];

  if(doSyst){
    TFile* getSys_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
    TH1F* getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_%s_%sCalo%sProjF%s", ppPbPbHist.c_str(), radString2[setNum], resSwap.c_str(), centHist.c_str()));
    for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
      sys1[iter] = getSysHist_p->GetBinContent(iter+1);
    }
    getSys_p->Close();
    delete getSys_p;
    makeSysError(sys1, sysClone1_p, isBottom, CNCR);

    if(strcmp(CNCR.c_str(), "") != 0 && !strcmp("Res", resSwap.c_str())){
      TFile* getSys5_1_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj_5_1.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
      TH1F* getSysHist5_1_p = (TH1F*)getSys5_1_p->Get(Form("syst_%s_%sCalo%sProj5_1%s", ppPbPbHist.c_str(), radString2[setNum], resSwap.c_str(), centHist.c_str()));
      for(Int_t iter = 0; iter < getSysHist5_1_p->GetNbinsX(); iter++){
	sys15_1[iter] = getSysHist5_1_p->GetBinContent(iter+1);
      }
      getSys5_1_p->Close();
      delete getSys5_1_p;

      TFile* getSys1_2_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj_1_2.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
      TH1F* getSysHist1_2_p = (TH1F*)getSys1_2_p->Get(Form("syst_%s_%sCalo%sProj1_2%s", ppPbPbHist.c_str(), radString2[setNum], resSwap.c_str(), centHist.c_str()));
      for(Int_t iter = 0; iter < getSysHist1_2_p->GetNbinsX(); iter++){
	sys11_2[iter] = getSysHist1_2_p->GetBinContent(iter+1);
      }
      getSys1_2_p->Close();
      delete getSys1_2_p;

      for(Int_t iter = 0; iter < 10; iter++){
	sys15_1[iter] += sys11_2[iter];

	if(hist1_p[1]->GetBinContent(iter+1) < hist1_p[2]->GetBinContent(iter+1)) hist1_p[1]->SetBinContent(iter+1, hist1_p[2]->GetBinContent(iter+1));
      }

      //     makeSysError(sys15_1, hist1_p[1], isBottom, CNCR);
    }
  }

  plotCanvas_p->cd(canvPos+1);

  TH1F* sysClone2_p = (TH1F*)hist2_p[6]->Clone();

  makeHistForPtStack(hist2_p, 2, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
  drawFullStack(hist2_p, 0, style, projMult.c_str(), canvPos + 1, isSub, CNCR.c_str(), isDoubleSub, true);

  Float_t sys2[10];
  Float_t sys25_1[10];
  Float_t sys21_2[10];

  if(doSyst){
    TFile* getSys_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
    TH1F* getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_%s_%sCalo%sProjF%s", ppPbPbHist.c_str(), radString2[setNum+1], resSwap.c_str(), centHist.c_str()));
    for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
      sys2[iter] = getSysHist_p->GetBinContent(iter+1);
    }
    getSys_p->Close();
    delete getSys_p;
    makeSysError(sys2, sysClone2_p, isBottom, CNCR);
    
    if(strcmp("", CNCR.c_str()) != 0 && !strcmp("Res", resSwap.c_str())){
      TFile* getSys5_1_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj_5_1.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
      TH1F* getSysHist5_1_p = (TH1F*)getSys5_1_p->Get(Form("syst_%s_%sCalo%sProj5_1%s", ppPbPbHist.c_str(), radString2[setNum+1], resSwap.c_str(), centHist.c_str()));
      for(Int_t iter = 0; iter < getSysHist5_1_p->GetNbinsX(); iter++){
	sys25_1[iter] = getSysHist5_1_p->GetBinContent(iter+1);
      }
      getSys5_1_p->Close();
      delete getSys5_1_p;

      TFile* getSys1_2_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj_1_2.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
      TH1F* getSysHist1_2_p = (TH1F*)getSys1_2_p->Get(Form("syst_%s_%sCalo%sProj1_2%s", ppPbPbHist.c_str(), radString2[setNum+1], resSwap.c_str(), centHist.c_str()));
      for(Int_t iter = 0; iter < getSysHist1_2_p->GetNbinsX(); iter++){
	sys21_2[iter] = getSysHist1_2_p->GetBinContent(iter+1);
      }
      getSys1_2_p->Close();
      delete getSys1_2_p;

      for(Int_t iter = 0; iter < 10; iter++){
	sys25_1[iter] += sys21_2[iter];

	if(hist2_p[1]->GetBinContent(iter+1) < hist2_p[2]->GetBinContent(iter+1)) hist2_p[1]->SetBinContent(iter+1, hist2_p[2]->GetBinContent(iter+1));
      }

      //      makeSysError(sys25_1, hist2_p[1], isBottom, CNCR);
    }
  }

  zeroLine_p->DrawClone();

  plotCanvas_p->cd(canvPos+2);

  TH1F* sysClone3_p = (TH1F*)hist3_p[6]->Clone();

  makeHistForPtStack(hist3_p, 3, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
  drawFullStack(hist3_p, 0, style, projMult.c_str(), canvPos + 2, isSub, CNCR.c_str(), isDoubleSub, true);

  Float_t sys3[10];
  Float_t sys35_1[10];
  Float_t sys31_2[10];

  if(doSyst){
    TFile* getSys_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
    TH1F* getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_%s_%sCalo%sProjF%s", ppPbPbHist.c_str(), radString2[setNum+2], resSwap.c_str(), centHist.c_str()));
    for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
      sys3[iter] = getSysHist_p->GetBinContent(iter+1);
    }
    getSys_p->Close();
    delete getSys_p;
    makeSysError(sys3, sysClone3_p, isBottom, CNCR);

    if(strcmp("", CNCR.c_str()) != 0 && !strcmp("Res", resSwap.c_str())){
      TFile* getSys5_1_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj_5_1.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
      TH1F* getSysHist5_1_p = (TH1F*)getSys5_1_p->Get(Form("syst_%s_%sCalo%sProj5_1%s", ppPbPbHist.c_str(), radString2[setNum+2], resSwap.c_str(), centHist.c_str()));
      for(Int_t iter = 0; iter < getSysHist5_1_p->GetNbinsX(); iter++){
	sys35_1[iter] = getSysHist5_1_p->GetBinContent(iter+1);
      }
      getSys5_1_p->Close();
      delete getSys5_1_p;

      TFile* getSys1_2_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj_1_2.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
      TH1F* getSysHist1_2_p = (TH1F*)getSys1_2_p->Get(Form("syst_%s_%sCalo%sProj1_2%s", ppPbPbHist.c_str(), radString2[setNum+2], resSwap.c_str(), centHist.c_str()));
      for(Int_t iter = 0; iter < getSysHist1_2_p->GetNbinsX(); iter++){
	sys31_2[iter] = getSysHist1_2_p->GetBinContent(iter+1);
      }
      getSys1_2_p->Close();
      delete getSys1_2_p;

      for(Int_t iter = 0; iter < 10; iter++){
	sys35_1[iter] += sys31_2[iter];

	if(hist3_p[1]->GetBinContent(iter+1) < hist3_p[2]->GetBinContent(iter+1)) hist3_p[1]->SetBinContent(iter+1, hist3_p[2]->GetBinContent(iter+1));
      }

      //      makeSysError(sys35_1, hist3_p[1], isBottom, CNCR);

    }
  }

  zeroLine_p->DrawClone();

  plotCanvas_p->cd(canvPos+3);

  TH1F* sysClone4_p = (TH1F*)hist4_p[6]->Clone();
  
  makeHistForPtStack(hist4_p, 4, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
  drawFullStack(hist4_p, 0, style, projMult.c_str(), canvPos + 3, isSub, CNCR.c_str(), isDoubleSub, true);

  Float_t sys4[10];
  Float_t sys45_1[10];
  Float_t sys41_2[10];

  if(doSyst){
    TFile* getSys_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
    TH1F* getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_%s_%sCalo%sProjF%s", ppPbPbHist.c_str(), radString2[setNum+3], resSwap.c_str(), centHist.c_str()));
    for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
      sys4[iter] = getSysHist_p->GetBinContent(iter+1);
    }
    getSys_p->Close();
    delete getSys_p;
    makeSysError(sys4, sysClone4_p, isBottom, CNCR);

    if(strcmp("", CNCR.c_str()) != 0 && !strcmp("Res", resSwap.c_str())){
      TFile* getSys5_1_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj_5_1.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
      TH1F* getSysHist5_1_p = (TH1F*)getSys5_1_p->Get(Form("syst_%s_%sCalo%sProj5_1%s", ppPbPbHist.c_str(), radString2[setNum+3], resSwap.c_str(), centHist.c_str()));
      for(Int_t iter = 0; iter < getSysHist5_1_p->GetNbinsX(); iter++){
	sys45_1[iter] = getSysHist5_1_p->GetBinContent(iter+1);
      }
      getSys5_1_p->Close();
      delete getSys5_1_p;

      TFile* getSys1_2_p = new TFile(Form("systDir/total_syst_%s_A%s_Calo%sProj_1_2.root", ppPbPbFile.c_str(), CNCR.c_str(), resSwap.c_str()), "READ");
      TH1F* getSysHist1_2_p = (TH1F*)getSys1_2_p->Get(Form("syst_%s_%sCalo%sProj1_2%s", ppPbPbHist.c_str(), radString2[setNum+3], resSwap.c_str(), centHist.c_str()));
      for(Int_t iter = 0; iter < getSysHist1_2_p->GetNbinsX(); iter++){
	sys41_2[iter] = getSysHist1_2_p->GetBinContent(iter+1);
      }
      getSys1_2_p->Close();
      delete getSys1_2_p;

      for(Int_t iter = 0; iter < 10; iter++){
	sys45_1[iter] += sys41_2[iter];

	if(hist4_p[1]->GetBinContent(iter+1) < hist4_p[2]->GetBinContent(iter+1)) hist4_p[1]->SetBinContent(iter+1, hist4_p[2]->GetBinContent(iter+1));
      }

      //      makeSysError(sys45_1, hist4_p[1], isBottom, CNCR);
    }
  }

  zeroLine_p->DrawClone();
    
  return;
}


void histToCanvas_BOTTOM(TCanvas* plotCanvas_p, TH1F* histPP_p[7], TH1F* hist1_p[7], TH1F* hist2_p[7], TH1F* hist3_p[7], TH1F* hist4_p[7], const Int_t setNum, const std::string CNCR, const std::string projMult, const std::string Tight, const Bool_t montecarlo, const Bool_t isHighPtTrk, Bool_t isDoubleSub = false, Bool_t forceSys = false)
{
  const char* mcLabel[4] = {"PYTHIA", "PYT.+HYD.", "PYT.+HYD.", "(P+H)-P"};
  const char* dataLabel[4] = {"pp 5.3 pb^{-1}", "PbPb" /* 150 #mub^{-1}"*/, "PbPb", "PbPb-pp"};
  const char* overLabel[4];
  Float_t overCoord[4] = {.84, .76, .90, .82};
  if(strcmp(CNCR.c_str(), "") != 0 && strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0){
    for(Int_t coordIter = 0; coordIter < 4; coordIter++){
      overCoord[coordIter] += .04;
    }
  }
  for(Int_t iter = 0; iter < 4; iter++){
    if(montecarlo && !forceSys)
      overLabel[iter] = mcLabel[iter];
    else
      overLabel[iter] = dataLabel[iter];
  }

  TLatex* label1_p = new TLatex();
  label1_p->SetNDC();
  label1_p->SetTextFont(43);
  label1_p->SetTextSizePixels(28);

  TLine* zeroLine_p;
  if(isAll(CNCR.c_str())) zeroLine_p = new TLine(0., 0., 2.0, 0.);
  else zeroLine_p = new TLine(0., 0., 0.5, 0.);

  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(1);
  zeroLine_p->DrawClone();



  Int_t panels;
  Int_t ppStart;
  const char* ppChar[2];

  if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")){
    panels = 10;
    ppStart = 7;
    ppChar[0] = "50-100%";
    ppChar[1] = "30-50%";
  }
  else{
    panels = 6;
    ppStart = 5;
    ppChar[0] = "30-100%";
    ppChar[1] = "0-30%";
  }

  if(strcmp(CNCR.c_str(), "") != 0 && strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0){
    for(Int_t coordIter = 0; coordIter < 4; coordIter++){
      overCoord[coordIter] -= .04;
    }
  }

  plotCanvas_p->cd(ppStart);

  TH1F* sysClone1_p = (TH1F*)hist1_p[6]->Clone();

  makeHistForPtStack(hist1_p, ppStart, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
  drawFullStack(hist1_p, 0, 24, projMult.c_str(), ppStart, true, CNCR.c_str(), isDoubleSub);
  drawNum(CNCR.c_str());

  histPP_p[5]->SetMarkerStyle(25);
  histPP_p[5]->SetLineStyle(2);

  if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")){
    label1_p->DrawLatex(.22, overCoord[2], overLabel[3]);
    label1_p->DrawLatex(.22, overCoord[3], ppChar[0]);
  }
  else{
    label1_p->DrawLatex(.60, overCoord[2], overLabel[3]);
    label1_p->DrawLatex(.60, overCoord[3], ppChar[0]);
  }

  if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")){
    label1_p->DrawLatex(.24, .38, "p_{T,1}>120 GeV/c");
    label1_p->DrawLatex(.24, .28, "p_{T,2}>50 GeV/c");
  }
  else{
    label1_p->DrawLatex(.44, .38, "p_{T,1}>120 GeV/c");
    label1_p->DrawLatex(.44, .28, "p_{T,2}>50 GeV/c");
  }



  Float_t sys1[10];
  Float_t sys15_1[10];
  Float_t sys11_2[10];

  if(!isDoubleSub && strcmp(CNCR.c_str(), "R_320") != 0 && strcmp(CNCR.c_str(), "RD_320") != 0 && strcmp(CNCR.c_str(), "RU_320") != 0  && strcmp(CNCR.c_str(), "R_340") != 0 && strcmp(CNCR.c_str(), "RD_340") != 0 && strcmp(CNCR.c_str(), "RU_340") != 0 && (strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0) && (!montecarlo || forceSys) && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum== 22 || setNum == 23 || setNum == 24 || setNum == 25 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31 || setNum == 32 || setNum == 33 || setNum == 34 || setNum == 35)){
    if(setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31){
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloResProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      if(!strcmp(CNCR.c_str(), "")) getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloResProjF_50100", radString2[setNum]));
      else getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloResProjF_30100", radString2[setNum]));
      
      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
        sys1[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;

      if(strcmp(CNCR.c_str(), "") != 0){
	TFile* getSys5_1_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloResProj_5_1.root", CNCR.c_str()), "READ");
	TH1F* getSysHist5_1_p;
	getSysHist5_1_p = (TH1F*)getSys5_1_p->Get(Form("syst_PbPb_pp_%sCaloResProj5_1_30100", radString2[setNum]));
	
	for(Int_t iter = 0; iter < getSysHist5_1_p->GetNbinsX(); iter++){
	  sys15_1[iter] = getSysHist5_1_p->GetBinContent(iter+1);
	}
	getSys5_1_p->Close();
	delete getSys5_1_p;

	TFile* getSys1_2_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloResProj_1_2.root", CNCR.c_str()), "READ");
	TH1F* getSysHist1_2_p;
	getSysHist1_2_p = (TH1F*)getSys1_2_p->Get(Form("syst_PbPb_pp_%sCaloResProj1_2_30100", radString2[setNum]));
	
	for(Int_t iter = 0; iter < getSysHist1_2_p->GetNbinsX(); iter++){
	  sys11_2[iter] = getSysHist1_2_p->GetBinContent(iter+1);
	}
	getSys1_2_p->Close();
	delete getSys1_2_p;

	for(Int_t iter = 0; iter < 10; iter++){
	  sys15_1[iter] += sys11_2[iter];
	  if(hist1_p[1]->GetBinContent(iter+1) < hist1_p[2]->GetBinContent(iter+1)) hist1_p[1]->SetBinContent(iter+1, hist1_p[2]->GetBinContent(iter+1));
	}
	//	makeSysError(sys15_1, hist1_p[1], true, CNCR);
      }
    }
    else{
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloSwapProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      if(!strcmp(CNCR.c_str(), "")) getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_50100", radString2[setNum]));
      else getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_30100", radString2[setNum]));

      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
        sys1[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;
    }
    
        makeSysError(sys1, sysClone1_p, true, CNCR);
  }

  zeroLine_p->DrawClone();


  plotCanvas_p->cd(ppStart+1);
  TH1F* sysClone2_p = (TH1F*)hist2_p[6]->Clone();

  makeHistForPtStack(hist2_p, ppStart+1, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
  drawFullStack(hist2_p, 0, 24, projMult.c_str(), ppStart+1, true, CNCR.c_str(), isDoubleSub);

  if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")){
    label1_p->DrawLatex(.05, overCoord[2], overLabel[3]);
    label1_p->DrawLatex(.05, overCoord[3], ppChar[1]);
  }
  else{
    label1_p->DrawLatex(.50, overCoord[2], overLabel[3]);
    label1_p->DrawLatex(.60, overCoord[3], ppChar[1]);
  }

  if(!strcmp(CNCR.c_str(), "C0") || !strcmp(CNCR.c_str(), "C1") || !strcmp(CNCR.c_str(), "C2") || !strcmp(CNCR.c_str(), "C3") || isAll(CNCR.c_str())){
    if(!strcmp(CNCR.c_str(), "R2") || !strcmp(CNCR.c_str(), "R2U") || !strcmp(CNCR.c_str(), "R2D")) label1_p->DrawLatex(.45, .38, "|#eta|_{1},|#eta|_{2}<0.8");
    else if(!strcmp(CNCR.c_str(), "RMID") || !strcmp(CNCR.c_str(), "RMIDU") || !strcmp(CNCR.c_str(), "RMIDD")) label1_p->DrawLatex(.45, .38, "|#eta|_{1},|#eta|_{2}<1.6");
    else label1_p->DrawLatex(.45, .38, "|#eta|_{1},|#eta|_{2}<0.6");
    label1_p->DrawLatex(.45, .28, "#Delta#phi_{1,2}>5#pi/6");
  }
  else{
    if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")){
      label1_p->DrawLatex(.05, .38, "|#eta|_{1},|#eta|_{2}<1.6");
      label1_p->DrawLatex(.05, .28, "#Delta#phi_{1,2}>5#pi/6");
    }
    else{
      label1_p->DrawLatex(.45, .38, "|#eta|_{1},|#eta|_{2}<1.6");
      label1_p->DrawLatex(.45, .28, "#Delta#phi_{1,2}>5#pi/6");
    }
  }

  Float_t sys2[10];
  Float_t sys25_1[10];
  Float_t sys21_2[10];

  if(!isDoubleSub && strcmp(CNCR.c_str(), "R_320") != 0 && strcmp(CNCR.c_str(), "RD_320") != 0 && strcmp(CNCR.c_str(), "RU_320") != 0  && strcmp(CNCR.c_str(), "R_340") != 0 && strcmp(CNCR.c_str(), "RD_340") != 0 && strcmp(CNCR.c_str(), "RU_340") != 0 && (strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0) && (!montecarlo || forceSys) && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum== 22 || setNum == 23 || setNum == 24 || setNum == 25 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31 || setNum == 32 || setNum == 33 || setNum == 34 || setNum == 35)){
    if(setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31){
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloResProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      if(!strcmp(CNCR.c_str(), "")) getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloResProjF_3050", radString2[setNum]));
      else getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloResProjF_030", radString2[setNum]));

      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
        sys2[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;

      if(strcmp(CNCR.c_str(), "") != 0){
	TFile* getSys5_1_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloResProj_5_1.root", CNCR.c_str()), "READ");
	TH1F* getSysHist5_1_p;
	getSysHist5_1_p = (TH1F*)getSys5_1_p->Get(Form("syst_PbPb_pp_%sCaloResProj5_1_030", radString2[setNum]));

	for(Int_t iter = 0; iter < getSysHist5_1_p->GetNbinsX(); iter++){
	  sys25_1[iter] = getSysHist5_1_p->GetBinContent(iter+1);
	}
	getSys5_1_p->Close();
	delete getSys5_1_p;

	TFile* getSys1_2_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloResProj_1_2.root", CNCR.c_str()), "READ");
	TH1F* getSysHist1_2_p;
	getSysHist1_2_p = (TH1F*)getSys1_2_p->Get(Form("syst_PbPb_pp_%sCaloResProj1_2_030", radString2[setNum]));

	for(Int_t iter = 0; iter < getSysHist1_2_p->GetNbinsX(); iter++){
	  sys21_2[iter] = getSysHist1_2_p->GetBinContent(iter+1);
	}
	getSys1_2_p->Close();
	delete getSys1_2_p;

	for(Int_t iter = 0; iter < 10; iter++){
	  sys25_1[iter] += sys21_2[iter];
	  if(hist2_p[1]->GetBinContent(iter+1) < hist2_p[2]->GetBinContent(iter+1)) hist2_p[1]->SetBinContent(iter+1, hist2_p[2]->GetBinContent(iter+1));
	}

	//	makeSysError(sys25_1, hist2_p[1], true, CNCR);
      }
    }
    else{
      TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloSwapProj.root", CNCR.c_str()), "READ");
      TH1F* getSysHist_p;
      if(!strcmp(CNCR.c_str(), "")) getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_3050", radString2[setNum]));
      else getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_030", radString2[setNum]));
      
      for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
        sys2[iter] = getSysHist_p->GetBinContent(iter+1);
      }
      getSys_p->Close();
      delete getSys_p;
    }
    
       makeSysError(sys2, sysClone2_p, true, CNCR);
  }

  zeroLine_p->DrawClone();

  if(!strcmp(CNCR.c_str(), "") || !strcmp(CNCR.c_str(), "C") || !strcmp(CNCR.c_str(), "NC")){
    plotCanvas_p->cd(9);

    TH1F* sysClone3_p = (TH1F*)hist3_p[6]->Clone();

    makeHistForPtStack(hist3_p, 9, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
    drawFullStack(hist3_p, 0, 24, projMult.c_str(), 9, true, CNCR.c_str(), isDoubleSub);
    label1_p->DrawLatex(.05, overCoord[2], overLabel[3]);
    label1_p->DrawLatex(.05, overCoord[3], "10-30%");
    //    label1_p->DrawLatex(.05, .38, "anti-k_{T} Calo R = 0.3");                                                                         
    Float_t sys3[10];

    if(!isDoubleSub && strcmp(CNCR.c_str(), "R_320") != 0 && strcmp(CNCR.c_str(), "RD_320") != 0 && strcmp(CNCR.c_str(), "RU_320") != 0  && strcmp(CNCR.c_str(), "R_340") != 0 && strcmp(CNCR.c_str(), "RD_340") != 0 && strcmp(CNCR.c_str(), "RU_340") != 0 && (strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0) && (!montecarlo || forceSys) && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum== 22 || setNum == 23 || setNum == 24 || setNum == 25 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31 || setNum == 32 || setNum == 33 || setNum == 34 || setNum == 35)){
      if(setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31){
	TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloResProj.root", CNCR.c_str()), "READ");
	TH1F* getSysHist_p;
	getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloResProjF_1030", radString2[setNum]));

	for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	  sys3[iter] = getSysHist_p->GetBinContent(iter+1);
	}
	getSys_p->Close();
	delete getSys_p;
      }
      else{
	TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloSwapProj.root", CNCR.c_str()), "READ");
	TH1F* getSysHist_p;
	getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_1030", radString2[setNum]));

	for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	  sys3[iter] = getSysHist_p->GetBinContent(iter+1);
	}
	getSys_p->Close();
	delete getSys_p;
      }

      makeSysError(sys3, sysClone3_p, true, CNCR);
    }


    zeroLine_p->DrawClone();

    plotCanvas_p->cd(10);

    TH1F* sysClone4_p = (TH1F*)hist4_p[6]->Clone();

    makeHistForPtStack(hist4_p, 10, Tight.c_str(), CNCR.c_str(), projMult.c_str(), isDoubleSub);
    drawFullStack(hist4_p, 0, 24, projMult.c_str(), 10, true, CNCR.c_str(), isDoubleSub);
    label1_p->DrawLatex(.05, overCoord[2], overLabel[3]);
    label1_p->DrawLatex(.05, overCoord[3], "0-10%");

    Float_t sys4[10];

    if(!isDoubleSub && strcmp(CNCR.c_str(), "R_320") != 0 && strcmp(CNCR.c_str(), "RD_320") != 0 && strcmp(CNCR.c_str(), "RU_320") != 0 && strcmp(CNCR.c_str(), "R_340") != 0 && strcmp(CNCR.c_str(), "RD_340") != 0 && strcmp(CNCR.c_str(), "RU_340") != 0 && (strcmp(CNCR.c_str(), "C") != 0 && strcmp(CNCR.c_str(), "NC") != 0) && (!montecarlo || forceSys) && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum== 22 || setNum == 23 || setNum == 24 || setNum == 25 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31 || setNum == 32 || setNum == 33 || setNum == 34 || setNum == 35)){
      if(setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31){
	TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloResProj.root", CNCR.c_str()), "READ");
	TH1F* getSysHist_p;
	getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloResProjF_010", radString2[setNum]));

	for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	  sys4[iter] = getSysHist_p->GetBinContent(iter+1);
	}
	getSys_p->Close();
	delete getSys_p;
      }
      else{
	TFile* getSys_p = new TFile(Form("systDir/total_syst_PbPb_pp_A%s_CaloSwapProj.root", CNCR.c_str()), "READ");
	TH1F* getSysHist_p;
	getSysHist_p = (TH1F*)getSys_p->Get(Form("syst_PbPb_pp_%sCaloSwapProjF_010", radString2[setNum]));

	for(Int_t iter = 0; iter < getSysHist_p->GetNbinsX(); iter++){
	  sys4[iter] = getSysHist_p->GetBinContent(iter+1);
	}
	getSys_p->Close();
	delete getSys_p;
      }

      makeSysError(sys4, sysClone4_p, true, CNCR);
    }

    zeroLine_p->DrawClone();
  }

  for(Int_t panelIter = 0; panelIter < panels; panelIter++){
    plotCanvas_p->cd(panelIter+1)->RedrawAxis();
  }

  return;
}


void makeMultStack(const char* filePbPbName, const char* fileTagPbPb, const char* outName, Int_t setNum, Bool_t montecarlo = false, const char* filePPName = "", const char* fileTagPP = "", const char* Tight = "")
{
  TFile* histPbPbFile_p = new TFile(filePbPbName, "READ");
  TFile* histPPFile_p;
  if(strcmp(fileTagPP, "") != 0)
    histPPFile_p = new TFile(filePPName, "READ");

  const char* cent[4] = {"50100", "3050", "1030", "010"};

  TH1F* histPbPb_p[4];
  TH1F* histPP_p;

  if(strcmp(fileTagPP, "") != 0){
    histPP_p = (TH1F*)histPPFile_p->Get(Form("r%sMultA%sFCorr_PP_%s_h", algType[setNum], Tight, fileTagPP));
    std::cout << Form("r%sMultA%sFCorr_PP_%s_h", algType[setNum], Tight, fileTagPP) << std::endl;
    histPP_p->SetYTitle("Hemisphere #Delta_{mult}");
    histPP_p->GetYaxis()->SetTitleOffset(2.2);
  }

  TCanvas* profPanel_p = new TCanvas(Form("r%sMultA%sFCorr_Stack_%s_c", algType[setNum], Tight, fileTagPbPb), Form("r%sMultA%s_Stack_%s_c", algType[setNum], Tight, fileTagPbPb), 4*300, 2*350);
  profPanel_p->Divide(4, 2, 0.0, 0.0);
  std::cout << "FourPanel Init" << std::endl;

  if(strcmp(fileTagPP, "") != 0){
    histPP_p->SetFillColor(kBlue);
    histPP_p->SetLineColor(kBlue);
    histPP_p->SetMarkerColor(kBlue);
    histPP_p->SetMarkerStyle(25);
  }

  TLine* zeroLine_p = new TLine(0., 0., 0.50, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(1);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  const char* cent2[4] = {"50-100%", "30-50%", "10-30%", "0-10%"};
  Float_t xCent[4] = {.35, .20, .20, .20};

  for(Int_t panelIter = 0; panelIter < 8; panelIter++){
    profPanel_p->cd(panelIter+1);

    if(panelIter < 4){
      if(strcmp(fileTagPP, "") != 0)
	histPP_p->DrawCopy("E1");

      histPbPb_p[panelIter] = (TH1F*)histPbPbFile_p->Get(Form("r%sMultA%sFCorr_%s_%s_h", algType[setNum], Tight, cent[panelIter], fileTagPbPb));

      std::cout << Form("r%sMultA%sFCorr_%s_%s_h", algType[setNum], Tight, cent[panelIter], fileTagPbPb) << std::endl;

      histPbPb_p[panelIter]->SetFillColor(kRed);
      histPbPb_p[panelIter]->SetLineColor(kRed);
      histPbPb_p[panelIter]->SetMarkerColor(kRed);
      histPbPb_p[panelIter]->SetMarkerStyle(28);
      histPbPb_p[panelIter]->DrawCopy("E1 SAME");

      label_p->DrawLatex(xCent[panelIter], .55, cent2[panelIter]);
    }
    else{
      histPbPb_p[panelIter - 4]->Add(histPP_p, -1);

      niceTH1(histPbPb_p[panelIter - 4], 29.999, -5, 505, 507);

      histPbPb_p[panelIter - 4]->SetMarkerStyle(20);
      histPbPb_p[panelIter - 4]->SetMarkerColor(1);
      histPbPb_p[panelIter - 4]->SetFillColor(1);
      histPbPb_p[panelIter - 4]->SetLineColor(1);

      histPbPb_p[panelIter - 4]->SetYTitle("PbPb-pp");
      histPbPb_p[panelIter - 4]->GetYaxis()->SetTitleOffset(2.2);
      histPbPb_p[panelIter - 4]->SetXTitle("A_{J}");
      histPbPb_p[panelIter - 4]->GetXaxis()->SetTitleOffset(1.6);
      histPbPb_p[panelIter - 4]->DrawCopy("E1");
      zeroLine_p->Draw("SAME");

      if(panelIter == 4){
	label_p->DrawLatex(.35, .85, "p_{T,1}>120 GeV/c");
	label_p->DrawLatex(.35, .75, "p_{T,2}>50 GeV/c");
	label_p->DrawLatex(.35, .65, "#Delta#phi_{1,2}>5#pi/6");
	label_p->DrawLatex(.35, .55, "|#eta_{1}|,|#eta_{2}|<1.6");
      }
      else if(panelIter == 5){
	label_p->DrawLatex(.20, .85, Form("anti-k_{T} %s R=%s", puVsString[setNum], radString[setNum]));
	label_p->DrawLatex(.20, .75, "|#eta_{trk}|<2.4");
	label_p->DrawLatex(.20, .65, "p_{T}^{trk}>0.5 GeV/c");
      }
    }

    profPanel_p->cd(panelIter+1)->RedrawAxis();
  }

  TLegend* legMult_p = new TLegend(0.45, 0.75, 0.75, 0.95);
  legMult_p->SetFillColor(0);
  legMult_p->SetFillStyle(0);
  legMult_p->SetTextFont(43);
  legMult_p->SetTextSizePixels(28);
  legMult_p->SetBorderSize(0);

  TH1F* dummHist_p = new TH1F("dummHist_p", "dummHist_p", 10, 0., 1.);
  dummHist_p->SetMarkerColor(kRed);
  dummHist_p->SetFillColor(kRed);
  dummHist_p->SetLineColor(kRed);
  dummHist_p->SetMarkerStyle(28);

  legMult_p->AddEntry(dummHist_p, "PbPb" /* 150 #mub^{-1}"*/, "P L");
  legMult_p->AddEntry(histPP_p, "pp 5.3 pb^{-1}", "P L");

  profPanel_p->cd(1);

  legMult_p->Draw("SAME");

  TFile* outFile_p = new TFile(outName, "UPDATE");
  profPanel_p->Write();
  claverCanvasSaving(profPanel_p, Form("pdfDir/%sMultA%sFCorrStack_%s", algType[setNum], Tight, fileTagPbPb), "pdf");
  outFile_p->Close();

  delete outFile_p;
  delete dummHist_p;
  delete legMult_p;
  delete label_p;
  delete zeroLine_p;
  delete profPanel_p;

  if(strcmp(fileTagPP, "") != 0){
    histPPFile_p->Close();
    delete histPPFile_p;
  }

  histPbPbFile_p->Close();
  delete histPbPbFile_p;
}


void makeImbPtStack(const char* filePbPbName, const char* fileTagPbPb, const char* outName, const char* gorr, Int_t setNum, const std::string projMult, const char* Corr = "", const char* CNCR = "", Bool_t montecarlo = false, const char* filePPName = "", const char* fileTagPP = "", const char* Tight = "", Bool_t doSpill = false, const char* filePbPbName2 = "", const char* fileTagPbPb2 = "", const char* filePPName2 = "", const char* fileTagPP2 = "", Bool_t isHighPtTrk = false)
{
  TFile* outFile_p = new TFile(outName, "UPDATE");

  TH1F* dumHist_p[7];
  for(Int_t dumIter = 0; dumIter < 7; dumIter++){
    dumHist_p[dumIter] = new TH1F(Form("dumHist_%d_h", dumIter), Form("dumHist_%d_h", dumIter), 10, 0, 1);
  }

  //  dumHist_p[0]->SetFillColor(kMagenta+4);
  dumHist_p[1]->SetFillColor(kBlue-9);
  dumHist_p[2]->SetFillColor(kYellow-9);
    dumHist_p[3]->SetFillColor(kOrange+1);
    dumHist_p[4]->SetFillColor(kGreen+3);
    dumHist_p[5]->SetFillColor(kRed+1);

  //  dumHist_p[3]->SetFillColor(kGray);
  //  dumHist_p[4]->SetFillColor(kGray+1);
  //  dumHist_p[5]->SetFillColor(kGray+2);

  dumHist_p[1]->SetMarkerStyle(25);
  dumHist_p[2]->SetMarkerStyle(28);
  dumHist_p[3]->SetMarkerStyle(24);

  TH1F* dumLineHist_p = new TH1F("dumLineHist_h", "dumLineHist_h", 10, 0, 1);
  dumLineHist_p->SetLineStyle(2);

  TLegend* legA_p = new TLegend(0.25, 0.22, 0.99, 0.88);
  TLegend* legB_p;
  TLegend* legC_p;

  TFile* histPbPbFile_p = new TFile(filePbPbName, "READ");
  TFile* histPPFile_p;
  if(strcmp(fileTagPP, "") != 0) histPPFile_p = new TFile(filePPName, "READ");
 
  //Grab hists for stack



  TH1F* hist1_p[7];
  TH1F* hist2_p[7];
  TH1F* hist3_p[7];
  TH1F* hist4_p[7];
  TH1F* histPP_p[7];

  grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_p);
  grabHistForStackPbPb(histPbPbFile_p, gorr, algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_p, hist2_p, hist3_p, hist4_p);

  TCanvas* profPanel_p;
  if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
    profPanel_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPP_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPP_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), 5*300, 2*350);
    makeMultiPanelCanvas(profPanel_p, 5, 2, 0.0, 0.0, 0.2, 0.2, 0.05);
    std::cout << "FivePanel Init" << std::endl;
  }
  else{
    profPanel_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPP_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPP_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), 300*3, 2*350);
    makeMultiPanelCanvas(profPanel_p, 3, 2, 0.0, 0.0, 0.2, 0.2, 0.01, CNCR);
    std::cout << "ThreePanel Init" << std::endl;
  }

  //Make legend

  if(isAll(CNCR)){
    legB_p = new TLegend(0.45, 0.1, 0.75, 0.40);
    legC_p = new TLegend(.20, .15, .45, .30);
    legC_p->SetFillColor(0);
    legC_p->SetFillStyle(0);
    legC_p->SetTextFont(43);
    legC_p->SetTextSizePixels(28);
    legC_p->SetBorderSize(0);
  }
  else legB_p = new TLegend(0.25, 0.1, 0.55, 0.40);

  legA_p->SetFillColor(0);
  legA_p->SetFillStyle(0);
  legA_p->SetTextFont(43);
  legA_p->SetTextSizePixels(28);
  legA_p->SetBorderSize(0);

  legB_p->SetFillColor(0);
  legB_p->SetFillStyle(0);
  legB_p->SetTextFont(43);
  legB_p->SetTextSizePixels(28);
  legB_p->SetBorderSize(0);

  histToCanvas_TOP(profPanel_p, histPP_p, hist1_p, hist2_p, hist3_p, hist4_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, false);

  profPanel_p->cd(1);
   legB_p->Draw();
  if(isAll(CNCR) && !strcmp(projMult.c_str(), "Proj")){
    profPanel_p->cd(2);
    legC_p->AddEntry(dumHist_p[1], "PbPb cumulative", "L");
       legC_p->Draw("SAME");
  }

  if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")) profPanel_p->cd(6);
  else profPanel_p->cd(4);
  legA_p->Draw("SAME");


  histPbPbFile_p->Close();
  delete histPbPbFile_p;
  histPbPbFile_p = new TFile(filePbPbName, "READ");

  if(strcmp(fileTagPP, "") != 0){
    histPPFile_p->Close();
    delete histPPFile_p;
    histPPFile_p = new TFile(filePPName, "READ");
  }

  grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_p);
  grabHistForStackPbPb(histPbPbFile_p, gorr, algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_p, hist2_p, hist3_p, hist4_p);



  for(Int_t histIter = 1; histIter < 7; histIter++){
    if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
      hist3_p[histIter]->Add(histPP_p[histIter], -1);
      hist4_p[histIter]->Add(histPP_p[histIter], -1);
    }
    hist1_p[histIter]->Add(histPP_p[histIter], -1);
    hist2_p[histIter]->Add(histPP_p[histIter], -1);
  }

  histToCanvas_BOTTOM(profPanel_p, histPP_p, hist1_p, hist2_p, hist3_p, hist4_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, false);

  legB_p->AddEntry(dumHist_p[1], "pp", "p");
  legB_p->AddEntry(dumHist_p[2], "PbPb", "p");
  if(isAll(CNCR)) legC_p->AddEntry(dumLineHist_p, "pp cumulative", "L");
  legB_p->AddEntry(dumHist_p[3], "PbPb - pp", "p");

  //  legA_p->AddEntry(dumHist_p[0], "0.0-0.5", "F");
  legA_p->AddEntry(dumHist_p[1], "0.5-1.0", "F");
  legA_p->AddEntry(dumHist_p[2], "1.0-2.0", "F");
  legA_p->AddEntry(dumHist_p[3], "2.0-4.0", "F");
  legA_p->AddEntry(dumHist_p[4], "4.0-8.0", "F");
  legA_p->AddEntry(dumHist_p[5], "8.0-300.0", "F");

  TFile* histPbPbFile2_p;
  TFile* histPPFile2_p;
  TH1F* hist1_2_p[7];
  TH1F* hist2_2_p[7];
  TH1F* hist3_2_p[7];
  TH1F* hist4_2_p[7];
  TH1F* histPP_2_p[7];
  TCanvas* profPanel2_p;
  TCanvas* profPanel3_p;
  TCanvas* profPanel4_p;
  TH1F* hist1_3_p[7];
  TH1F* hist2_3_p[7];
  TH1F* hist3_3_p[7];
  TH1F* hist4_3_p[7];
  TH1F* histPP_3_p[7];

  if(strcmp(filePbPbName2, "") != 0  && false){
    histPbPbFile_p->Close();
    delete histPbPbFile_p;
    histPbPbFile_p = new TFile(filePbPbName, "READ");

    if(strcmp(fileTagPP, "") != 0){
      histPPFile_p->Close();
      delete histPPFile_p;
      histPPFile_p = new TFile(filePPName, "READ");
    }

    grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_p);
    grabHistForStackPbPb(histPbPbFile_p, gorr, algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_p, hist2_p, hist3_p, hist4_p);

    histPbPbFile2_p = new TFile(filePbPbName2, "READ");    
    histPPFile2_p = new TFile(filePPName2, "READ");    

    grabHistForStackPP(histPPFile2_p, gorr, algTypePP[setNum], projMult, CNCR, Corr, fileTagPP2, histPP_2_p);
    grabHistForStackPbPb(histPbPbFile2_p, gorr, algType[setNum], projMult, CNCR, Corr, fileTagPbPb2, hist1_2_p, hist2_2_p, hist3_2_p, hist4_2_p);

    for(Int_t histIter = 0; histIter < 7; histIter++){
      histPP_p[histIter]->Add(histPP_2_p[histIter], -1);
      hist1_p[histIter]->Add(hist1_2_p[histIter], -1);
      hist2_p[histIter]->Add(hist2_2_p[histIter], -1);
      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
	hist3_p[histIter]->Add(hist3_2_p[histIter], -1);
	hist4_p[histIter]->Add(hist4_2_p[histIter], -1);
      }
    }

    if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
      profPanel2_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPPSubOldNew_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPPSubOldNew_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), 5*300, 2*350);
      makeMultiPanelCanvas(profPanel2_p, 5, 2, 0.0, 0.0, 0.2, 0.2, 0.05);
      std::cout << "FivePanel Init" << std::endl;
    }
    else{
      profPanel2_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPPSubOldNew_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPPSubOldNew_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), 300*3, 2*350);
      makeMultiPanelCanvas(profPanel2_p, 3, 2, 0.0, 0.0, 0.2, 0.2, 0.01, CNCR);
      std::cout << "ThreePanel Init" << std::endl;
    }
    histToCanvas_TOP(profPanel2_p, histPP_p, hist1_p, hist2_p, hist3_p, hist4_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, true);

    profPanel2_p->cd(1);
      legB_p->Draw();
    if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")) profPanel2_p->cd(6);
    else profPanel2_p->cd(4);
    legA_p->Draw("SAME");

    histPbPbFile_p->Close();
    delete histPbPbFile_p;
    histPbPbFile_p = new TFile(filePbPbName, "READ");

    histPbPbFile2_p->Close();
    delete histPbPbFile2_p;
    histPbPbFile2_p = new TFile(filePbPbName2, "READ");

    if(strcmp(fileTagPP, "") != 0){
      histPPFile_p->Close();
      delete histPPFile_p;
      histPPFile_p = new TFile(filePPName, "READ");

      histPPFile2_p->Close();
      delete histPPFile2_p;
      histPPFile2_p = new TFile(filePPName2, "READ");
    }

    grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_p);
    grabHistForStackPbPb(histPbPbFile_p, gorr, algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_p, hist2_p, hist3_p, hist4_p);
    for(Int_t histIter = 0; histIter < 7; histIter++){
      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
	hist3_p[histIter]->Add(histPP_p[histIter], -1);
	hist4_p[histIter]->Add(histPP_p[histIter], -1);
      }
      hist1_p[histIter]->Add(histPP_p[histIter], -1);
      hist2_p[histIter]->Add(histPP_p[histIter], -1);
    }

    grabHistForStackPP(histPPFile2_p, gorr, algTypePP[setNum], projMult, CNCR, Corr, fileTagPP2, histPP_2_p);
    grabHistForStackPbPb(histPbPbFile2_p, gorr, algType[setNum], projMult, CNCR, Corr, fileTagPbPb2, hist1_2_p, hist2_2_p, hist3_2_p, hist4_2_p);
    for(Int_t histIter = 0; histIter < 7; histIter++){
      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
	hist3_2_p[histIter]->Add(histPP_2_p[histIter], -1);
	hist4_2_p[histIter]->Add(histPP_2_p[histIter], -1);
      }
      hist1_2_p[histIter]->Add(histPP_2_p[histIter], -1);
      hist2_2_p[histIter]->Add(histPP_2_p[histIter], -1);
    }
    for(Int_t histIter = 0; histIter < 7; histIter++){
      hist1_p[histIter]->Add(hist1_2_p[histIter], -1);
      hist2_p[histIter]->Add(hist2_2_p[histIter], -1);
      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
	hist3_p[histIter]->Add(hist3_2_p[histIter], -1);
	hist4_p[histIter]->Add(hist4_2_p[histIter], -1);
      }
    }

    histToCanvas_BOTTOM(profPanel2_p, histPP_p, hist1_p, hist2_p, hist3_p, hist4_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, true);
  }

  if(montecarlo && (setNum == 5 || setNum == 11 || setNum == 17 || setNum == 19 || setNum == 23 || setNum == 29 || setNum == 31 || setNum == 33) && !strcmp(gorr, "r")){
    std::cout << setNum << std::endl;

    histPbPbFile_p->Close();
    delete histPbPbFile_p;
    histPbPbFile_p = new TFile(filePbPbName, "READ");

    if(strcmp(fileTagPP, "") != 0){
      histPPFile_p->Close();
      delete histPPFile_p;
      histPPFile_p = new TFile(filePPName, "READ");
    }

    grabHistForStackPP(histPPFile_p, "r", algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_p);
    grabHistForStackPbPb(histPbPbFile_p, "r", algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_p, hist2_p, hist3_p, hist4_p);

    grabHistForStackPP(histPPFile_p, "g", algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_2_p);
    grabHistForStackPbPb(histPbPbFile_p, "g", algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_2_p, hist2_2_p, hist3_2_p, hist4_2_p);

    for(Int_t histIter = 0; histIter < 7; histIter++){
      histPP_p[histIter]->Add(histPP_2_p[histIter], -1);
      hist1_p[histIter]->Add(hist1_2_p[histIter], -1);
      hist2_p[histIter]->Add(hist2_2_p[histIter], -1);

      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
        hist3_p[histIter]->Add(hist3_2_p[histIter], -1);
        hist4_p[histIter]->Add(hist4_2_p[histIter], -1);
      }
    }

    /*
    std::cout << std::endl;
    std::cout << "Tracking" << std::endl;
    std::cout << "PP" << std::endl;
    for(Int_t iter = 0; iter < 10; iter++){
      std::cout << histPP_p[6]->GetBinContent(iter+1) << std::endl;
    }

    std::cout << std::endl;
    std::cout << "30100" << std::endl;
    for(Int_t iter = 0; iter < 10; iter++){
      std::cout << hist1_p[6]->GetBinContent(iter+1) << std::endl;
    }

    std::cout << std::endl;
    std::cout << "030" << std::endl;
    for(Int_t iter = 0; iter < 10; iter++){
      std::cout << hist2_p[6]->GetBinContent(iter+1) << std::endl;
    }
    */

    if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
      profPanel3_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPPSubRVs_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPPSubRVs_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), 5*300, 2*350);
      makeMultiPanelCanvas(profPanel3_p, 5, 2, 0.0, 0.0, 0.2, 0.2, 0.05);
      std::cout << "FivePanel Init" << std::endl;
    }
    else{
      profPanel3_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPPSubRVs_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPPSubRVs_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), 300*3, 2*350);
      makeMultiPanelCanvas(profPanel3_p, 3, 2, 0.0, 0.0, 0.2, 0.2, 0.01, CNCR);
      std::cout << "ThreePanel Init" << std::endl;
    }

    histToCanvas_TOP(profPanel3_p, histPP_p, hist1_p, hist2_p, hist3_p, hist4_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, true);

    profPanel3_p->cd(1);
    legB_p->Draw();
    if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")) profPanel3_p->cd(6);
    else profPanel3_p->cd(4);
    legA_p->Draw("SAME");

    histPbPbFile_p->Close();
    delete histPbPbFile_p;
    histPbPbFile_p = new TFile(filePbPbName, "READ");

    if(strcmp(fileTagPP, "") != 0){
      histPPFile_p->Close();
      delete histPPFile_p;
      histPPFile_p = new TFile(filePPName, "READ");
    }

    grabHistForStackPP(histPPFile_p, "r", algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_p);
    grabHistForStackPbPb(histPbPbFile_p, "r", algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_p, hist2_p, hist3_p, hist4_p);
    for(Int_t histIter = 0; histIter < 7; histIter++){
      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
        hist3_p[histIter]->Add(histPP_p[histIter], -1);
        hist4_p[histIter]->Add(histPP_p[histIter], -1);
      }
      hist1_p[histIter]->Add(histPP_p[histIter], -1);
      hist2_p[histIter]->Add(histPP_p[histIter], -1);
    }


    grabHistForStackPP(histPPFile_p, "g", algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_2_p);
    grabHistForStackPbPb(histPbPbFile_p, "g", algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_2_p, hist2_2_p, hist3_2_p, hist4_2_p);
    for(Int_t histIter = 0; histIter < 7; histIter++){
      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
        hist3_2_p[histIter]->Add(histPP_2_p[histIter], -1);
        hist4_2_p[histIter]->Add(histPP_2_p[histIter], -1);
      }
      hist1_2_p[histIter]->Add(histPP_2_p[histIter], -1);
      hist2_2_p[histIter]->Add(histPP_2_p[histIter], -1);
    }


    for(Int_t histIter = 0; histIter < 7; histIter++){
      hist1_p[histIter]->Add(hist1_2_p[histIter], -1);
      hist2_p[histIter]->Add(hist2_2_p[histIter], -1);
      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
        hist3_p[histIter]->Add(hist3_2_p[histIter], -1);
        hist4_p[histIter]->Add(hist4_2_p[histIter], -1);
      }
    }
    /*
    std::cout << std::endl;
    std::cout << "30100 Sub" << std::endl;
    for(Int_t iter = 0; iter < 10; iter++){
      std::cout << hist1_p[6]->GetBinContent(iter+1) << std::endl;
    }

    std::cout << std::endl;
    std::cout << "030 Sub" << std::endl;
    for(Int_t iter = 0; iter < 10; iter++){
      std::cout << hist2_p[6]->GetBinContent(iter+1) << std::endl;
    }
    */

    histToCanvas_BOTTOM(profPanel3_p, histPP_p, hist1_p, hist2_p, hist3_p, hist4_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, true);
  }

  if(montecarlo && (setNum == 5 || setNum == 11 || setNum == 17 || setNum == 19 || setNum == 23 || setNum == 29  || setNum == 31 || setNum == 33) && !strcmp(gorr, "g")){
    std::cout << setNum << std::endl;
    histPbPbFile_p->Close();
    delete histPbPbFile_p;
    histPbPbFile_p = new TFile(filePbPbName, "READ");

    if(strcmp(fileTagPP, "") != 0){
      histPPFile_p->Close();
      delete histPPFile_p;
      histPPFile_p = new TFile(filePPName, "READ");
    }

    grabHistForStackPP(histPPFile_p, "g", algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_p);
    grabHistForStackPbPb(histPbPbFile_p, "g", algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_p, hist2_p, hist3_p, hist4_p);

    grabHistForStackPP(histPPFile_p, "g", truthType[setNum], projMult, CNCR, Corr, fileTagPP, histPP_2_p);
    grabHistForStackPbPb(histPbPbFile_p, "g", truthType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_2_p, hist2_2_p, hist3_2_p, hist4_2_p);

    for(Int_t histIter = 0; histIter < 7; histIter++){
      histPP_p[histIter]->Add(histPP_2_p[histIter], -1);
      hist1_p[histIter]->Add(hist1_2_p[histIter], -1);
      hist2_p[histIter]->Add(hist2_2_p[histIter], -1);
      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
        hist3_p[histIter]->Add(hist3_2_p[histIter], -1);
        hist4_p[histIter]->Add(hist4_2_p[histIter], -1);
      }
    }
    
     
    //previous spillover preserved temporarily
    /*
    if(!montecarlo && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 22|| setNum == 23 || setNum == 24 || setNum == 25)){
      histPbPbFile2_p = new TFile(filePbPbName2, "READ");
      histPPFile2_p = new TFile(filePPName2, "READ");

      grabHistForStackPP(histPPFile2_p, "r", algTypePP[setNum], projMult, CNCR, Corr, Tight, fileTagPP2, histPP_3_p);
      grabHistForStackPbPb(histPbPbFile2_p, "r", algType[setNum], projMult, CNCR, Corr, Tight, fileTagPbPb2, hist1_3_p, hist2_3_p, hist3_3_p, hist4_3_p);

      for(Int_t histIter = 0; histIter < 7; histIter++){
	histPP_3_p[histIter]->Add(histPP_p[histIter], -1);
	hist1_3_p[histIter]->Add(hist1_p[histIter], -1);
	hist2_3_p[histIter]->Add(hist2_p[histIter], -1);
	if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
	  hist3_3_p[histIter]->Add(hist3_p[histIter], -1);
	  hist4_3_p[histIter]->Add(hist4_p[histIter], -1);
	}
      }

      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
	profPanel2_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPP_SPILL_%s_c", "r", algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPP_SPILL_%s_c", "r", algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb2), 5*300, 2*350);
	makeMultiPanelCanvas(profPanel2_p, 5, 2, 0.0, 0.0, 0.2, 0.2, 0.05);
	std::cout << "FivePanel Init" << std::endl;
      }
      else{
	profPanel2_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPP_SPILL_%s_c", "r", algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPP_SPILL_%s_c", "r", algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb2), 300*3, 2*350);
	makeMultiPanelCanvas(profPanel2_p, 3, 2, 0.0, 0.0, 0.2, 0.2, 0.01, CNCR);
	std::cout << "ThreePanel Init" << std::endl;
      }

      histToCanvas_TOP(profPanel2_p, histPP_3_p, hist1_3_p, hist2_3_p, hist3_3_p, hist4_3_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, false, true);

      profPanel2_p->cd(1);
      legB_p->Draw();
      if(isAll(CNCR) && !strcmp(projMult.c_str(), "Proj")){
	profPanel2_p->cd(2);
	legC_p->Draw("SAME");
      }

      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")) profPanel2_p->cd(6);
      else profPanel2_p->cd(4);
      legA_p->Draw("SAME");

      histPbPbFile2_p->Close();
      delete histPbPbFile2_p;
      histPbPbFile2_p = new TFile(filePbPbName2, "READ");

      histPPFile2_p->Close();
      delete histPPFile2_p;
      histPPFile2_p = new TFile(filePPName2, "READ");

      grabHistForStackPP(histPPFile2_p, "r", algTypePP[setNum], projMult, CNCR, Corr, Tight, fileTagPP2, histPP_3_p);
      grabHistForStackPbPb(histPbPbFile2_p, "r", algType[setNum], projMult, CNCR, Corr, Tight, fileTagPbPb2, hist1_3_p, hist2_3_p, hist3_3_p, hist4_3_p);

      for(Int_t histIter = 0; histIter < 7; histIter++){
        histPP_3_p[histIter]->Add(histPP_p[histIter], -1);
        hist1_3_p[histIter]->Add(hist1_p[histIter], -1);
        hist2_3_p[histIter]->Add(hist2_p[histIter], -1);
        if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
          hist3_3_p[histIter]->Add(hist3_p[histIter], -1);
          hist4_3_p[histIter]->Add(hist4_p[histIter], -1);
        }

        hist1_3_p[histIter]->Add(histPP_3_p[histIter], -1);
        hist2_3_p[histIter]->Add(histPP_3_p[histIter], -1);
        if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
          hist3_3_p[histIter]->Add(histPP_3_p[histIter], -1);
          hist4_3_p[histIter]->Add(histPP_3_p[histIter], -1);
        }
      }

      histToCanvas_BOTTOM(profPanel2_p, histPP_3_p, hist1_3_p, hist2_3_p, hist3_3_p, hist4_3_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, false, true);
    }
    */
    /*
    std::cout << std::endl;
    std::cout << "Jet" << std::endl;
    std::cout << "PP" << std::endl;
    for(Int_t iter = 0; iter < 10; iter++){
      std::cout << histPP_p[6]->GetBinContent(iter+1) << std::endl;
    }

    std::cout << std::endl;
    std::cout << "30100" << std::endl;
    for(Int_t iter = 0; iter < 10; iter++){
      std::cout << hist1_p[6]->GetBinContent(iter+1) << std::endl;
    }

    std::cout << std::endl;
    std::cout << "030" << std::endl;
    for(Int_t iter = 0; iter < 10; iter++){
      std::cout << hist2_p[6]->GetBinContent(iter+1) << std::endl;
    }
    */
  
    if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
      profPanel4_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPPSubgTVs_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPPSubgTVs_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), 5*300, 2*350);
      makeMultiPanelCanvas(profPanel4_p, 5, 2, 0.0, 0.0, 0.2, 0.2, 0.05);
      std::cout << "FivePanel Init" << std::endl;
    }
    else{
      profPanel4_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPPSubgTVs_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPPSubgTVs_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), 300*3, 2*350);
      makeMultiPanelCanvas(profPanel4_p, 3, 2, 0.0, 0.0, 0.2, 0.2, 0.01, CNCR);
      std::cout << "ThreePanel Init" << std::endl;
    }

    histToCanvas_TOP(profPanel4_p, histPP_p, hist1_p, hist2_p, hist3_p, hist4_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, true);
 
    profPanel4_p->cd(1);
    legB_p->Draw();
    if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")) profPanel4_p->cd(6);
    else profPanel4_p->cd(4);
    legA_p->Draw("SAME");

    histPbPbFile_p->Close();
    delete histPbPbFile_p;
    histPbPbFile_p = new TFile(filePbPbName, "READ");

    if(strcmp(fileTagPP, "") != 0){
      histPPFile_p->Close();
      delete histPPFile_p;
      histPPFile_p = new TFile(filePPName, "READ");
    }

    grabHistForStackPP(histPPFile_p, "g", algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_p);
    grabHistForStackPbPb(histPbPbFile_p, "g", algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_p, hist2_p, hist3_p, hist4_p);
    for(Int_t histIter = 0; histIter < 7; histIter++){
      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
        hist3_p[histIter]->Add(histPP_p[histIter], -1);
        hist4_p[histIter]->Add(histPP_p[histIter], -1);
      }
      hist1_p[histIter]->Add(histPP_p[histIter], -1);
      hist2_p[histIter]->Add(histPP_p[histIter], -1);
    }

    grabHistForStackPP(histPPFile_p, "g", truthType[setNum], projMult, CNCR, Corr, fileTagPP, histPP_2_p);
    grabHistForStackPbPb(histPbPbFile_p, "g", truthType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_2_p, hist2_2_p, hist3_2_p, hist4_2_p);
    for(Int_t histIter = 0; histIter < 7; histIter++){
      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
        hist3_2_p[histIter]->Add(histPP_2_p[histIter], -1);
        hist4_2_p[histIter]->Add(histPP_2_p[histIter], -1);
      }
      hist1_2_p[histIter]->Add(histPP_2_p[histIter], -1);
      hist2_2_p[histIter]->Add(histPP_2_p[histIter], -1);
    }


    for(Int_t histIter = 0; histIter < 7; histIter++){
      hist1_p[histIter]->Add(hist1_2_p[histIter], -1);
      hist2_p[histIter]->Add(hist2_2_p[histIter], -1);
      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
        hist3_p[histIter]->Add(hist3_2_p[histIter], -1);
        hist4_p[histIter]->Add(hist4_2_p[histIter], -1);
      }
    }

    /*
    std::cout << std::endl;
    std::cout << "30100 Sub" << std::endl;
    for(Int_t iter = 0; iter < 10; iter++){
      std::cout << hist1_p[6]->GetBinContent(iter+1) << std::endl;
    }

    std::cout << std::endl;
    std::cout << "030 Sub" << std::endl;
    for(Int_t iter = 0; iter < 10; iter++){
      std::cout << hist2_p[6]->GetBinContent(iter+1) << std::endl;
    }
    */

    histToCanvas_BOTTOM(profPanel4_p, histPP_p, hist1_p, hist2_p, hist3_p, hist4_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, true);
  }
  
  if(!montecarlo && doSpill && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 22|| setNum == 23 || setNum == 24 || setNum == 25 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31 || setNum == 32 || setNum == 33 || setNum == 34 || setNum == 35)){
    std::cout << setNum << std::endl;
    std::cout << "CHECK" << std::endl;
    histPbPbFile_p->Close();
    delete histPbPbFile_p;
    histPbPbFile_p = new TFile(filePbPbName, "READ");

    if(strcmp(fileTagPP, "") != 0){
      histPPFile_p->Close();
      delete histPPFile_p;
      histPPFile_p = new TFile(filePPName, "READ");
    }

    std::cout << "CHECK" << std::endl;

    grabHistForStackPP(histPPFile_p, "r", algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_p);
    grabHistForStackPbPb(histPbPbFile_p, "r", algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_p, hist2_p, hist3_p, hist4_p);

    if(setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 28 || setNum == 29 || setNum == 30 || setNum == 31){
      histPPFile2_p = new TFile(Form("spilloverFits/spillover_CaloResProj_A%s_pp.root", CNCR), "READ");
      histPbPbFile2_p = new TFile(Form("spilloverFits/spillover_CaloResProj_A%s.root", CNCR), "READ");


      std::cout << Form("spilloverFits/spillover_CaloResProj_A%s_pp.root", CNCR) << std::endl;
      std::cout << Form("spilloverFits/spillover_CaloResProj_A%s.root", CNCR) << std::endl;
      
      grabSpillForStackPP(histPPFile2_p, "Res", radString2[setNum], histPP_3_p);
      grabSpillForStackPbPb(histPbPbFile2_p, "Res", radString2[setNum], CNCR, hist1_3_p, hist2_3_p, hist3_3_p, hist4_3_p);
    }
    else if(setNum == 22 || setNum == 23 || setNum == 24 || setNum == 25 || setNum == 32 || setNum == 33 || setNum == 34 || setNum == 35){
      histPPFile2_p = new TFile(Form("spilloverFits/spillover_CaloSwapProj_A%s_pp.root", CNCR), "READ");
      histPbPbFile2_p = new TFile(Form("spilloverFits/spillover_CaloSwapProj_A%s.root", CNCR), "READ");
      
      grabSpillForStackPP(histPPFile2_p, "Swap", radString2[setNum], histPP_3_p);
      grabSpillForStackPbPb(histPbPbFile2_p, "Swap", radString2[setNum], CNCR, hist1_3_p, hist2_3_p, hist3_3_p, hist4_3_p);
    }

    for(Int_t iter = 1; iter < 7; iter++){
      histPP_p[iter]->Add(histPP_3_p[iter], -1);

      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
        hist3_p[iter]->Add(hist3_3_p[iter], -1);
        hist4_p[iter]->Add(hist4_3_p[iter], -1);
      }
      hist1_p[iter]->Add(hist1_3_p[iter], -1);
      hist2_p[iter]->Add(hist2_3_p[iter], -1);
    }
    
    if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
      profPanel2_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPP_SPILL_%s_c", "r", algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPP_SPILL_%s_c", "r", algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb2), 5*300, 2*350);
      makeMultiPanelCanvas(profPanel2_p, 5, 2, 0.0, 0.0, 0.2, 0.2, 0.05);
      std::cout << "FivePanel Init" << std::endl;
    }
    else{
      profPanel2_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPP_SPILL_%s_c", "r", algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPP_SPILL_%s_c", "r", algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb2), 300*3, 2*350);
      makeMultiPanelCanvas(profPanel2_p, 3, 2, 0.0, 0.0, 0.2, 0.2, 0.01, CNCR);
      std::cout << "ThreePanel Init" << std::endl;
    }
    
    histToCanvas_TOP(profPanel2_p, histPP_p, hist1_p, hist2_p, hist3_p, hist4_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, false, true);     

    profPanel2_p->cd(1);
    legB_p->Draw();
    if(isAll(CNCR) && !strcmp(projMult.c_str(), "Proj")){
      profPanel2_p->cd(2);
      legC_p->Draw("SAME");
    }
    
    if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")) profPanel2_p->cd(6);
    else profPanel2_p->cd(4);
    legA_p->Draw("SAME");  

    histPbPbFile_p->Close();
    delete histPbPbFile_p;
    histPbPbFile_p = new TFile(filePbPbName, "READ");

    if(strcmp(fileTagPP, "") != 0){
      histPPFile_p->Close();
      delete histPPFile_p;
      histPPFile_p = new TFile(filePPName, "READ");
    }

    grabHistForStackPP(histPPFile_p, "r", algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP_p);
    grabHistForStackPbPb(histPbPbFile_p, "r", algType[setNum], projMult, CNCR, Corr, fileTagPbPb, hist1_p, hist2_p, hist3_p, hist4_p);

    for(Int_t iter = 1; iter < 7; iter++){
      histPP_p[iter]->Add(histPP_3_p[iter], -1);

      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
        hist3_p[iter]->Add(hist3_3_p[iter], -1);
        hist4_p[iter]->Add(hist4_3_p[iter], -1);
      }
      hist1_p[iter]->Add(hist1_3_p[iter], -1);
      hist2_p[iter]->Add(hist2_3_p[iter], -1);

      if(!strcmp(CNCR, "") || !strcmp(CNCR, "C") || !strcmp(CNCR, "NC")){
        hist3_p[iter]->Add(histPP_p[iter], -1);
        hist4_p[iter]->Add(histPP_p[iter], -1);
      }
      hist1_p[iter]->Add(histPP_p[iter], -1);
      hist2_p[iter]->Add(histPP_p[iter], -1);
    }

    histToCanvas_BOTTOM(profPanel2_p, histPP_p, hist1_p, hist2_p, hist3_p, hist4_p, setNum, CNCR, projMult, Tight, montecarlo, isHighPtTrk, false, true);
  }
      

  outFile_p->cd();
  profPanel_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(profPanel_p, Form("pdfDir/%s%s%sA%s%s%sPTStack_%s", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), "pdf");
  if(doSpill && !montecarlo && (setNum == 16 || setNum == 17 || setNum == 18 || setNum == 19 || setNum == 22|| setNum == 23 || setNum == 24 || setNum == 25)) claverCanvasSaving(profPanel2_p, Form("pdfDir/%s%s%sA%s%s%sPTStack_SPILL_%s", "r", algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), "pdf");
  if(montecarlo && (setNum == 5 || setNum == 11 || setNum == 17 || setNum == 19 || setNum == 23 || setNum == 29 || setNum == 31 || setNum == 33) && !strcmp(gorr, "r")) claverCanvasSaving(profPanel3_p, Form("pdfDir/%s%s%sA%s%s%sPTStackSubRVs_%s", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), "pdf");
  if(montecarlo && (setNum == 5 || setNum == 11 || setNum == 17 || setNum == 19 || setNum == 23 || setNum == 29 || setNum == 31 || setNum == 33) && !strcmp(gorr, "g")) claverCanvasSaving(profPanel4_p, Form("pdfDir/%s%s%sA%s%s%sPTStackSubgTVs_%s", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, Tight, fileTagPbPb), "pdf");

  if(montecarlo && (setNum == 5 || setNum == 11 || setNum == 17  || setNum == 19 || setNum == 23 || setNum == 29 || setNum == 31 || setNum == 33) && !strcmp(gorr, "g")) delete profPanel4_p;

  if(montecarlo && (setNum == 5 || setNum == 11 || setNum == 17 || setNum == 19 || setNum == 23 || setNum == 29 || setNum == 31 || setNum == 33) && !strcmp(gorr, "r")) delete profPanel3_p;

  if(strcmp(filePbPbName2, "") != 0 && (setNum == 15 || setNum == 16 || setNum == 17 || setNum == 18) && !strcmp(gorr, "g")) delete profPanel2_p;

  if(isAll(CNCR))  delete legC_p;
  delete legB_p;

  delete profPanel_p;

  if(strcmp(fileTagPP, "") != 0){
    histPPFile_p->Close();
    delete histPPFile_p;
  }

  histPbPbFile_p->Close();
  delete histPbPbFile_p;

  delete legA_p;
  delete dumLineHist_p;
  dumLineHist_p = 0;

  for(Int_t dumIter = 0; dumIter < 7; dumIter++){
    delete dumHist_p[dumIter];
    dumHist_p[dumIter] = 0;
  }

  outFile_p->Close();
  delete outFile_p;

  return;
}


void makeImbPtRStack(const char* filePbPbName, const char* fileTagPbPb, const char* outName, const char* gorr, Int_t setNum, const std::string projMult, const char* Corr = "", const char* CNCR = "", const std::string cent = "030", Bool_t montecarlo = false, const char* filePPName = "", const char* fileTagPP = "")
{
  TFile* outFile_p = new TFile(outName, "UPDATE");

  TH1F* dumHist_p[7];
  for(Int_t dumIter = 0; dumIter < 7; dumIter++){
    dumHist_p[dumIter] = new TH1F(Form("dumHist_%d_h", dumIter), Form("dumHist_%d_h", dumIter), 10, 0, 1);
  }

  //  dumHist_p[0]->SetFillColor(kMagenta+4);
  dumHist_p[1]->SetFillColor(kBlue-9);
  dumHist_p[2]->SetFillColor(kYellow-9);
    dumHist_p[3]->SetFillColor(kOrange+1);
    dumHist_p[4]->SetFillColor(kGreen+3);
    dumHist_p[5]->SetFillColor(kRed+1);

  //  dumHist_p[3]->SetFillColor(kGray);
  //  dumHist_p[4]->SetFillColor(kGray+1);
  //  dumHist_p[5]->SetFillColor(kGray+2);

  dumHist_p[1]->SetMarkerStyle(25);
  dumHist_p[2]->SetMarkerStyle(28);
  dumHist_p[3]->SetMarkerStyle(24);

  TH1F* dumLineHist_p = new TH1F("dumLineHist_h", "dumLineHist_h", 10, 0, 1);
  dumLineHist_p->SetLineStyle(2);

  TLegend* legA_p = new TLegend(0.25, 0.22, 0.99, 0.44);
  TLegend* legD_p = new TLegend(0.25, 0.22, 0.99, 0.44);
  TLegend* legB_p;
  TLegend* legC_p;

  if(isAll(CNCR)){
    legB_p = new TLegend(0.45, 0.1, 0.75, 0.40);
    legC_p = new TLegend(.20, .15, .45, .30);
    legC_p->SetFillColor(0);
    legC_p->SetFillStyle(0);
    legC_p->SetTextFont(43);
    legC_p->SetTextSizePixels(28);
    legC_p->SetBorderSize(0);
  }
  else legB_p = new TLegend(0.25, 0.1, 0.55, 0.40);

  legA_p->SetFillColor(0);
  legA_p->SetFillStyle(0);
  legA_p->SetTextFont(43);
  legA_p->SetTextSizePixels(28);
  legA_p->SetBorderSize(0);

  legD_p->SetFillColor(0);
  legD_p->SetFillStyle(0);
  legD_p->SetTextFont(43);
  legD_p->SetTextSizePixels(28);
  legD_p->SetBorderSize(0);

  legB_p->SetFillColor(0);
  legB_p->SetFillStyle(0);
  legB_p->SetTextFont(43);
  legB_p->SetTextSizePixels(28);
  legB_p->SetBorderSize(0);


  TFile* histPbPbFile_p = new TFile(filePbPbName, "READ");
  TFile* histPPFile_p = new TFile(filePPName, "READ");

  TH1F* histPbPb1_p[7];
  TH1F* histPbPb2_p[7];
  TH1F* histPbPb3_p[7];
  TH1F* histPbPb4_p[7];

  TH1F* histPP1_p[7];
  TH1F* histPP2_p[7];
  TH1F* histPP3_p[7];
  TH1F* histPP4_p[7];

  const Float_t xCoord1[4] = {.38, .20, .20, .20};
  const std::string mcStr1[3] = {"PYT.", "PYT.+HYD.", "(P+H)-P"};
  const std::string dataStr1[3] = {"pp", "PbPb", "PbPb-pp"};

  std::string dataMCStr1[3];
  if(montecarlo){
    for(Int_t iter = 0; iter < 3; iter++){
      dataMCStr1[iter] = mcStr1[iter];
    }
  }
  else{
    for(Int_t iter = 0; iter < 3; iter++){
      dataMCStr1[iter] = dataStr1[iter];
    }
  }

  grabHistForStackPbPbCent(histPbPbFile_p, gorr, algType[setNum], projMult, CNCR, Corr, cent, fileTagPbPb, histPbPb1_p);
  grabHistForStackPbPbCent(histPbPbFile_p, gorr, algType[setNum+1], projMult, CNCR, Corr, cent, fileTagPbPb, histPbPb2_p);
  grabHistForStackPbPbCent(histPbPbFile_p, gorr, algType[setNum+2], projMult, CNCR, Corr, cent, fileTagPbPb, histPbPb3_p);
  grabHistForStackPbPbCent(histPbPbFile_p, gorr, algType[setNum+3], projMult, CNCR, Corr, cent, fileTagPbPb, histPbPb4_p);

  grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP1_p);
  grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum+1], projMult, CNCR, Corr, fileTagPP, histPP2_p);
  grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum+2], projMult, CNCR, Corr, fileTagPP, histPP3_p);
  grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum+3], projMult, CNCR, Corr, fileTagPP, histPP4_p);


  TCanvas* profPanel_p = new TCanvas(Form("%s%s%sA%s%s%sPTStackPP_RDEP%s_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, "", cent.c_str(), fileTagPbPb), Form("%s%s%sA%s%s%sPTStackPP_RDEP%s_%s_c", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, "", cent.c_str(), fileTagPbPb), 300*4, 3*350);
  profPanel_p->Divide(4, 3, 0.0, 0.0);
  std::cout << "Four Panel Init" << std::endl;

  std::cout << cent << std::endl;
  std::cout << CNCR << std::endl;

  histToCanvasR(profPanel_p, 1, histPP1_p, histPP2_p, histPP3_p, histPP4_p, setNum, CNCR, projMult, montecarlo, "", false);

  TLatex* label1_p = new TLatex();
  label1_p->SetNDC();
  label1_p->SetTextFont(43);
  label1_p->SetTextSizePixels(28);

  for(Int_t iter = 0; iter < 4; iter++){
    profPanel_p->cd(iter+1);
    label1_p->DrawLatex(xCoord1[iter], .92, Form("%s R=0.%d", dataMCStr1[0].c_str(), iter+2));
    if(iter == 3) label1_p->DrawLatex(xCoord1[iter], .84, "CMS Preliminary");
  }

  profPanel_p->cd(1);
  legB_p->Draw();
  if(isAll(CNCR) && !strcmp(projMult.c_str(), "Proj")){
    profPanel_p->cd(2);
    legC_p->AddEntry(dumHist_p[1], Form("%s cumulative", dataMCStr1[1].c_str()), "L");
    legC_p->Draw("SAME");
  }


  std::cout << "ELLO1" << std::endl;
  histToCanvasR(profPanel_p, 5, histPbPb1_p, histPbPb2_p, histPbPb3_p, histPbPb4_p, setNum, CNCR, projMult, montecarlo, cent, false);

  std::cout << "ELLO2" << std::endl;

  for(Int_t iter = 0; iter < 4; iter++){
    profPanel_p->cd(iter+5);
    label1_p->DrawLatex(xCoord1[iter], .92, Form("%s R=0.%d", dataMCStr1[1].c_str(), iter+2));
    label1_p->DrawLatex(xCoord1[iter]+.25, .84, Form("%s%%", cent.c_str()));
  }

  if(isAll(CNCR) || isInc(CNCR) || isU(CNCR) || isD(CNCR)){
    profPanel_p->cd(5);
    histPP1_p[6]->DrawCopy("HIST SAME C");
    
    profPanel_p->cd(6);
    histPP2_p[6]->DrawCopy("HIST SAME C");
    
    profPanel_p->cd(7);
    histPP3_p[6]->DrawCopy("HIST SAME C");
    
    profPanel_p->cd(8);
    histPP4_p[6]->DrawCopy("HIST SAME C");
  }
    
  histPPFile_p->Close();
  histPbPbFile_p->Close();

  delete histPPFile_p;
  delete histPbPbFile_p;

  histPbPbFile_p = new TFile(filePbPbName, "READ");
  histPPFile_p = new TFile(filePPName, "READ");

  grabHistForStackPbPbCent(histPbPbFile_p, gorr, algType[setNum], projMult, CNCR, Corr, cent, fileTagPbPb, histPbPb1_p);
  grabHistForStackPbPbCent(histPbPbFile_p, gorr, algType[setNum+1], projMult, CNCR, Corr, cent, fileTagPbPb, histPbPb2_p);
  grabHistForStackPbPbCent(histPbPbFile_p, gorr, algType[setNum+2], projMult, CNCR, Corr, cent, fileTagPbPb, histPbPb3_p);
  grabHistForStackPbPbCent(histPbPbFile_p, gorr, algType[setNum+3], projMult, CNCR, Corr, cent, fileTagPbPb, histPbPb4_p);
  
  grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum], projMult, CNCR, Corr, fileTagPP, histPP1_p);
  grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum+1], projMult, CNCR, Corr, fileTagPP, histPP2_p);
  grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum+2], projMult, CNCR, Corr, fileTagPP, histPP3_p);
  grabHistForStackPP(histPPFile_p, gorr, algTypePP[setNum+3], projMult, CNCR, Corr, fileTagPP, histPP4_p);

  for(Int_t iter = 0; iter < 7; iter++){
    histPbPb1_p[iter]->Add(histPP1_p[iter], -1);
    histPbPb2_p[iter]->Add(histPP2_p[iter], -1);
    histPbPb3_p[iter]->Add(histPP3_p[iter], -1);
    histPbPb4_p[iter]->Add(histPP4_p[iter], -1);
  }

  histToCanvasR(profPanel_p, 9, histPbPb1_p, histPbPb2_p, histPbPb3_p, histPbPb4_p, setNum, CNCR, projMult, montecarlo, cent, true);

  for(Int_t iter = 0; iter < 4; iter++){
    profPanel_p->cd(iter+9);
    label1_p->DrawLatex(xCoord1[iter], .92, Form("%s R=0.%d", dataMCStr1[2].c_str(), iter+2));
  }

  profPanel_p->cd(9);
  legA_p->Draw("SAME");
  profPanel_p->cd(10);
  legD_p->Draw("SAME");
  profPanel_p->cd(11);


  if(isInc(CNCR)) label1_p->DrawLatex(xCoord1[1], .35, "A_{J} Inclusive");
  else if(isU(CNCR)) label1_p->DrawLatex(xCoord1[1], .35, "A_{J}>0.22");
  else if(isD(CNCR)) label1_p->DrawLatex(xCoord1[1], .35, "A_{J}<0.22");
  label1_p->DrawLatex(xCoord1[1], .25, puVsString[setNum]);


  //Quark gluon label etc. here
  profPanel_p->cd(12);
  label1_p->DrawLatex(xCoord1[1], .25, "3rd Jet p_{T}>50 GeV/c");

  legB_p->AddEntry(dumHist_p[1], Form("%s", dataMCStr1[0].c_str()), "p");
  legB_p->AddEntry(dumHist_p[2], Form("%s", dataMCStr1[1].c_str()), "p");
  if(isAll(CNCR)) legC_p->AddEntry(dumLineHist_p, Form("%s cumulative", dataMCStr1[0].c_str()), "L");
  legB_p->AddEntry(dumHist_p[3], Form("%s", dataMCStr1[2].c_str()), "p");

  //  legA_p->AddEntry(dumHist_p[0], "0.0-0.5", "F");

  legA_p->AddEntry(dumHist_p[1], "0.5-1.0", "F");
  legA_p->AddEntry(dumHist_p[2], "1.0-2.0", "F");
  legD_p->AddEntry(dumHist_p[3], "2.0-4.0", "F");
  legD_p->AddEntry(dumHist_p[4], "4.0-8.0", "F");
  legD_p->AddEntry(dumHist_p[5], "8.0-300.0", "F");


  outFile_p->cd();
  profPanel_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(profPanel_p, Form("pdfDir/%s%s%sA%s%s%sPTStack_RDEP%s_%s", gorr, algType[setNum], projMult.c_str(), CNCR, Corr, "", cent.c_str(), fileTagPbPb), "pdf");

  delete profPanel_p;
  histPPFile_p->Close();
  delete histPPFile_p;

  histPbPbFile_p->Close();
  delete histPbPbFile_p;

  outFile_p->Close();
  delete outFile_p;

  return;
}


void makeDiJetPlots(const char* filePbPbName, const char* fileTagPbPb, const char* outName, Bool_t montecarlo = false, const char* filePPName = "", const char* fileTagPP = "", Bool_t doSpill = false, const char* filePbPbName2 = "", const char* fileTagPbPb2 = "", const char* filePPName2 = "", const char* fileTagPP2 = "", Bool_t isHighPtTrk = false)
{
  TH1::SetDefaultSumw2();

  Int_t jetAlgMax = 7;
  
  if(montecarlo)
    jetAlgMax = 8;
  
  const char* corr[2] = {"", "Corr"};
  const char* CNCR[52] = {"", "C", "NC", "C0", "C1", "C2", "C3", "R", "RD", "RU", "R_320", "RD_320", "RU_320", "R_340", "RD_340", "RU_340", "Eta", "EtaD", "EtaU", "Phi", "PhiD", "PhiU", "RCut", "RCutD", "RCutU", "RCutEta", "RCutEtaD", "RCutEtaU", "RCutPhi", "RCutPhiD", "RCutPhiU", "EtaCut", "EtaCutD", "EtaCutU", "PhiCut", "PhiCutD", "PhiCutU", "RFOR", "RFORD", "RFORU", "RFORMID", "RFORMIDD", "RFORMIDU", "RFORFOR", "RFORFORD", "RFORFORU", "RMID", "RMIDD", "RMIDU", "R2", "R2U", "R2D"};
  const char* Tight[2] = {"", "Tight"};

  for(Int_t algIter = 3; algIter < nSumAlg; algIter++){
    if(!montecarlo && algIter == nSumAlg - 4) break;
    for(Int_t tightIter = 0; tightIter < 1; tightIter++){
      for(Int_t corrIter = 1; corrIter < 2; corrIter++){
	for(Int_t CNCRIter = 0; CNCRIter < 39; CNCRIter++){
	  if(CNCRIter == 3 || CNCRIter == 4 || CNCRIter == 5 || CNCRIter == 6) continue;

	  if(CNCRIter > 15) continue;

	  if(algIter < nSumAlg - 4) makeImbPtStack(filePbPbName, fileTagPbPb, outName, "r", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], montecarlo, filePPName, fileTagPP, Tight[tightIter], doSpill, filePbPbName2, fileTagPbPb2, filePPName2, fileTagPP2, isHighPtTrk);


	  if(algIter == 16 || algIter == 22 || algIter == 28 || algIter == 32){
	    if(!strcmp("", CNCR[CNCRIter]) || !strcmp("C", CNCR[CNCRIter]) || !strcmp("NC", CNCR[CNCRIter])){
	      makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "r", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "010", montecarlo, filePPName, fileTagPP);
	      if(!strcmp("", CNCR[CNCRIter])) makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "r", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "030", montecarlo, filePPName, fileTagPP);
	      makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "r", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "50100", montecarlo, filePPName, fileTagPP);
	      if(!strcmp("", CNCR[CNCRIter])) makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "r", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "30100", montecarlo, filePPName, fileTagPP);
	    }
	    else{
	      makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "r", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "030", montecarlo, filePPName, fileTagPP);
	      makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "r", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "30100", montecarlo, filePPName, fileTagPP);
	    }
	  }

	  if(algIter == 36){
	    if(!strcmp("", CNCR[CNCRIter]) || !strcmp("C", CNCR[CNCRIter]) || !strcmp("NC", CNCR[CNCRIter])){
	      makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "g", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "010", montecarlo, filePPName, fileTagPP);
	      if(!strcmp("", CNCR[CNCRIter])) makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "g", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "030", montecarlo, filePPName, fileTagPP);
	      makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "g", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "50100", montecarlo, filePPName, fileTagPP);
	      if(!strcmp("", CNCR[CNCRIter])) makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "g", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "30100", montecarlo, filePPName, fileTagPP);
	    }
	    else{
	      makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "g", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "030", montecarlo, filePPName, fileTagPP);
	      makeImbPtRStack(filePbPbName, fileTagPbPb, outName, "g", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], "30100", montecarlo, filePPName, fileTagPP);
	    }
	  }


	  //	    if(CNCRIter < 31 && CNCRIter > 0) makeImbPtStack(filePbPbName, fileTagPbPb, outName, "r", algIter, "Mult", corr[corrIter], CNCR[CNCRIter], montecarlo, filePPName, fileTagPP, Tight[tightIter], filePbPbName2, fileTagPbPb2, filePPName2, fileTagPP2, isHighPtTrk);
	

	  if(montecarlo && corrIter > 0 && (algIter == 5 || algIter == 11 || algIter == 17 || algIter == 19 || algIter == 23 || algIter == 29 || algIter == 31 || algIter == 33 ||  algIter >= nSumAlg-4)){
	    makeImbPtStack(filePbPbName, fileTagPbPb, outName, "g", algIter, "Proj", corr[corrIter], CNCR[CNCRIter], montecarlo, filePPName, fileTagPP, Tight[tightIter], false, filePbPbName2, fileTagPbPb2, filePPName2, fileTagPP2, isHighPtTrk);

	    //	    if(CNCRIter < 31) makeImbPtStack(filePbPbName, fileTagPbPb, outName, "g", alvgIter, "Mult", corr[corrIter], CNCR[CNCRIter], montecarlo, filePPName, fileTagPP, Tight[tightIter], filePbPbName2, fileTagPbPb2, filePPName2, fileTagPP2, isHighPtTrk);
	  }
	  
	
	}
      }
    
      //      makeMultStack(filePbPbName, fileTagPbPb, outName, algIter, montecarlo, filePPName, fileTagPP);
      //      makeMultStack(filePbPbName, fileTagPbPb, outName, algIter, montecarlo, filePPName, fileTagPP, "Tight");

    }

  }
  return;
}


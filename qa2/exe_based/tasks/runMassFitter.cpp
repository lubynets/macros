// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file runMassFitter.C
/// \brief HFInvMassFitter class steering macro
///
/// \author Zhen Zhang <zhenz@cern.ch>
/// \author Mingyu Zhang <mingyu.zang@cern.ch>
/// \author Xinye Peng  <xinye.peng@cern.ch>
/// \author Biao Zhang <biao.zhang@cern.ch>

#if !defined(__CINT__) || defined(__CLING__)

#include "HFInvMassFitter.hpp"

#include <iostream> // std::cout
#include <string>   // std::string
#include <vector>   // std::vector

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>

// if .h file not found, please include your local rapidjson/document.h and rapidjson/filereadstream.h here
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#endif

using namespace rapidjson;

int runMassFitter(const TString& configFileName = "config_massfitter.json");

template <typename ValueType>
void readArray(const Value& jsonArray, std::vector<ValueType>& output)
{
  for (auto it = jsonArray.Begin(); it != jsonArray.End(); it++) {
    auto value = it->template Get<ValueType>();
    output.emplace_back(value);
  }
}

void parseStringArray(const Value& jsonArray, std::vector<std::string>& output)
{
  size_t arrayLength = jsonArray.Size();
  for (size_t i = 0; i < arrayLength; i++) {
    if (jsonArray[i].IsString()) {
      output.emplace_back(jsonArray[i].GetString());
    }
  }
}

void divideCanvas(TCanvas* c, int nSliceVarBins);
void setHistoStyle(TH1* histo, Color_t color = kBlack, Size_t markerSize = 1);

int runMassFitter(const TString& configFileName)
{
  // load config
  FILE* configFile = fopen(configFileName.Data(), "r");
  if (!configFile) {
    std::cerr << "ERROR: Missing configuration json file: " << configFileName << std::endl;
    return -1;
  }

  Document config;
  char readBuffer[65536];
  FileReadStream is(configFile, readBuffer, sizeof(readBuffer));
  config.ParseStream(is);
  fclose(configFile);

  Bool_t isMc = config["IsMC"].GetBool();
  TString inputFileName = config["InFileName"].GetString();
  TString reflFileName = config["ReflFileName"].GetString();
  TString outputFileName = config["OutFileName"].GetString();
  TString particleName = config["Particle"].GetString();

  std::vector<std::string> inputHistoName;
  std::vector<std::string> promptHistoName;
  std::vector<std::string> fdHistoName;
  std::vector<std::string> reflHistoName;
  std::vector<std::string> promptSecPeakHistoName;
  std::vector<std::string> fdSecPeakHistoName;
  TString sliceVarName;
  TString sliceVarUnit;
  std::vector<double> sliceVarMin;
  std::vector<double> sliceVarMax;
  std::vector<double> massMin;
  std::vector<double> massMax;
  std::vector<double> fixSigmaManual;
  std::vector<int> nRebin;
  std::vector<int> bkgFuncConfig;
  std::vector<int> sgnFuncConfig;

  const Value& inputHistoNameValue = config["InputHistoName"];
  parseStringArray(inputHistoNameValue, inputHistoName);

  const Value& promptHistoNameValue = config["PromptHistoName"];
  parseStringArray(promptHistoNameValue, promptHistoName);

  const Value& fdHistoNameValue = config["FDHistoName"];
  parseStringArray(fdHistoNameValue, fdHistoName);

  const Value& reflHistoNameValue = config["ReflHistoName"];
  parseStringArray(reflHistoNameValue, reflHistoName);

  const Value& promptSecPeakHistoNameValue = config["PromptSecPeakHistoName"];
  parseStringArray(promptSecPeakHistoNameValue, promptSecPeakHistoName);

  const Value& fdSecPeakHistoNameValue = config["FDSecPeakHistoName"];
  parseStringArray(fdSecPeakHistoNameValue, fdSecPeakHistoName);

  bool fixSigma = config["FixSigma"].GetBool();
  std::string sigmaFile = config["SigmaFile"].GetString();

  bool fixMean = config["FixMean"].GetBool();
  std::string meanFile = config["MeanFile"].GetString();

  const Value& fixSigmaManualValue = config["FixSigmaManual"];
  readArray(fixSigmaManualValue, fixSigmaManual);
  
  sliceVarName = config["SliceVarName"].GetString();
  sliceVarUnit = config["SliceVarUnit"].GetString();
  
  const Value& sliceVarMinValue = config["SliceVarMin"];
  readArray(sliceVarMinValue, sliceVarMin);

  const Value& sliceVarMaxValue = config["SliceVarMax"];
  readArray(sliceVarMaxValue, sliceVarMax);

  const Value& massMinValue = config["MassMin"];
  readArray(massMinValue, massMin);

  const Value& massMaxValue = config["MassMax"];
  readArray(massMaxValue, massMax);

  const Value& rebinValue = config["Rebin"];
  readArray(rebinValue, nRebin);

  bool includeSecPeak = config["InclSecPeak"].GetBool();
  bool useLikelihood = config["UseLikelihood"].GetBool();

  const Value& bkgFuncValue = config["BkgFunc"];
  readArray(bkgFuncValue, bkgFuncConfig);

  const Value& sgnFuncValue = config["SgnFunc"];
  readArray(sgnFuncValue, sgnFuncConfig);

  bool enableRefl = config["EnableRefl"].GetBool();

  const unsigned int nSliceVarBins = sliceVarMin.size();
  int bkgFunc[nSliceVarBins], sgnFunc[nSliceVarBins];
  double sliceVarLimits[nSliceVarBins + 1];

  for (unsigned int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    sliceVarLimits[iSliceVar] = sliceVarMin[iSliceVar];
    sliceVarLimits[iSliceVar + 1] = sliceVarMax[iSliceVar];

    if (bkgFuncConfig[iSliceVar] == 0) {
      bkgFunc[iSliceVar] = HFInvMassFitter::Expo;
    } else if (bkgFuncConfig[iSliceVar] == 1) {
      bkgFunc[iSliceVar] = HFInvMassFitter::Poly1;
    } else if (bkgFuncConfig[iSliceVar] == 2) {
      bkgFunc[iSliceVar] = HFInvMassFitter::Poly2;
    } else if (bkgFuncConfig[iSliceVar] == 3) {
      bkgFunc[iSliceVar] = HFInvMassFitter::Pow;
    } else if (bkgFuncConfig[iSliceVar] == 4) {
      bkgFunc[iSliceVar] = HFInvMassFitter::PowExpo;
    } else if (bkgFuncConfig[iSliceVar] == 5) {
      bkgFunc[iSliceVar] = HFInvMassFitter::Poly3;
    } else if (bkgFuncConfig[iSliceVar] == 6) {
      bkgFunc[iSliceVar] = HFInvMassFitter::NoBkg;
    } else {
      std::cerr << "ERROR: only Expo, Poly1, Poly2, Pow and PowEx background "
              "functions supported! Exit"
           << std::endl;
      return -1;
    }

    if (sgnFuncConfig[iSliceVar] == 0) {
      sgnFunc[iSliceVar] = HFInvMassFitter::SingleGaus;
    } else if (sgnFuncConfig[iSliceVar] == 1) {
      sgnFunc[iSliceVar] = HFInvMassFitter::DoubleGaus;
    } else if (sgnFuncConfig[iSliceVar] == 2) {
      sgnFunc[iSliceVar] = HFInvMassFitter::DoubleGausSigmaRatioPar;
    } else {
      std::cerr << "ERROR: only SingleGaus, DoubleGaus and DoubleGausSigmaRatioPar signal "
              "functions supported! Exit"
           << std::endl;
      return -1;
    }
  }

  TString massAxisTitle = "";
  double massPDG;
  if (particleName == "Dplus") {
    massAxisTitle = "#it{M}(K#pi#pi) (GeV/#it{c}^{2})";
    massPDG = TDatabasePDG::Instance()->GetParticle("D+")->Mass();
  } else if (particleName == "D0") {
    massAxisTitle = "#it{M}(K#pi) (GeV/#it{c}^{2})";
    massPDG = TDatabasePDG::Instance()->GetParticle("D0")->Mass();
  } else if (particleName == "Ds") {
    massAxisTitle = "#it{M}(KK#pi) (GeV/#it{c}^{2})";
    massPDG = TDatabasePDG::Instance()->GetParticle("D_s+")->Mass();
  } else if (particleName == "LcToPKPi") {
    massAxisTitle = "#it{M}(pK#pi) (GeV/#it{c}^{2})";
    massPDG = TDatabasePDG::Instance()->GetParticle("Lambda_c+")->Mass();
  } else if (particleName == "LcToPK0s") {
    massAxisTitle = "#it{M}(pK^{0}_{s}) (GeV/#it{c}^{2})";
    massPDG = TDatabasePDG::Instance()->GetParticle("Lambda_c+")->Mass();
  } else if (particleName == "Dstar") {
    massAxisTitle = "#it{M}(pi^{+}) (GeV/#it{c}^{2})";
    massPDG = TDatabasePDG::Instance()->GetParticle("D*+")->Mass();
  } else {
    std::cerr << "ERROR: only Dplus, D0, Ds, LcToPKPi, LcToPK0s and Dstar particles supported! Exit" << std::endl;
    return -1;
  }

  // load inv-mass histograms
  auto inputFile = TFile::Open(inputFileName.Data());
  if (!inputFile || !inputFile->IsOpen()) {
    return -1;
  }

  TFile* inputFileRefl = nullptr;
  if (enableRefl) {
    inputFileRefl = TFile::Open(reflFileName.Data());
    if (!inputFileRefl || !inputFileRefl->IsOpen()) {
      return -1;
    }
  }

  TH1* hMassSgn[nSliceVarBins];
  TH1* hMassRefl[nSliceVarBins];
  TH1* hMass[nSliceVarBins];

  for (unsigned int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    if (!isMc) {
      hMass[iSliceVar] = inputFile->Get<TH1>(inputHistoName[iSliceVar].data());
      if (enableRefl) {
        hMassRefl[iSliceVar] = inputFileRefl->Get<TH1>(reflHistoName[iSliceVar].data());
        hMassSgn[iSliceVar] = inputFileRefl->Get<TH1>(fdHistoName[iSliceVar].data());
        hMassSgn[iSliceVar]->Add(inputFileRefl->Get<TH1>(promptHistoName[iSliceVar].data()));
        if (!hMassRefl[iSliceVar]) {
          std::cerr << "ERROR: MC reflection histogram not found! Exit!" << std::endl;
          return -1;
        }
        if (!hMassSgn[iSliceVar]) {
          std::cerr << "ERROR: MC prompt or FD histogram not found! Exit!" << std::endl;
          return -1;
        }
      }
    } else {
      hMass[iSliceVar] = inputFile->Get<TH1>(promptHistoName[iSliceVar].data());
      hMass[iSliceVar]->Add(inputFile->Get<TH1>(fdHistoName[iSliceVar].data()));
      if (includeSecPeak) {
        hMass[iSliceVar]->Add(inputFile->Get<TH1>(promptSecPeakHistoName[iSliceVar].data()));
        hMass[iSliceVar]->Add(inputFile->Get<TH1>(fdSecPeakHistoName[iSliceVar].data()));
      }
    }
    if (!hMass[iSliceVar]) {
      std::cerr << "ERROR: input histogram for fit not found! Exit!" << std::endl;
      return -1;
    }
    hMass[iSliceVar]->SetDirectory(nullptr);
  }
  inputFile->Close();

  // define output histos
  auto hRawYields = new TH1D("hRawYields", ";" + sliceVarName + "(" + sliceVarUnit + ");raw yield",
                             nSliceVarBins, sliceVarLimits);
  auto hRawYieldsCounted = new TH1D("hRawYieldsCounted", ";" + sliceVarName + "(" + sliceVarUnit + ");raw yield via bin count",
                             nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSigma = new TH1D(
    "hRawYieldsSigma", ";" + sliceVarName + "(" + sliceVarUnit + ");width (GeV/#it{c}^{2})",
    nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSigmaRatio = new TH1D(
    "hRawYieldsSigmaRatio",
    ";" + sliceVarName + "(" + sliceVarUnit + ");ratio #sigma_{1}/#sigma_{2}", nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSigma2 = new TH1D(
    "hRawYieldsSigma2", ";" + sliceVarName + "(" + sliceVarUnit + ");width (GeV/#it{c}^{2})",
    nSliceVarBins, sliceVarLimits);
  auto hRawYieldsMean = new TH1D(
    "hRawYieldsMean", ";" + sliceVarName + "(" + sliceVarUnit + ");mean (GeV/#it{c}^{2})",
    nSliceVarBins, sliceVarLimits);
  auto hRawYieldsFracGaus2 = new TH1D(
    "hRawYieldsFracGaus2",
    ";" + sliceVarName + "(" + sliceVarUnit + ");second-gaussian fraction", nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSignificance = new TH1D(
    "hRawYieldsSignificance",
    ";" + sliceVarName + "(" + sliceVarUnit + ");significance (3#sigma)", nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSgnOverBkg =
    new TH1D("hRawYieldsSgnOverBkg", ";" + sliceVarName + "(" + sliceVarUnit + ");S/B (3#sigma)",
             nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSignal =
    new TH1D("hRawYieldsSignal", ";" + sliceVarName + "(" + sliceVarUnit + ");Signal (3#sigma)",
             nSliceVarBins, sliceVarLimits);
  auto hRawYieldsBkg =
    new TH1D("hRawYieldsBkg", ";" + sliceVarName + "(" + sliceVarUnit + ");Background (3#sigma)",
             nSliceVarBins, sliceVarLimits);
  auto hRawYieldsChiSquare =
    new TH1D("hRawYieldsChiSquare",
             ";" + sliceVarName + "(" + sliceVarUnit + ");#chi^{2}/#it{ndf}", nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSecondPeak = new TH1D(
    "hRawYieldsSecondPeak", ";" + sliceVarName + "(" + sliceVarUnit + ");raw yield second peak",
    nSliceVarBins, sliceVarLimits);
  auto hRawYieldsMeanSecondPeak =
    new TH1D("hRawYieldsMeanSecondPeak",
             ";" + sliceVarName + "(" + sliceVarUnit + ");mean second peak (GeV/#it{c}^{2})",
             nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSigmaSecondPeak =
    new TH1D("hRawYieldsSigmaSecondPeak",
             ";" + sliceVarName + "(" + sliceVarUnit + ");width second peak (GeV/#it{c}^{2})",
             nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSignificanceSecondPeak =
    new TH1D("hRawYieldsSignificanceSecondPeak",
             ";" + sliceVarName + "(" + sliceVarUnit + ");signficance second peak (3#sigma)",
             nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSigmaRatioSecondFirstPeak =
    new TH1D("hRawYieldsSigmaRatioSecondFirstPeak",
             ";" + sliceVarName + "(" + sliceVarUnit + ");width second peak / width first peak",
             nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSoverBSecondPeak = new TH1D(
    "hRawYieldsSoverBSecondPeak",
    ";" + sliceVarName + "(" + sliceVarUnit + ");S/B second peak (3#sigma)", nSliceVarBins, sliceVarLimits);
  auto hRawYieldsSignalSecondPeak = new TH1D(
    "hRawYieldsSignalSecondPeak",
    ";" + sliceVarName + "(" + sliceVarUnit + ");Signal second peak (3#sigma)", nSliceVarBins, sliceVarLimits);
  auto hRawYieldsBkgSecondPeak =
    new TH1D("hRawYieldsBkgSecondPeak",
             ";" + sliceVarName + "(" + sliceVarUnit + ");Background second peak (3#sigma)",
             nSliceVarBins, sliceVarLimits);
  auto hReflectionOverSignal =
    new TH1D("hReflectionOverSignal", ";" + sliceVarName + "(" + sliceVarUnit + ");Refl/Signal",
             nSliceVarBins, sliceVarLimits);

  const Int_t nConfigsToSave = 6;
  auto hFitConfig = new TH2F("hfitConfig", "Fit Configurations", nConfigsToSave, 0, 6, nSliceVarBins, sliceVarLimits);
  const char* hFitConfigXLabel[nConfigsToSave] = {"mass min", "mass max", "rebin num", "fix sigma", "bkg func", "sgn func"};
  hFitConfig->SetStats(0);
  hFitConfig->LabelsDeflate("X");
  hFitConfig->LabelsDeflate("Y");
  hFitConfig->LabelsOption("v");
  for (int i = 0; i < nConfigsToSave; i++) {
    hFitConfig->GetXaxis()->SetBinLabel(i + 1, hFitConfigXLabel[i]);
  }

  setHistoStyle(hRawYields);
  setHistoStyle(hRawYieldsCounted);
  setHistoStyle(hRawYieldsSigma);
  setHistoStyle(hRawYieldsSigma2);
  setHistoStyle(hRawYieldsMean);
  setHistoStyle(hRawYieldsFracGaus2);
  setHistoStyle(hRawYieldsSignificance);
  setHistoStyle(hRawYieldsSgnOverBkg);
  setHistoStyle(hRawYieldsSignal);
  setHistoStyle(hRawYieldsBkg);
  setHistoStyle(hRawYieldsChiSquare);
  setHistoStyle(hRawYieldsSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsMeanSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsSigmaSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsSignificanceSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsSigmaRatioSecondFirstPeak, kRed + 1);
  setHistoStyle(hRawYieldsSoverBSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsSignalSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsBkgSecondPeak, kRed + 1);
  setHistoStyle(hReflectionOverSignal, kRed + 1);

  TH1* hSigmaToFix = nullptr;
  if (fixSigma) {
    if (fixSigmaManual.empty()) {
      auto inputFileSigma = TFile::Open(sigmaFile.data());
      if (!inputFileSigma) {
        return -2;
      }
      hSigmaToFix = inputFileSigma->Get<TH1>("hRawYieldsSigma");
      hSigmaToFix->SetDirectory(0);
      if (static_cast<unsigned int>(hSigmaToFix->GetNbinsX()) != nSliceVarBins) {
        std::cout << "WARNING: Different number of bins for this analysis and histo for fix sigma!" << std::endl;
      }
      inputFileSigma->Close();
    }
  }

  TH1* hMeanToFix = nullptr;
  if (fixMean) {
    auto inputFileMean = TFile::Open(meanFile.data());
    if (!inputFileMean) {
      return -3;
    }
    hMeanToFix = inputFileMean->Get<TH1>("hRawYieldsMean");
    hMeanToFix->SetDirectory(0);
    if (static_cast<unsigned int>(hMeanToFix->GetNbinsX()) != nSliceVarBins) {
      std::cout << "WARNING: Different number of bins for this analysis and histo for fix mean" << std::endl;
    }
    inputFileMean->Close();
  }

  // fit histograms

  TH1* hMassForFit[nSliceVarBins];
  TH1* hMassForRefl[nSliceVarBins];
  TH1* hMassForSgn[nSliceVarBins];

  Int_t canvasSize[2] = {1920, 1080};
  if (nSliceVarBins == 1) {
    canvasSize[0] = 500;
    canvasSize[1] = 500;
  }

  Int_t nCanvasesMax = 20; // do not put more than 20 bins per canvas to make them visible
  const Int_t nCanvases = ceil(static_cast<float>(nSliceVarBins) / nCanvasesMax);
  TCanvas *canvasMass[nCanvases], *canvasResiduals[nCanvases], *canvasRefl[nCanvases];
  for (int iCanvas = 0; iCanvas < nCanvases; iCanvas++) {
    int nPads = (nCanvases == 1) ? nSliceVarBins : nCanvasesMax;
    canvasMass[iCanvas] = new TCanvas(Form("canvasMass%d", iCanvas), Form("canvasMass%d", iCanvas),
                                      canvasSize[0], canvasSize[1]);
    divideCanvas(canvasMass[iCanvas], nPads);

    canvasResiduals[iCanvas] =
      new TCanvas(Form("canvasResiduals%d", iCanvas), Form("canvasResiduals%d", iCanvas), canvasSize[0], canvasSize[1]);
    divideCanvas(canvasResiduals[iCanvas], nPads);
    canvasRefl[iCanvas] = new TCanvas(Form("canvasRefl%d", iCanvas), Form("canvasRefl%d", iCanvas),
                                      canvasSize[0], canvasSize[1]);
    divideCanvas(canvasRefl[iCanvas], nPads);
  }

  for (unsigned int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    Int_t iCanvas = floor(static_cast<float>(iSliceVar) / nCanvasesMax);

    hMassForFit[iSliceVar] = static_cast<TH1*>(hMass[iSliceVar]->Rebin(nRebin[iSliceVar]));
    TString ptTitle =
      Form("%0.2f < " + sliceVarName + " < %0.2f " + sliceVarUnit, sliceVarMin[iSliceVar], sliceVarMax[iSliceVar]);
    hMassForFit[iSliceVar]->SetTitle(Form("%s;%s;Counts per %0.1f MeV/#it{c}^{2}",
                                    ptTitle.Data(), massAxisTitle.Data(),
                                    hMassForFit[iSliceVar]->GetBinWidth(1) * 1000));
    hMassForFit[iSliceVar]->SetName(Form("MassForFit%d", iSliceVar));

    if (enableRefl) {
      hMassForRefl[iSliceVar] = static_cast<TH1*>(hMassRefl[iSliceVar]->Rebin(nRebin[iSliceVar]));
      hMassForSgn[iSliceVar] = static_cast<TH1*>(hMassSgn[iSliceVar]->Rebin(nRebin[iSliceVar]));
    }

    Double_t reflOverSgn = 0;
    double markerSize = 1.;
    if (nSliceVarBins > 15) {
      markerSize = 0.5;
    }

    if (isMc) {
      HFInvMassFitter* massFitter;
      massFitter = new HFInvMassFitter(hMassForFit[iSliceVar], massMin[iSliceVar], massMax[iSliceVar], HFInvMassFitter::NoBkg, sgnFunc[iSliceVar]);
      massFitter->setInitialGaussianMean(massPDG);
      massFitter->setParticlePdgMass(massPDG);
      massFitter->setBoundGaussianMean(massPDG, 0.8*massPDG, 1.2*massPDG);
      massFitter->doFit();

      if (nSliceVarBins > 1) {
        canvasMass[iCanvas]->cd(iSliceVar - nCanvasesMax * iCanvas + 1);
      } else {
        canvasMass[iCanvas]->cd();
      }

      massFitter->drawFit(gPad);

      Double_t rawYield = massFitter->getRawYield();
      Double_t rawYieldErr = massFitter->getRawYieldError();
      Double_t rawYieldCounted = massFitter->getRawYieldCounted();
      Double_t rawYieldCountedErr = massFitter->getRawYieldCountedError();

      Double_t sigma = massFitter->getSigma();
      Double_t sigmaErr = massFitter->getSigmaUncertainty();
      Double_t mean = massFitter->getMean();
      Double_t meanErr = massFitter->getMeanUncertainty();
      Double_t reducedChiSquare = massFitter->getChiSquareOverNDF();

      hRawYields->SetBinContent(iSliceVar + 1, rawYield);
      hRawYields->SetBinError(iSliceVar + 1, rawYieldErr);
      hRawYieldsCounted->SetBinContent(iSliceVar + 1, rawYieldCounted);
      hRawYieldsCounted->SetBinError(iSliceVar + 1, rawYieldCountedErr);
      hRawYieldsSigma->SetBinContent(iSliceVar + 1, sigma);
      hRawYieldsSigma->SetBinError(iSliceVar + 1, sigmaErr);
      hRawYieldsMean->SetBinContent(iSliceVar + 1, mean);
      hRawYieldsMean->SetBinError(iSliceVar + 1, meanErr);
      hRawYieldsChiSquare->SetBinContent(iSliceVar + 1, reducedChiSquare);
      hRawYieldsChiSquare->SetBinError(iSliceVar + 1, 0.);
    } else {
      HFInvMassFitter* massFitter;
      massFitter = new HFInvMassFitter(hMassForFit[iSliceVar], massMin[iSliceVar], massMax[iSliceVar],
                                       bkgFunc[iSliceVar], sgnFunc[iSliceVar]);
      massFitter->setInitialGaussianMean(massPDG);
      massFitter->setParticlePdgMass(massPDG);
      massFitter->setBoundGaussianMean(massPDG, 0.8*massPDG, 1.2*massPDG);
      if (useLikelihood) {
        massFitter->setUseLikelihoodFit();
      }
      if (fixMean) {
        massFitter->setFixGaussianMean(hMeanToFix->GetBinContent(iSliceVar + 1));
      }
      if (fixSigma) {
        if (fixSigmaManual.empty()) {
          massFitter->setFixGaussianSigma(hSigmaToFix->GetBinContent(iSliceVar + 1));
          std::cout << "*****************************"
               << "\n"
               << "FIXED SIGMA: " << hSigmaToFix->GetBinContent(iSliceVar + 1) << "\n"
               << "*****************************" << std::endl;
        } else if (!fixSigmaManual.empty()) {
          massFitter->setFixGaussianSigma(fixSigmaManual[iSliceVar]);
          std::cout << "*****************************"
               << "\n"
               << "FIXED SIGMA: " << fixSigmaManual[iSliceVar] << "\n"
               << "*****************************" << std::endl;
        } else {
          std::cout << "WARNING: impossible to fix sigma! Wrong fix sigma file or value!" << std::endl;
        }
      }

      if (enableRefl) {
        reflOverSgn = hMassForSgn[iSliceVar]->Integral(hMassForSgn[iSliceVar]->FindBin(massMin[iSliceVar] * 1.0001), hMassForSgn[iSliceVar]->FindBin(massMax[iSliceVar] * 0.999));
        reflOverSgn = hMassForRefl[iSliceVar]->Integral(hMassForRefl[iSliceVar]->FindBin(massMin[iSliceVar] * 1.0001), hMassForRefl[iSliceVar]->FindBin(massMax[iSliceVar] * 0.999)) / reflOverSgn;
        massFitter->setFixReflOverSgn(reflOverSgn);
        massFitter->setTemplateReflections(hMassRefl[iSliceVar], HFInvMassFitter::DoubleGaus);
      }

      massFitter->doFit();

      double rawYield = massFitter->getRawYield();
      double rawYieldErr = massFitter->getRawYieldError();
      double rawYieldCounted = massFitter->getRawYieldCounted();
      double rawYieldCountedErr = massFitter->getRawYieldCountedError();
      double sigma = massFitter->getSigma();
      double sigmaErr = massFitter->getSigmaUncertainty();
      double mean = massFitter->getMean();
      double meanErr = massFitter->getMeanUncertainty();
      double reducedChiSquare = massFitter->getChiSquareOverNDF();
      double significance = massFitter->getSignificance();
      double significanceErr = massFitter->getSignificanceError();
      double bkg = massFitter->getBkgYield();
      double bkgErr = massFitter->getBkgYieldError();

      hRawYields->SetBinContent(iSliceVar + 1, rawYield);
      hRawYields->SetBinError(iSliceVar + 1, rawYieldErr);
      hRawYieldsCounted->SetBinContent(iSliceVar + 1, rawYieldCounted);
      hRawYieldsCounted->SetBinError(iSliceVar + 1, rawYieldCountedErr);
      hRawYieldsSigma->SetBinContent(iSliceVar + 1, sigma);
      hRawYieldsSigma->SetBinError(iSliceVar + 1, sigmaErr);
      hRawYieldsMean->SetBinContent(iSliceVar + 1, mean);
      hRawYieldsMean->SetBinError(iSliceVar + 1, meanErr);
      hRawYieldsSignificance->SetBinContent(iSliceVar + 1, significance);
      hRawYieldsSignificance->SetBinError(iSliceVar + 1, significanceErr);
      hRawYieldsSgnOverBkg->SetBinContent(iSliceVar + 1, rawYield / bkg);
      hRawYieldsSgnOverBkg->SetBinError(iSliceVar + 1, rawYield / bkg * std::sqrt(rawYieldErr / rawYield * rawYieldErr / rawYield + bkgErr / bkg * bkgErr / bkg));
      hRawYieldsSignal->SetBinContent(iSliceVar + 1, rawYield);
      hRawYieldsSignal->SetBinError(iSliceVar + 1, rawYieldErr);
      hRawYieldsBkg->SetBinContent(iSliceVar + 1, bkg);
      hRawYieldsBkg->SetBinError(iSliceVar + 1, bkgErr);
      hRawYieldsChiSquare->SetBinContent(iSliceVar + 1, reducedChiSquare);
      hRawYieldsChiSquare->SetBinError(iSliceVar + 1, 1.e-20);
      if (enableRefl) {
        hReflectionOverSignal->SetBinContent(iSliceVar + 1, reflOverSgn);
      }

      if (enableRefl) {
        if (nSliceVarBins > 1) {
          canvasRefl[iCanvas]->cd(iSliceVar - nCanvasesMax * iCanvas + 1);
        } else {
          canvasRefl[iCanvas]->cd();
        }
        massFitter->drawReflection(gPad);
        canvasRefl[iCanvas]->Modified();
        canvasRefl[iCanvas]->Update();
      }

      if (nSliceVarBins > 1) {
        canvasMass[iCanvas]->cd(iSliceVar - nCanvasesMax * iCanvas + 1);
      } else {
        canvasMass[iCanvas]->cd();
      }
      massFitter->drawFit(gPad);
      canvasMass[iCanvas]->Modified();
      canvasMass[iCanvas]->Update();

      if (nSliceVarBins > 1) {
        canvasResiduals[iCanvas]->cd(iSliceVar - nCanvasesMax * iCanvas + 1);
      } else {
        canvasResiduals[iCanvas]->cd();
      }
      massFitter->drawResidual(gPad);
      canvasResiduals[iCanvas]->Modified();
      canvasResiduals[iCanvas]->Update();
    }

    hFitConfig->SetBinContent(1, iSliceVar + 1, massMin[iSliceVar]);
    hFitConfig->SetBinContent(2, iSliceVar + 1, massMax[iSliceVar]);
    hFitConfig->SetBinContent(3, iSliceVar + 1, nRebin[iSliceVar]);
    if (fixSigma) {
      if (fixSigmaManual.empty()) {
        hFitConfig->SetBinContent(4, iSliceVar + 1, hSigmaToFix->GetBinContent(iSliceVar + 1));
      } else {
        hFitConfig->SetBinContent(4, iSliceVar + 1, fixSigmaManual[iSliceVar]);
      }
    }
    hFitConfig->SetBinContent(5, iSliceVar + 1, bkgFuncConfig[iSliceVar]);
    hFitConfig->SetBinContent(6, iSliceVar + 1, sgnFuncConfig[iSliceVar]);
  }

  // save output histograms
  TFile outputFile(outputFileName.Data(), "recreate");
  for (int iCanvas = 0; iCanvas < nCanvases; iCanvas++) {
    canvasMass[iCanvas]->Write();
    if (!isMc) {
      canvasResiduals[iCanvas]->Write();
      canvasRefl[iCanvas]->Write();
    }
  }

  for (unsigned int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    hMass[iSliceVar]->Write();
  }
  hRawYields->Write();
  hRawYieldsCounted->Write();
  hRawYieldsSigma->Write();
  hRawYieldsMean->Write();
  hRawYieldsSignificance->Write();
  hRawYieldsSgnOverBkg->Write();
  hRawYieldsSignal->Write();
  hRawYieldsBkg->Write();
  hRawYieldsChiSquare->Write();
  hRawYieldsSigma2->Write();
  hRawYieldsFracGaus2->Write();
  hRawYieldsSecondPeak->Write();
  hRawYieldsMeanSecondPeak->Write();
  hRawYieldsSigmaSecondPeak->Write();
  hRawYieldsSignificanceSecondPeak->Write();
  hRawYieldsSigmaRatioSecondFirstPeak->Write();
  hRawYieldsSoverBSecondPeak->Write();
  hRawYieldsSignalSecondPeak->Write();
  hRawYieldsBkgSecondPeak->Write();
  hFitConfig->Write();

  outputFile.Close();

  outputFileName.ReplaceAll(".root", ".pdf");
  TString outputFileNameResidual = outputFileName;
  outputFileNameResidual.ReplaceAll(".pdf", "_Residuals.pdf");
  for (int iCanvas = 0; iCanvas < nCanvases; iCanvas++) {
    if (iCanvas == 0 && nCanvases > 1) {
      canvasMass[iCanvas]->SaveAs(Form("%s[", outputFileName.Data()));
    }
    canvasMass[iCanvas]->SaveAs(outputFileName.Data());
    if (iCanvas == nCanvases - 1 && nCanvases > 1) {
      canvasMass[iCanvas]->SaveAs(Form("%s]", outputFileName.Data()));
    }
    if (!isMc) {
      if (iCanvas == 0 && nCanvases > 1) {
        canvasResiduals[iCanvas]->SaveAs(Form("%s[", outputFileNameResidual.Data()));
      }
      canvasResiduals[iCanvas]->SaveAs(outputFileNameResidual.Data());
      if (iCanvas == nCanvases - 1 && nCanvases > 1) {
        canvasResiduals[iCanvas]->SaveAs(Form("%s]", outputFileNameResidual.Data()));
      }
    }
  }
  return 0;
}

void setHistoStyle(TH1* histo, Color_t color, Size_t markerSize)
{
  histo->SetStats(kFALSE);
  histo->SetMarkerSize(markerSize);
  histo->SetMarkerStyle(20);
  histo->SetLineWidth(2);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}

void divideCanvas(TCanvas* canvas, int nSliceVarBins)
{
  if (nSliceVarBins < 2) {
    canvas->cd();
  } else if (nSliceVarBins == 2 || nSliceVarBins == 3) {
    canvas->Divide(nSliceVarBins, 1);
  } else if (nSliceVarBins == 4 || nSliceVarBins == 6 || nSliceVarBins == 8) {
    canvas->Divide(nSliceVarBins / 2, 2);
  } else if (nSliceVarBins == 5 || nSliceVarBins == 7) {
    canvas->Divide((nSliceVarBins + 1) / 2, 2);
  } else if (nSliceVarBins == 9 || nSliceVarBins == 12 || nSliceVarBins == 15) {
    canvas->Divide(nSliceVarBins / 3, 3);
  } else if (nSliceVarBins == 10 || nSliceVarBins == 11) {
    canvas->Divide(4, 3);
  } else if (nSliceVarBins == 13 || nSliceVarBins == 14) {
    canvas->Divide(5, 3);
  } else if (nSliceVarBins > 15 && nSliceVarBins <= 20 && nSliceVarBins % 4 == 0) {
    canvas->Divide(nSliceVarBins / 4, 4);
  } else if (nSliceVarBins > 15 && nSliceVarBins <= 20 && nSliceVarBins % 4 != 0) {
    canvas->Divide(5, 4);
  } else if (nSliceVarBins == 21) {
    canvas->Divide(7, 3);
  } else if (nSliceVarBins > 21 && nSliceVarBins <= 25) {
    canvas->Divide(5, 5);
  } else if (nSliceVarBins > 25 && nSliceVarBins % 2 == 0) {
    canvas->Divide(nSliceVarBins / 2, 2);
  } else {
    canvas->Divide((nSliceVarBins + 1) / 2, 2);
  }
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./runMassFitter configFileName" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string configFileName = argv[1];

  runMassFitter(configFileName);

  return 0;
}
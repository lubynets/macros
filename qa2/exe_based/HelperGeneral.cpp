//
// Created by oleksii on 21.07.25.
//

#include "HelperGeneral.hpp"

#include <TH1.h>
#include <TROOT.h>

#include <fstream>
#include <string>
#include <string_view>

bool HelperGeneral::string_to_bool(const std::string& str) {
  if(str == "true") return true;
  else if(str == "false") return false;
  else throw std::runtime_error("string_to_bool(): argument must be either true or false");
}

void HelperGeneral::LoadMacro(const std::string& macroName) {
  TString currentMacroPath = __FILE__;
  TString directory = currentMacroPath(0, currentMacroPath.Last('/'));
  gROOT->Macro( directory + "/" + macroName );
}

void HelperGeneral::CD(TFile* file, const std::string& dirName) {
  if(file == nullptr) throw std::runtime_error("Helper::CD() - file is nullptr");

  if(file->GetDirectory(dirName.c_str()) == nullptr) file->mkdir(dirName.c_str());
  file->cd(dirName.c_str());
}

std::map<std::string_view, int> HelperGeneral::MapTHnSparseAxesIndices(const THnSparse* histo) {
  std::map<std::string_view, int> result;
  const int nDims = histo->GetNdimensions();
  for(int iDim=0; iDim<nDims; ++iDim) {
    result.insert({histo->GetAxis(iDim)->GetTitle(), iDim});
  }
  return result;
}

void HelperGeneral::CheckTAxisForRanges(const TAxis& axis, const std::vector<double>& ranges) {
  const int nBins = axis.GetNbins();
  for(const auto& range : ranges) {
    bool ok{false};
    for(int iBin=1; iBin<=nBins+1; ++iBin) {
      const float edge = axis.GetBinLowEdge(iBin);
      if(std::fabs(edge - range) < 1e-4) {
        ok = true;
        break;
      }
    }
    if(!ok) {
      throw std::runtime_error("HelperGeneral::CheckTAxisForRanges() - the range " + std::to_string(range) + " is missing");
    }
  }
}

void HelperGeneral::SetTHnSparseAxisRanges(THnSparse* histo, int axisNum, float lo, float hi) {
  constexpr double tolerance = 1e-6;

  if(std::fabs(lo+999)<tolerance && std::fabs(hi+999)<tolerance) {
    histo->GetAxis(axisNum)->SetRange();
    return;
  }

  if(lo >= hi) throw std::runtime_error("SetTHnSparseAxisRanges(): lo >= hi");

  const TAxis* axis = histo->GetAxis(axisNum);
  int binLo{-999}, binHi{-999};
  for(int iBin=1, nBins=axis->GetNbins(); iBin<=nBins; ++iBin) {
    const float binLowEdge = axis->GetBinLowEdge(iBin);
    const float binUpEdge = axis->GetBinUpEdge(iBin);
    if(std::fabs(binLowEdge - lo)<tolerance) binLo = iBin;
    if(std::fabs(binUpEdge - hi)<tolerance) binHi = iBin;
    if(binLo != -999 && binHi != -999) break;
  }
  if(binLo == -999 || binHi == -999) throw std::runtime_error("SetTHnSparseAxisRanges(): binLo == -999 || binHi == -999");
  histo->GetAxis(axisNum)->SetRange(binLo, binHi);
}

double HelperGeneral::InterpolateTH1SuppressWarning(const TH1* h, double value) {
  double result;
  if (value <= h->GetBinLowEdge(1) || value >= h->GetBinLowEdge(h->GetNbinsX() + 1)) result = 0.;
  else
    result = h->Interpolate(value);
  return result;
}

std::string HelperGeneral::ReadNthLine(const std::string& fileName) {
  if(fileName.find(':') == std::string::npos) return fileName;

  std::string result;
  const size_t colonPosition = fileName.find(':');
  const std::string fileListName = fileName.substr(0, colonPosition);
  const std::string fileLineNumberStr = fileName.substr(colonPosition + 1);
  const int fileLineNumberInt = std::stoi(fileLineNumberStr);

  std::ifstream fileList(fileListName);
  if (!fileList.is_open()) throw std::runtime_error("HelperGeneral::ReadNthLine() - the fileList " + fileListName + " is missing!");

  for(size_t iLine=0; iLine<fileLineNumberInt; ++iLine) {
    if (!std::getline(fileList, result)) throw std::runtime_error("HelperGeneral::ReadNthLine() - the EOF of fileList " + fileListName + " reached before line " + fileLineNumberStr);
  }

  return result;
}

void HelperGeneral::MkDirBash(const std::string& dirName) {
  const auto status = std::system(("mkdir -p " + dirName).c_str());
  if(status != 0) {
    throw std::runtime_error("HelperGeneral::MkDirBash() - could not create directory " + dirName);
  }
}

TFile* HelperGeneral::OpenFileWithNullptrCheck(const std::string& fileName, const std::string& option) {
  TFile* file = TFile::Open(fileName.c_str(), option.c_str());
  if(file == nullptr) {
    throw std::runtime_error("HelperGeneral::OpenFileWithNullptrCheck() - file " + fileName + " is missing");
  }
  return file;
}

void HelperGeneral::ReplaceSubstrInStr(std::string& s, const std::string& from, const std::string& to) {
  if (from.empty())
    return;

  std::size_t pos = 0;

  while ((pos = s.find(from, pos)) != std::string::npos) {
    s.replace(pos, from.size(), to);
    pos += to.size();
  }
}

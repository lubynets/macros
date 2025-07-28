//
// Created by oleksii on 21.07.25.
//

#include "HelperGeneral.hpp"

#include <TH1.h>
#include <TROOT.h>

#include <iostream>

std::vector<std::pair<std::string, std::string>> HelperGeneral::FindCuts(TFile* fileIn, std::string name_start, bool printCuts) {
  if (name_start.back() != '_') name_start.push_back('_');
  std::vector<std::pair<std::string, std::string>> result;

  auto getSubstringBeforeLastSlash = [](const std::string& input) {
      size_t pos = input.rfind('/');  // Find the last occurrence of the slash

      if (pos == std::string::npos) {
        return std::string();  // No slash found, return an empty string
      }

      return input.substr(0, pos);  // Return the substring before the last slash
  };

  const std::string name_start_before_slash = getSubstringBeforeLastSlash(name_start);
  auto lok = fileIn->GetDirectory(name_start_before_slash.c_str())->GetListOfKeys();
  const int nDirs = lok->GetEntries();
  for (int iDir = 0; iDir < nDirs; iDir++) {
    const std::string dirName = lok->At(iDir)->GetName();
    const std::string fullDirName = !name_start_before_slash.empty() ? name_start_before_slash + "/" + dirName : dirName;
    if (fullDirName.substr(0, name_start.size()) != name_start) continue;
    std::pair<std::string, std::string> cutPair;
    bool isFirstCutRead{false};
    bool isGoodCut{true};
    for (int iChar = name_start.size(); iChar < fullDirName.size(); iChar++) {
      char letter = fullDirName.at(iChar);
      if(!(std::isdigit(letter) || letter == '.' || letter == '_')) {
        isGoodCut = false;
        break;
      }
      if (letter != '_') {
        if (!isFirstCutRead) cutPair.first.push_back(letter);
        else cutPair.second.push_back(letter);
      } else {
        isFirstCutRead = true;
      }
    }
    if(isGoodCut) result.emplace_back(cutPair);
  }

  auto stofCompare = [](const std::pair<std::string, std::string>& a, const std::pair<std::string, std::string>& b) {
      return atof(a.first.c_str()) < atof(b.first.c_str());
  };

  std::sort(result.begin(), result.end(), stofCompare);

  if (result.empty()) {
    throw std::runtime_error("FindCuts(): " + name_start + " cuts are not present");
  }

  if(printCuts) {
    std::cout << "Slice cuts are:\n";
    for(auto& r : result) std::cout << r.first << "\t" << r.second << "\n";
  }

  return result;
}

bool HelperGeneral::string_to_bool(const std::string& str) {
  if(str == "true") return true;
  else if(str == "false") return false;
  else throw std::runtime_error("string_to_bool(): argument must be either true or false");
}

TFile* HelperGeneral::OpenFileWithNullptrCheck(const std::string& fileName, const std::string& option) {
  TFile* file = TFile::Open(fileName.c_str(), option.c_str());
  if(file == nullptr) {
    throw std::runtime_error("OpenFileWithNullptrCheck() - file " + fileName + " is missing");
  }
  return file;
}

void HelperGeneral::PrintInfoOnTF1(const TF1* f) {
  std::cout << "Name = " << f->GetName() << "\n";
  std::cout << "Title = " << f->GetTitle() << "\n";
  const int nPar = f->GetNpar();
  std::cout << "NPar = " << nPar << "\n";
  std::cout << "NFreePar = " << f->GetNumberFreeParameters() << "\n";
  for(int iPar = 0; iPar<nPar; iPar++) {
    double min, max;
    f->GetParLimits(iPar, min, max);
    std::cout << f->GetParName(iPar) << "\t" << f->GetParameter(iPar) << " +- " << f->GetParError(iPar) << "\t(" << min << "; " << max << ")\n";
  }
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

void HelperGeneral::CheckHistogramsForXaxisIdentity(const TH1* h1, const TH1* h2) {
  if(h1->GetNbinsX() != h2->GetNbinsX()) {
    throw std::runtime_error("HelperGeneral::CheckHistogramsForXaxisIdentity(): nBinsX do not match for " + static_cast<std::string>(h1->GetName()) + " and " + h2->GetName());
  }
  const int nBins = h1->GetNbinsX();
  for(int iBin=1; iBin<=nBins; iBin++) {
    if(std::abs(h1->GetBinCenter(iBin) - h2->GetBinCenter(iBin)) > 1e-6) {
      throw std::runtime_error("HelperGeneral::CheckHistogramsForXaxisIdentity(): bins do not coincide for " + static_cast<std::string>(h1->GetName()) + " and " + h2->GetName());
    }
  }
}
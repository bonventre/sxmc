#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TDirectory.h>
#include <TNtupleD.h>
#include <TMath.h>

#include <sxmc/utils.h>

unsigned nint(float nexpected) {
  return TMath::Nint(nexpected);
}


double get_ntuple_entry(TNtupleD* nt, int i, std::string field) {
  double v;
  nt->SetBranchAddress(field.c_str(), &v);
  assert(i < nt->GetEntries());
  nt->GetEvent(i);
  nt->ResetBranchAddresses();
  return v;
}


std::vector<double> get_correlation_matrix(TNtupleD* nt) {
  int nentries = nt->GetEntries();

  // Get list of branch names
  std::vector<std::string> names;
  for (int i=0; i<nt->GetListOfBranches()->GetEntries(); i++) {
    std::string name = nt->GetListOfBranches()->At(i)->GetName();
    if (name == "likelihood") {
      continue;
    }
    names.push_back(name);
  }

  std::vector<double> matrix(names.size() * names.size());

  // Convert the ntuple to a vector, calculating means as we go
  std::vector<double> table(names.size() * nentries);
  std::vector<double> means(names.size(), 0);
  for (int i=0; i<nentries; i++) {
    for (size_t j=0; j<names.size(); j++) {
      double v = get_ntuple_entry(nt, i, names.at(j));
      table.at(j + i * names.size()) = v;
      means.at(j) += v;
    }
  }

  // Sums to means
  for (size_t i=0; i<names.size(); i++) {
    means.at(i) /= nentries;
  }

  // Compute correlations
  for (size_t i=0; i<names.size(); i++) {
    for (size_t j=i; j<names.size(); j++) {
      double t = 0;
      double dx2 = 0;
      double dy2 = 0;
      for (int k=0; k<nentries; k++) {
        double x1 = table.at(i + k * names.size()) - means.at(i);
        double x2 = table.at(j + k * names.size()) - means.at(j);
        t += x1 * x2;
        dx2 += x1 * x1;
        dy2 += x2 * x2;
      }
      matrix.at(i * names.size() + j) = t / TMath::Sqrt(dx2 * dy2);
    }
  }

  return matrix;
}


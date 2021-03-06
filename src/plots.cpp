#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TNtuple.h>
#include <TLegend.h>

#include <sxmc/plots.h>

const int ncolors = 27;
const int colors[27] = {kRed,      kGreen,    kBlue,      kMagenta,
                        kCyan,     kYellow,   kOrange,    kViolet+2,
                        kRed+2,    kGreen+2,  kBlue+2,    kMagenta+2,
                        kCyan+2,   kYellow+2, kOrange+2,  kRed-7,
                        kGreen-7,  kBlue-7,   kMagenta-7, kCyan-7,
                        kYellow-7, kOrange-7, kRed-6,     kAzure+1,
                        kTeal+1,   kSpring-9, kAzure-9};



SpectralPlot::SpectralPlot(int _line_width, float _xmin, float _xmax,
                           float _ymin, float _ymax, bool _logy,
                           std::string _title,
                           std::string _xtitle,
                           std::string _ytitle)
    : logy(_logy), line_width(_line_width),
      xmin(_xmin), xmax(_xmax), ymin(_ymin), ymax(_ymax),
      title(_title), xtitle(_xtitle), ytitle(_ytitle) {
  this->c = new TCanvas();

  if (this->logy) {
    this->c->SetLogy();
  }

  this->legend = new TLegend(0.85, 0.15, 0.985, 0.95);
  this->legend->SetFillColor(kWhite);
}


SpectralPlot::SpectralPlot(const SpectralPlot& o) {
  this->logy = o.logy;
  this->line_width = o.line_width;
  this->xmin = o.xmin;
  this->xmax = o.xmax;
  this->ymin = o.ymin;
  this->ymax = o.ymax;
  this->title = o.title;
  this->xtitle = o.xtitle;
  this->ytitle = o.ytitle;
  for (size_t i=0; i<o.histograms.size(); i++) {
    this->histograms.push_back(o.histograms[i]);
  }
  this->c = new TCanvas();

  if (o.logy) {
    this->c->SetLogy();
  }
  this->legend = (TLegend*) o.legend->Clone("");
}


SpectralPlot::~SpectralPlot() {
  for (size_t i=0; i<this->histograms.size(); i++) {
    delete this->histograms[i];
  }
  delete this->legend;
}


void SpectralPlot::add(TH1* _h, std::string title, std::string options) {
  std::string name = "__" + std::string(title);
  TH1* h = dynamic_cast<TH1*>(_h->Clone(name.c_str()));
  h->SetDirectory(NULL);

  h->SetLineWidth(this->line_width);
  h->SetTitle(this->title.c_str());
  h->SetXTitle(xtitle.c_str());
  h->SetYTitle(ytitle.c_str());

  this->legend->AddEntry(h, title.c_str());

  if (!h || h->Integral() == 0) {
    return;
  }

  this->histograms.push_back(h);

  if (this->histograms.size() == 1) {
    h->SetAxisRange(this->ymin, this->ymax, "Y");
    h->SetAxisRange(this->xmin, this->xmax, "X");
    h->GetXaxis()->SetLabelFont(132);
    h->GetXaxis()->SetTitleFont(132);
    h->GetYaxis()->SetLabelFont(132);
    h->GetYaxis()->SetTitleFont(132);
    if (this->logy) {
      this->c->SetLogy();
    }
    this->c->cd();
    h->DrawClone(options.c_str());
  }
  else {
    this->c->cd();
    h->DrawClone(("same " + options).c_str());
  }
  this->c->Update();
}


void SpectralPlot::save(std::string filename) {
  this->c->cd();
  this->legend->SetTextFont(132);
  this->legend->Draw();
  this->c->Update();
  this->c->SaveAs(filename.c_str(), "q");
}


TH1* SpectralPlot::make_like(TH1* h, std::string name) {
  TH1* hnew = dynamic_cast<TH1*>(h->Clone(name.c_str()));
  hnew->Reset();
  return hnew;
}


void plot_fit(std::map<std::string, Interval> best_fit, float live_time,
              std::vector<Signal> signals,
              std::vector<Systematic> systematics,
              std::vector<Observable> observables,
              std::vector<float> data, std::vector<int> weights,
              std::string output_path) {

  std::vector<std::string> categories;
  categories.push_back("full");
  for (size_t i=0;i<signals.size();i++){
    if (signals[i].category.size() > 0){
      if (std::find(categories.begin(),categories.end(),signals[i].category) == categories.end())
        categories.push_back(signals[i].category);
    }
  }

  std::map<std::string, std::vector<SpectralPlot> > all_plots;
  std::map<std::string, std::vector<TH1D*> > all_totals;

  for (size_t j=0;j<categories.size();j++){
    std::vector<SpectralPlot> plots;

    // Set up plots for each observable
    for (size_t i=0; i<observables.size(); i++) {
      Observable* o = &observables[i];
      std::stringstream ytitle;
      ytitle << "Counts/" << std::setprecision(3)
        << (o->upper - o->lower) / o->bins << " " << o->units
        << "/" << live_time << " y";
      plots.push_back(SpectralPlot(2, o->lower, o->upper, 1e-2, 1e6,
            true, "", o->title, ytitle.str().c_str()));
    }
    all_plots[categories[j]] = plots;

    std::vector<TH1D*> totals(observables.size(), NULL);
    all_totals[categories[j]] = totals;
  }

//  std::vector<TH1D*> external_total(observables.size(), NULL);
//  std::vector<TH1D*> cosmogenic_total(observables.size(), NULL);
//  std::vector<TH1D*> fit_total(observables.size(), NULL);

  // Extract best-fit parameter values
  std::vector<float> params;
  for (size_t i=0; i<signals.size(); i++) {
    params.push_back(best_fit[signals[i].name].point_estimate);
  }
  for (size_t i=0; i<systematics.size(); i++) {
    params.push_back(best_fit[systematics[i].name].point_estimate);
  }

  hemi::Array<unsigned> norms_buffer(signals.size(), true);
  norms_buffer.writeOnlyHostPtr();

  hemi::Array<double> param_buffer(params.size(), true);
  for (size_t i=0; i<params.size(); i++) {
    param_buffer.writeOnlyHostPtr()[i] = params[i];
  }

  for (size_t i=0; i<signals.size(); i++) {
    pdfz::EvalHist* phist = \
      dynamic_cast<pdfz::EvalHist*>(signals[i].histogram);

    phist->SetParameterBuffer(&param_buffer, signals.size());
    phist->SetNormalizationBuffer(&norms_buffer, i);

    TH1* hpdf_nd = phist->CreateHistogram();
    hpdf_nd->Scale(params[i] / hpdf_nd->Integral());

    std::vector<TH1D*> hpdf(observables.size(), NULL);
    if (hpdf_nd->IsA() == TH1D::Class()) {
      hpdf[0] = dynamic_cast<TH1D*>(hpdf_nd);
    }
    else if (hpdf_nd->IsA() == TH2D::Class()) {
      hpdf[0] = dynamic_cast<TH2D*>(hpdf_nd)->ProjectionX("hpdf_x");
      hpdf[1] = dynamic_cast<TH2D*>(hpdf_nd)->ProjectionY("hpdf_y");
    }
    else if (hpdf_nd->IsA() == TH3D::Class()) {
      hpdf[0] = dynamic_cast<TH3D*>(hpdf_nd)->ProjectionX("hpdf_x");
      hpdf[1] = dynamic_cast<TH3D*>(hpdf_nd)->ProjectionY("hpdf_y");
      hpdf[2] = dynamic_cast<TH3D*>(hpdf_nd)->ProjectionZ("hpdf_y");
    }

    std::string n = signals[i].name;

    for (size_t j=0; j<observables.size(); j++) {
      hpdf[j]->SetLineColor(colors[i % ncolors]);
      if (all_totals["full"][j] == NULL) {
        std::string hfname = "fit_total_" + signals[i].name;
        all_totals["full"][j] = (TH1D*) hpdf[j]->Clone(hfname.c_str());
      } else {
        if (hpdf[j] && hpdf[j]->Integral() > 0) {
          all_totals["full"][j]->Add(hpdf[j]);
        }
      }

      if (signals[i].category.size() == 0){
        if (hpdf[j] && hpdf[j]->Integral() > 0) {
          all_plots["full"][j].add(hpdf[j], signals[i].title, "hist");
        }
      }else{
        for (size_t k=0;k<categories.size();k++){
          if (signals[i].category == categories[k]){
            all_plots[categories[k]][j].add(hpdf[j], signals[i].title, "hist");
            std::string hname = "et" + signals[i].name + observables[j].name;
            if (all_totals[categories[k]][j] == NULL) {
              all_totals[categories[k]][j] = (TH1D*) hpdf[j]->Clone(hname.c_str());
            } else {
              if (hpdf[j] && hpdf[j]->Integral() > 0) {
                all_totals[categories[k]][j]->Add((TH1D*) hpdf[j]->Clone(hname.c_str()));
              }
            }
          }
        } // end loop over categories
      }
    } // end loop over observables
  } // end loop over signals

  for (size_t i=0; i<observables.size(); i++) {
    TH1D* hdata = \
      (TH1D*) SpectralPlot::make_like(all_plots["full"][i].histograms[0], "hdata");

    hdata->SetMarkerStyle(20);
    hdata->SetLineColor(kBlack);

    for (size_t idata=0; idata<data.size() / observables.size(); idata++) {
      hdata->Fill(data[idata * observables.size() + i],weights[idata]);
    }

    for (size_t j=0;j<categories.size();j++){
      if (categories[j] == "full")
        continue;
      if (all_totals[categories[j]][i] != NULL){
        all_totals[categories[j]][i]->SetLineColor(colors[j+1]);
        TH1D* t = (TH1D*) all_totals[categories[j]][i]->Clone(categories[j].c_str());
        t->SetLineStyle(2);
        all_plots[categories[j]][i].add(t,"Total", "hist");
        all_plots["full"][i].add(all_totals[categories[j]][i],categories[j],"hist");
      }
    }

    if (all_totals["full"][i] != NULL) {
      all_totals["full"][i]->SetLineColor(kRed);
      all_plots["full"][i].add(all_totals["full"][i], "Fit", "hist");
    }

    all_plots["full"][i].add(hdata, "Fake Data");

    all_plots["full"][i].save(output_path + observables[i].name + "_spectrum_full.pdf");
    for (size_t j=0;j<categories.size();j++){
      if (categories[j] == "full")
        continue;
      all_plots[categories[j]][i].save(output_path + observables[i].name + "_spectrum_" + categories[j] + ".pdf");
    }
  }
}


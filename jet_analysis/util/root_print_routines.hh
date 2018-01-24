
#ifndef ROOT_PRINT_ROUTINES_HH
#define ROOT_PRINT_ROUTINES_HH

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>
#include <vector>
#include <string>

#include "boost/filesystem.hpp"

struct histogramOpts {
  
};

template<typename H>
void PrettyPrint1D(H* h,
                   std::string hist_title,
                   std::string output_loc,
                   std::string output_name,
                   std::string canvas_title,
                   std::string x_axis_label,
                   std::string y_axis_label,
                   bool logx,
                   bool logy,
                   bool do_legend,
                   std::string legend_title,
                   double y_max = 9999,
                   double y_min = -9999) {
  // we assume the output location exists, so create
  // the final output string that will be used for pdf creation
  std::string canvas_name = output_loc + "/" + output_name + ".pdf";
  
  // and axis labels, and title
  h->SetTitle(canvas_title.c_str());
  h->GetXaxis()->SetTitle(x_axis_label.c_str());
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetTitleSize(0.075);
  h->GetXaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->SetTitle(y_axis_label.c_str());
  h->GetYaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetTitleSize(0.075);
  h->GetYaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->CenterTitle(true);
  
  // generate a canvas
  TCanvas c;
  c.SetLeftMargin(0.12);
  c.SetBottomMargin(0.15);
  
  h->SetLineColor(1);
  h->SetMarkerColor(1);
  h->SetMarkerSize(1);
  h->SetLineWidth(2);
  h->SetMarkerStyle(21);
  
  h->Draw();
  
  if (do_legend) {
    TLegend* leg = new TLegend(0.65, 0.65, .88, .88);
    leg->SetTextSize(0.04);
    leg->SetHeader(legend_title.c_str());
    leg->AddEntry(h, hist_title.c_str(), "lep");
  }
  
  if (logx)
    c.SetLogx();
  if (logy)
    c.SetLogy();
  c.SaveAs(canvas_name.c_str());
}

template<typename H>
void Overlay1D(const std::vector<H*>& h,
               std::vector<std::string> hist_titles,
               std::string output_loc,
               std::string output_name,
               std::string canvas_title,
               std::string x_axis_label,
               std::string y_axis_label,
               bool logx,
               bool logy,
               bool do_legend,
               std::string legend_title,
               double y_max = 9999,
               double y_min = -9999) {
  
  // first, check that there is a name for each histogram
  if (h.size() != hist_titles.size()) {
    std::cerr << "incorrect number of histogram names for given set of histograms" << std::endl;
    std::cerr << "for canvas: " << canvas_title << ", exiting" << std::endl;
    return;
  }
  
  // we assume the output location exists, so create
  // the final output string that will be used for pdf creation
  std::string canvas_name = output_loc + "/" + output_name + ".pdf";
  
  // find proper bounds on y axis, given the height of every histogram
  double y_low = 0;
  double y_high = 0;
  for (int i = 0; i  < h.size(); ++i) {
    if (h[i]->GetMaximum() > y_high)
      y_high = h[i]->GetMaximum();
    if (h[i]->GetMinimum() < y_low)
      y_low = h[i]->GetMaximum();
  }
  
  // give a small boundary on top and bottom
  double diff = y_high - y_low;
  y_high += diff * 0.05;
  
  if (!logy)
    y_low -= diff * 0.05;
  else
    y_low = 0.0 + y_low*0.5;
  
  if (y_max != 9999)
    y_high = y_max;
  if (y_min != -9999)
    y_low = y_min;
  
  // set this for the y range for the first histogram
  if (!logy)
    h[0]->GetYaxis()->SetRangeUser(y_low, y_high);
  
  // and axis labels, and title
  h[0]->SetTitle(canvas_title.c_str());
  h[0]->GetXaxis()->SetTitle(x_axis_label.c_str());
  h[0]->GetXaxis()->SetLabelSize(0.06);
  h[0]->GetXaxis()->SetTitleSize(0.075);
  h[0]->GetXaxis()->SetTitleOffset(0.80);
  h[0]->GetYaxis()->SetTitle(y_axis_label.c_str());
  h[0]->GetYaxis()->SetLabelSize(0.06);
  h[0]->GetYaxis()->SetTitleSize(0.075);
  h[0]->GetYaxis()->SetTitleOffset(0.80);
  h[0]->GetYaxis()->CenterTitle(true);
  
  
  
  // generate a canvas
  TCanvas c;
  c.SetLeftMargin(0.12);
  c.SetBottomMargin(0.15);
  
  // pick a set of colors to use
  int chooseColor[11] = {kBlack, kRed, kBlue, kGreen, kCyan, kMagenta,
                         kOrange, kYellow, kRed+2, kGreen+3, kBlue-7};
  
  // print histograms, giving them some nominal settings to differentiate them
  for (int i = 0; i < h.size(); ++i) {
    h[i]->SetLineColor(chooseColor[i%11]);
    h[i]->SetMarkerColor(chooseColor[i%11]);
    h[i]->SetMarkerSize(1);
    h[i]->SetLineWidth(2);
    h[i]->SetMarkerStyle(21);
    
    if (i == 0)
      h[i]->Draw();
    else
      h[i]->Draw("SAME");
  }
  
  if (do_legend) {
    TLegend* leg = new TLegend(0.65, 0.65, .88, .88);
    leg->SetTextSize(0.04);
    leg->SetHeader(legend_title.c_str());
    for (int i = 0; i < h.size(); ++i) {
      leg->AddEntry(h[i], hist_titles[i].c_str(), "lep");
    }
    leg->Draw();
  }
  if (logx)
    c.SetLogx();
  if (logy)
    c.SetLogy();
  c.SaveAs(canvas_name.c_str());
}

template<typename H>
void Overlay1D(H* h1,
               H* h2,
               std::string h1_title,
               std::string h2_title,
               std::string output_loc,
               std::string output_name,
               std::string canvas_title,
               std::string x_axis_label,
               std::string y_axis_label,
               bool logx,
               bool logy,
               bool do_legend,
               std::string legend_title) {
  // we assume the output location exists, so create
  // the final output string that will be used for pdf creation
  std::string canvas_name = output_loc + "/" + output_name + ".pdf";
  
  // find proper bounds on y axis, given the height of every histogram
  double y_low = 0;
  double y_high = 0;
  if (h1->GetMaximum() > y_high)
    y_high = h1->GetMaximum();
  if (h1->GetMinimum() < y_low)
    y_low = h1->GetMinimum();
  if (h2->GetMaximum() > y_high)
    y_high = h2->GetMaximum();
  if (h2->GetMinimum() < y_low)
    y_low = h2->GetMinimum();
  
  
  // give a small boundary on top and bottom
  double diff = y_high - y_low;
  y_high += diff * 0.05;
  if (!logy)
    y_low -= diff * 0.05;
  else
    y_low = 0.0 + y_low*0.5;
  
  // set this for the y range for the first histogram
  //h1->GetYaxis()->SetRangeUser(y_low, y_high);
  
  // and axis labels, and title
  h1->SetTitle(canvas_title.c_str());
  h1->GetXaxis()->SetTitle(x_axis_label.c_str());
  h1->GetXaxis()->SetTitleSize(0.075);
  h1->GetXaxis()->SetTitleOffset(0.80);
  h1->GetXaxis()->SetLabelSize(0.06);
  h1->GetYaxis()->SetTitle(y_axis_label.c_str());
  h1->GetYaxis()->SetTitleSize(0.075);
  h1->GetYaxis()->SetTitleOffset(0.80);
  h1->GetYaxis()->SetLabelSize(0.06);
  
  h1->SetLineWidth(2);
  h1->SetLineColor(kBlack);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(21);
  h2->SetLineWidth(2);
  h2->SetLineColor(kRed);
  h2->SetMarkerSize(1);
  h2->SetMarkerColor(kRed);
  h2->SetMarkerStyle(22);

  TCanvas c;
  c.SetLeftMargin(0.12);
  c.SetBottomMargin(0.15);
  h1->Draw();
  h2->Draw("SAME");
  
  if (do_legend) {
    TLegend* leg = new TLegend(0.65, 0.65, .88, .88);
    leg->SetHeader(legend_title.c_str());
    leg->AddEntry(h1, h1_title.c_str(), "lep");
    leg->AddEntry(h2, h2_title.c_str(), "lep");
    leg->Draw();
  }
  if (logx)
    c.SetLogx();
  if (logy)
    c.SetLogy();
  c.SaveAs(canvas_name.c_str());
}

template<typename H>
void Print2DSimple(H* h,
                   std::string output_loc,
                   std::string output_name,
                   std::string canvas_title,
                   std::string x_axis_label,
                   std::string y_axis_label,
                   bool logx,
                   bool logy,
                   bool logz,
                   std::string opt = "COLZ") {
  // we assume the output location exists, so create
  // the final output string that will be used for pdf creation
  std::string canvas_name = output_loc + "/" + output_name + ".pdf";
  
  // and axis labels, and title
  h->SetTitle(canvas_title.c_str());
  h->GetXaxis()->SetTitle(x_axis_label.c_str());
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetTitleSize(0.075);
  h->GetXaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->SetTitle(y_axis_label.c_str());
  h->GetYaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetTitleSize(0.075);
  h->GetYaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->CenterTitle(true);
  
  TCanvas c;
  c.SetLeftMargin(0.12);
  c.SetBottomMargin(0.15);
  if (logx)
    c.SetLogx();
  if (logy)
    c.SetLogy();
  if (logz)
    c.SetLogz();
  
  h->Draw(opt.c_str());
  c.SaveAs(canvas_name.c_str());
}

#endif // ROOT_PRINT_ROUTINES_HH

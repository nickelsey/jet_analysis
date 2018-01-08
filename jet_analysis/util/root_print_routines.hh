#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>
#include <vector>
#include <string>

#include "boost/filesystem.hpp"

void Overlay1D(const std::vector<TH1D*>& h,
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
  h[0]->GetYaxis()->SetTitle(y_axis_label.c_str());
  
  
  // generate a canvas
  TCanvas c;
  
  // pick a set of colors to use
  int chooseColor[11] = {kBlack, kRed, kBlue, kGreen, kCyan, kYellow,
                        kMagenta, kOrange, kRed+2, kGreen+3,
                        kBlue-7};
  
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
    TLegend* leg = new TLegend(0.6, 0.6, .88, .88);
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

void Overlay1D(TH1D* h1,
               TH1D* h2,
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
  h1->GetYaxis()->SetTitle(y_axis_label.c_str());
  
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
  h1->Draw();
  h2->Draw("SAME");
  
  if (do_legend) {
    TLegend* leg = new TLegend(0.6, 0.6, .88, .88);
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

void Print2DSimple(TH2D* h,
                   std::string output_loc,
                   std::string output_name,
                   std::string canvas_title,
                   std::string x_axis_label,
                   std::string y_axis_label,
                   bool logx,
                   bool logy,
                   bool logz) {
  // we assume the output location exists, so create
  // the final output string that will be used for pdf creation
  std::string canvas_name = output_loc + "/" + output_name + ".pdf";
  
  h->SetTitle(canvas_title.c_str());
  h->GetXaxis()->SetTitle(x_axis_label.c_str());
  h->GetYaxis()->SetTitle(y_axis_label.c_str());
  
  TCanvas c;
  if (logx)
    c.SetLogx();
  if (logy)
    c.SetLogy();
  if (logz)
    c.SetLogz();
  
  h->Draw("COLZ");
  c.SaveAs(canvas_name.c_str());
}

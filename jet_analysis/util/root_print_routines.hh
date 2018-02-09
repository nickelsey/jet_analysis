
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
                   int color = 1,
                   double y_min = -9999,
                   double y_max = 9999) {
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
  h->GetYaxis()->SetLabelSize(0.055);
  h->GetYaxis()->SetTitleSize(0.075);
  h->GetYaxis()->SetTitleOffset(0.78);
  h->GetYaxis()->CenterTitle(true);
  
  // generate a canvas
  TCanvas c;
  c.SetLeftMargin(0.12);
  c.SetBottomMargin(0.15);
  
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerSize(1);
  h->SetLineWidth(2);
  h->SetMarkerStyle(21);
  if (y_max != 9999 && y_min != -9999) {
    h->GetYaxis()->SetRangeUser(y_min, y_max);
  }
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
  h[0]->GetYaxis()->SetLabelSize(0.055);
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
  h1->GetYaxis()->SetLabelSize(0.055);
  
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
  h->GetYaxis()->SetLabelSize(0.055);
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

// print three histograms and their ratio
template<typename H>
void printWithRatio3(H* h1,
                     H* h2,
                     H* h3,
                     std::string h1_title,
                     std::string h2_title,
                     std::string h3_title,
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
  
  
  // force error calculation
  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();
  
  // Define the Canvas
  TCanvas *c = new TCanvas("c", "canvas", 800, 800);
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetLeftMargin(0.12);
  //pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  h1->SetStats(0);          // No statistics on upper plot
  h1->SetTitle("");
  if (logx) pad1->SetLogx();
  if (logy) pad1->SetLogy();
  
  // set styles for h1, h2, h3
  h1->SetTitle(canvas_title.c_str());
  h1->GetXaxis()->SetTitle("");
  h1->GetYaxis()->SetTitle(y_axis_label.c_str());
  h1->GetYaxis()->SetTitleSize(0.075);
  h1->GetYaxis()->SetTitleOffset(0.80);
  h1->GetYaxis()->SetLabelSize(0.055);
  
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
  h2->SetLineWidth(2);
  h2->SetLineColor(kBlue);
  h2->SetMarkerSize(1);
  h2->SetMarkerColor(kBlue);
  h2->SetMarkerStyle(22);
  
  h1->Draw();               // Draw h1
  h2->Draw("same");         // Draw h2 on top of h1
  h3->Draw("same");         // Draw h3 on top of h1 & h2
  
  if (do_legend) {
  TLegend* leg = new TLegend( 0.6, 0.7, 0.9, 0.9 );
    leg->AddEntry( h1, h1_title.c_str(), "lep" );
    leg->AddEntry( h2, h2_title.c_str(), "lep" );
    leg->AddEntry( h3, h3_title.c_str(), "lep" );
    leg->Draw();
  }
  // Do not draw the Y axis label on the upper plot and redraw a small
  // axis instead, in order to avoid the first label (0) to be clipped.
  //h1->GetYaxis()->SetLabelOffset( 10 );
  //h1->GetYaxis()->SetLabelSize(0.);
  //TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
  //axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  //axis->SetLabelSize(15);
  //axis->Draw();
  //
  // lower plot will be in pad
  c->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.4);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetLeftMargin(0.12);
  if ( logx ) pad2->SetLogx();
  //pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  
  // Define the ratio plot h1/h2
  TH1D *hR = (TH1D*)h1->Clone("hR");
  hR->SetLineColor(h1->GetLineColor());
  hR->SetMarkerColor(h1->GetLineColor());
  hR->SetMarkerStyle(h2->GetMarkerStyle());
  hR->SetMinimum(0.8);  // Define Y ..
  hR->SetMaximum(1.35); // .. range
  hR->Sumw2();
  hR->SetStats(0);      // No statistics on lower plot
  hR->Divide(h3);
  hR->Draw("ep");       // Draw the ratio plot
  
  // Define the ratio plot h1/h3
  TH1D *hR1 = (TH1D*)h1->Clone("hR1");
  hR1->SetLineColor(h1->GetLineColor());
  hR1->SetMarkerColor(h1->GetLineColor());
  hR1->SetMarkerStyle(h1->GetMarkerStyle());
  hR1->SetMinimum(0.8);  // Define Y ..
  hR1->SetMaximum(1.35); // .. range
  hR1->Sumw2();
  hR1->SetStats(0);      // No statistics on lower plot
  hR1->Divide(h3);
  hR1->Draw("same");       // Draw the ratio plot
  
  // Define the ratio plot h2/h3
  TH1D *hR2 = (TH1D*)h2->Clone("hR2");
  hR2->SetLineColor(h2->GetLineColor());
  hR2->SetMarkerColor(h2->GetLineColor());
  hR2->SetMarkerStyle(h2->GetMarkerStyle());
  hR2->SetMinimum(0.8);  // Define Y ..
  hR2->SetMaximum(1.35); // .. range
  hR2->Sumw2();
  hR2->SetStats(0);      // No statistics on lower plot
  hR2->Divide(h3);
  hR2->Draw("same");       // Draw the ratio plot
  
  // Y axis h1 plot settings
//  h1->GetYaxis()->SetTitleSize(20);
//  h1->GetYaxis()->SetTitleFont(43);
//  h1->GetYaxis()->SetTitleOffset(1.5);
  //h1->GetYaxis()->SetTitleOffset(1.55);
  
  // Ratio plot (h3) settings
  hR->SetTitle(""); // Remove the ratio title
  
  // Y axis ratio plot settings
  hR->GetYaxis()->SetTitle("ratio ");
  hR->GetXaxis()->SetTitle(x_axis_label.c_str() );
  hR->GetYaxis()->SetNdivisions(505);
  hR->GetYaxis()->SetTitleSize(20);
  hR->GetYaxis()->SetTitleFont(43);
  hR->GetYaxis()->SetTitleOffset(1.55);
  //hR->GetYaxis()->SetLabelFont(44); // Absolute font size in pixel (precision 3)
  //hR->GetYaxis()->SetLabelSize(15);
  hR->GetYaxis()->SetTitleSize(0.1);
  //hR->GetYaxis()->SetTitleOffset(0.80);
  hR->GetYaxis()->SetLabelSize(0.1);
  
  // X axis ratio plot settings
  hR->GetXaxis()->SetTitleSize(20);
  hR->GetXaxis()->SetTitleFont(43);
  //hR->GetXaxis()->SetTitleOffset(3.);
  //hR->GetXaxis()->SetLabelFont(44); // Absolute font size in pixel (precision 3)
  //hR->GetXaxis()->SetLabelSize(15);
  hR->GetXaxis()->SetTitleSize(0.075);
  //hR->GetXaxis()->SetTitleOffset(1);
  hR->GetXaxis()->SetLabelSize(0.1);
  
  
  c->SaveAs(canvas_name.c_str());
  delete c;
  
}

#endif // ROOT_PRINT_ROUTINES_HH

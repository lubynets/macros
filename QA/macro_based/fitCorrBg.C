class classFuncBkgWithTemplate : public TNamed {

    public:
    classFuncBkgWithTemplate(TH1D* histo_template_input, int polDegreeInput, bool addSignalInput) : TNamed(), addSignal(false), histo_template(nullptr), polDegree(2) { // constructor
        
        /// get the template
        if(histo_template_input) {
            std::cout << "[classFuncBkgWithTemplate] template read correctly" << std::endl;
            histo_template = dynamic_cast<TH1D*>(histo_template_input->Clone());
            histo_template->Smooth(100000);
        }

        /// polynomial background
        polDegree = polDegreeInput;

        /// add signal
        addSignal = addSignalInput;

        std::cout << "[classFuncBkgWithTemplate] Constructor:" << std::endl;
        std::cout << "[classFuncBkgWithTemplate]    - addSignal: " << addSignal << std::endl;
        std::cout << "[classFuncBkgWithTemplate]    - histo_template: " << histo_template << std::endl;
        std::cout << "[classFuncBkgWithTemplate]    - polDegree: " << polDegree << std::endl;
    }

    double operator()(double *x, double *par) {

        float return_value = 0;
        int npars = 0;
        // the template
        if(histo_template) {
            return_value += par[npars]*histo_template->GetBinContent(histo_template->FindBin(x[0]));
            npars++;
        }

        /// add signal
        if (addSignal) {
            return_value += (par[npars+0]/(TMath::Sqrt(2*TMath::Pi())*par[npars+2])) * TMath::Exp( -0.5*(x[0]-par[npars+1])*(x[0]-par[npars+1])/(par[npars+2]*par[npars+2]) );
            npars += 3;
        }

        /// polynomial background
        switch (polDegree)
          {
          case 0:
              return_value += par[npars];
              npars += 1;
              break;
          case 1:
              return_value += par[npars] + par[npars+1]*x[0];
              npars += 2;
              break;
          case 2:
              return_value += par[npars] + par[npars+1]*x[0] + par[npars+2]*x[0]*x[0];
              npars += 3;
              break;
          case 3:
              return_value += par[npars] + par[npars+1]*x[0] + par[npars+2]*x[0]*x[0] + par[npars+3]*x[0]*x[0]*x[0];
              npars += 4;
              break;
          case 4:
              return_value += par[npars] + par[npars+1]*x[0] + par[npars+2]*x[0]*x[0] + par[npars+3]*x[0]*x[0]*x[0] + par[npars+4]*x[0]*x[0]*x[0]*x[0];
              npars += 4;
              break;
          
          default:
          std::cout << "[classFuncBkgWithTemplatel] no correct polynomial degree defined: " << polDegree << " (must be within 0 and 4)" << std::endl;
              break;
          }
        
        return return_value;
    }

    private:
        bool addSignal;
        TH1D* histo_template;
        int polDegree;
};
///____________________________________________________
///____________________________________________________
///____________________________________________________
void fitCorrBg() {
    const std::string filename = "/home/oleksii/alidir/working/cutVar/data/input/mass_bdt_qa_thn.HL.data.HF_LHC23_pass4_Thin_2P3PDstar.522578.root";
    const std::string histoname = "pT_1_20/T_0.90_1.60/hM_NPgt0.00";
    const std::string filename_template = "/home/oleksii/alidir/working/correlBG/corrBkgLc.canvaser.HL.mc.HF_LHC24h1b_All.522675.mainOnly.root";
    const std::string histoname_template = "TBin4/histos_bkgSum_TBin4";

    gStyle->SetCanvasPreferGL(true);

    const float massLc = 2.286; // GeV/c2
    const float minX = 2.1;
    const float maxX = 2.45;

    /// retrieve the histogram
    TFile* file = TFile::Open(filename.c_str());
    TH1D* histogram = nullptr;
    file->GetObject(histoname.c_str(), histogram);
    TH1D* histogram_2 = dynamic_cast<TH1D*>(histogram->Clone());
    
    //////////////////////////////////
    ///                            ///
    ///   1st fit - gauss + pol2   ///
    ///                            ///
    //////////////////////////////////
    
    /// canvas
    TCanvas* can1 = new TCanvas("can1", "no template", 1200, 600);
    can1->Divide(2, 1);
    float minY = histogram->GetMinimum();
    float maxY = histogram->GetMaximum();
    TH1F* hFrame_can1 = can1->cd(1)->DrawFrame(minX, minY*0.2, maxX, maxY*1.1, histogram->GetTitle());
    hFrame_can1->GetXaxis()->SetTitle(histogram->GetXaxis()->GetTitle()); 
    hFrame_can1->GetYaxis()->SetTitle(histogram->GetYaxis()->GetTitle());
    histogram->Draw("same");
    hFrame_can1->GetXaxis()->UnZoom();
    gPad->SetTicks();

    /// fit function
    std::string str_signal = "([0]/(TMath::Sqrt(2*TMath::Pi())*[2])) * TMath::Exp( -0.5*(x[0]-[1])*(x[0]-[1])/([2]*[2]) )";
    std::string str_bkg = "pol3";
    TF1* fit_noTemplate = new TF1("fit_noTemplate", Form("%s + %s(3)", str_signal.c_str(), str_bkg.c_str()), minX, maxX);
    TF1* func_signal_1 = new TF1("func_signal_1", str_signal.c_str(), minX, maxX);
    TF1* func_bkg_1 = new TF1("func_bkg_1", str_bkg.c_str(), minX, maxX);
    
    /// pre-fit background
    TF1* prefit_pol3 = new TF1("prefit_pol3", "pol3", 0, 3);
    histogram->Fit(prefit_pol3, "NQR", "", minX, maxX);

    /// set parameters
    fit_noTemplate->SetParameters(/*histogram->GetBinContent(histogram->FindBin(massLc))*/ 1000, massLc, 0.01, prefit_pol3->GetParameter(0), prefit_pol3->GetParameter(1), prefit_pol3->GetParameter(2), prefit_pol3->GetParameter(3));
    
    /// fit
    histogram->Fit(fit_noTemplate, "RN", "", minX, maxX);
    fit_noTemplate->SetLineColor(kBlue+1);
    func_bkg_1->Draw("same");
    fit_noTemplate->Draw("same");
    const float binwidth = histogram->GetBinWidth(1);
    const float signal_1 = fit_noTemplate->GetParameter(0) / binwidth;
    const float signal_error_1 = fit_noTemplate->GetParError(0) / binwidth;
    const float mean_1 = fit_noTemplate->GetParameter(1);
    const float mean_error_1 = fit_noTemplate->GetParError(1);
    const float sigma_1 = fit_noTemplate->GetParameter(2);
    const float sigma_error_1 = fit_noTemplate->GetParError(2);

    TPaveText* txt_1 = new TPaveText(minX, maxY*0.20, maxX, maxY*0.40);
    txt_1->SetTextFont(42);
    txt_1->SetTextSize(0.045);
    txt_1->SetFillStyle(0);
    txt_1->SetBorderSize(0);
    txt_1->AddText(Form("#mu: %.6f #pm %.6f GeV/#it{c}^{2}", mean_1, mean_error_1));
    txt_1->AddText(Form("#sigma: %.6f #pm %.6f GeV/#it{c}^{2}", sigma_1, sigma_error_1));
    txt_1->AddText(Form("Signal: %.0f #pm %.0f", signal_1, signal_error_1));
    txt_1->Draw();

    /// signal and background functions
    func_signal_1->SetParameters(fit_noTemplate->GetParameter(0), fit_noTemplate->GetParameter(1), fit_noTemplate->GetParameter(2));
    func_bkg_1->SetParameters(fit_noTemplate->GetParameter(3), fit_noTemplate->GetParameter(4), fit_noTemplate->GetParameter(5), fit_noTemplate->GetParameter(6));

    /// residuals
    TH1D* histo_residuals_noTempl = dynamic_cast<TH1D*>(histogram->Clone());
    auto list = histo_residuals_noTempl->GetListOfFunctions();
    std::vector<float> values = {};
    for(int bin=1; bin<=histo_residuals_noTempl->GetNbinsX(); bin++) {
        const float center = histo_residuals_noTempl->GetBinCenter(bin);
        if(center < minX || center > maxX) {
            continue;
        }
        histo_residuals_noTempl->SetBinContent(bin, histogram->GetBinContent(bin) - func_bkg_1->Eval(center));
        values.push_back(histo_residuals_noTempl->GetBinContent(bin));
    }
    float minY_residuals = *std::min_element(values.begin(), values.end());
    float maxY_residuals = *std::max_element(values.begin(), values.end());
    TH1F* hFrame_residuals_can1 = can1->cd(2)->DrawFrame(2.1, minY_residuals*1.1, 2.45, maxY_residuals*1.1, histo_residuals_noTempl->GetTitle());
    hFrame_residuals_can1->GetXaxis()->SetTitle(histogram->GetXaxis()->GetTitle()); 
    hFrame_residuals_can1->GetYaxis()->SetTitle("residuals");
    histo_residuals_noTempl->Draw("same");
    func_signal_1->SetLineWidth(2);
    func_signal_1->SetLineColor(kRed+1);
    func_signal_1->SetFillColor(func_signal_1->GetLineColor());
    //TH1D* h_signal = dynamic_cast<TH1D*>(func_signal_1->GetHistogram());
    //h_signal->SetFillColorAlpha(h_signal->GetLineColor(), 0.2);
    //h_signal->SetFillStyle(1001);
    //h_signal->DrawClone("samehisto");
    func_signal_1->Draw("same");
    hFrame_residuals_can1->GetXaxis()->UnZoom();
    gPad->SetTicks();


    ////////////////////////////////////////////
    ///                                      ///
    ///   2nd fit - gauss + pol2 + template  ///
    ///                                      ///
    ////////////////////////////////////////////

    /// canvas
    TCanvas* can2 = new TCanvas("can2", "with template", 1200, 600);
    can2->Divide(2, 1);
    TH1F* hFrame_can2 = can2->cd(1)->DrawFrame(minX, minY*0.2, maxX, maxY*1.1, histogram_2->GetTitle());
    hFrame_can2->GetXaxis()->SetTitle(histogram_2->GetXaxis()->GetTitle());
    hFrame_can2->GetYaxis()->SetTitle(histogram_2->GetYaxis()->GetTitle());
    histogram_2->Draw("same");
    hFrame_can2->GetXaxis()->UnZoom();
    gPad->SetTicks();

    /// template histogram
    TFile* file_template = TFile::Open(filename_template.c_str(), "read");
    TH1D* histo_template = nullptr;
    file_template->GetObject(histoname_template.c_str(), histo_template);

    /// fit function
    classFuncBkgWithTemplate* fFit_bkg_template = new classFuncBkgWithTemplate(histo_template, 3, false);
    classFuncBkgWithTemplate* fFit_template = new classFuncBkgWithTemplate(histo_template, 3, true);
    int npars = 1 /*template*/ + 3 /*signal*/ + 4 /*pol3*/;
    //int npars = 1 /*template*/ + 3 /*signal*/ + 3 /*pol2*/;
    std::cout << "npars = " << npars << std::endl;
    TF1* fit_bkg_template = new TF1("fit_bkg_template", fFit_bkg_template, minX, maxX, npars-3);
    TF1* fit_template = new TF1("fit_template", fFit_template, minX, maxX, npars);

    /// prefit with background only
    fit_bkg_template->SetParameters(5000, func_bkg_1->GetParameter(0), func_bkg_1->GetParameter(1), func_bkg_1->GetParameter(2), func_bkg_1->GetParameter(3));
    histogram_2->Fit(fit_bkg_template, "R", "", minX, maxX);

    /// set parameters
    fit_template->SetParameters(fit_bkg_template->GetParameter(0), func_signal_1->GetParameter(0), func_signal_1->GetParameter(1), func_signal_1->GetParameter(2), fit_bkg_template->GetParameter(1), fit_bkg_template->GetParameter(2), fit_bkg_template->GetParameter(3), fit_bkg_template->GetParameter(4));
    //fit_template->SetParLimits(0, 0., 1000000000.);
    //fit_template->SetParLimits(1, 0., 1000000000.);

    /// fit
    histogram_2->Fit(fit_template, "RL", "", minX, maxX);
    const float signal_2 = fit_template->GetParameter(1) / binwidth;
    const float signal_error_2 = fit_template->GetParError(1) / binwidth;
    const float mean_2 = fit_template->GetParameter(2);
    const float mean_error_2 = fit_template->GetParError(2);
    const float sigma_2 = fit_template->GetParameter(3);
    const float sigma_error_2 = fit_template->GetParError(3);

    TPaveText* txt_2 = new TPaveText(minX, maxY*0.20, maxX, maxY*0.40);
    txt_2->SetTextFont(42);
    txt_2->SetTextSize(0.045);
    txt_2->SetFillStyle(0);
    txt_2->SetBorderSize(0);
    txt_2->AddText(Form("#mu: %.6f #pm %.6f GeV/#it{c}^{2}", mean_2, mean_error_2));
    txt_2->AddText(Form("#sigma: %.6f #pm %.6f GeV/#it{c}^{2}", sigma_2, sigma_error_2));
    txt_2->AddText(Form("Signal: %.0f #pm %.0f", signal_2, signal_error_2));
    txt_2->Draw();

    /// signal and background functions
    TF1* func_signal_2 = new TF1("func_signal_2", str_signal.c_str(), minX, maxX);
    func_signal_2->SetParameters(fit_template->GetParameter(1), fit_template->GetParameter(2), fit_template->GetParameter(3));
    fit_bkg_template->SetParameters(fit_template->GetParameter(0), fit_template->GetParameter(4), fit_template->GetParameter(5), fit_template->GetParameter(6), fit_template->GetParameter(7));

    /// residuals
    TH1D* histo_residuals_templ = dynamic_cast<TH1D*>(histogram_2->Clone());
    auto list_2 = histo_residuals_templ->GetListOfFunctions();
    std::vector<float> values_2 = {};
    for(int bin=1; bin<=histo_residuals_templ->GetNbinsX(); bin++) {
        const float center = histo_residuals_templ->GetBinCenter(bin);
        if(center < minX || center > maxX) {
            continue;
        }
        histo_residuals_templ->SetBinContent(bin, histogram_2->GetBinContent(bin) - fit_bkg_template->Eval(center));
        values_2.push_back(histo_residuals_templ->GetBinContent(bin));
    }
    float minY_residuals_2 = *std::min_element(values_2.begin(), values_2.end());
    float maxY_residuals_2 = *std::max_element(values_2.begin(), values_2.end());
    TH1F* hFrame_residuals_can2 = can2->cd(2)->DrawFrame(2.1, minY_residuals_2*1.1, 2.45, maxY_residuals_2*1.1, histo_residuals_noTempl->GetTitle());
    hFrame_residuals_can2->GetXaxis()->SetTitle(histogram_2->GetXaxis()->GetTitle()); 
    hFrame_residuals_can2->GetYaxis()->SetTitle("residuals");
    histo_residuals_templ->Draw("same");
    func_signal_2->SetLineWidth(2);
    func_signal_2->SetLineColor(kRed+1);
    func_signal_2->SetFillColor(func_signal_2->GetLineColor());
    //TH1D* h_signal = dynamic_cast<TH1D*>(func_signal_2->GetHistogram());
    //h_signal->SetFillColorAlpha(h_signal->GetLineColor(), 0.2);
    //h_signal->SetFillStyle(1001);
    //h_signal->DrawClone("samehisto");
    func_signal_2->Draw("same");
    hFrame_residuals_can2->GetXaxis()->UnZoom();
    gPad->SetTicks();

    std::cout << "===>  signal w/o template / signal with template = " << signal_1 << " / " << signal_2 << " = " << signal_1/signal_2 << std::endl;

    can1->SaveAs("wo_template.pdf");
    can2->SaveAs("with_template.pdf");

    return;
}

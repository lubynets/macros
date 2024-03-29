TH2F* Get2DLayerFromTH3(TH3F* histo3D, TString axis, int bin);
void SetAxesNames(TH2F* histo, TString xaxisname, TString yaxisname);
std::string StringBinNumber(int number);

void TH3ToTH2Complex()
{
  gStyle -> SetOptStat(0);
  
  TFile* file_fit = TFile::Open("/home/user/cbmdir/working/qna/fits/out.fitter.apr20.dcmqgsm.nopid.lightcuts1.set4.XX.root");
  TFile* file_mcfit = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcfitter.apr20.dcmqgsm.nopid.lightcuts1.set4.XX.root");
  TFile* file_mcv1 = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcv1.apr20.dcmqgsm.nopid.lightcuts1.set4.XX.root");
  
  bool is_ave = false;
  
  struct axis
  {
    std::string name_;
    std::string title_;
    std::string id_;
    int another_first_;
    int another_second_;
    int nbins_;
  };
  
  std::vector<axis> axes
  {
    {"centrality", "Centrality, %", "x", 1, 2, -1},
    {"rapidity",   "y_{LAB}",       "y", 0, 2, -1},
    {"pT",         "p_{T}, GeV/c",  "z", 0, 1, -1}
  };
  
  struct infotype
  {
    std::string name_;
    std::string folder_;
    bool is_mcfitter_;
    bool is_mcv1_;
    bool not4ave_;
  };
  
  std::vector<infotype> infotypes
  {
    {"hsignal",        "sgnl",            true,  true , false},
    {"hbckgr_0",       "bckgr/intercept", true,  false, false},
    {"hbckgr_1",       "bckgr/slope",     true,  false, false},
    {"hfit_chi2ndf",   "chi2ndf",         false, false, true },
    {"hentries_sgnl",  "Nentries/sgnl",   false, false, true },
    {"hentries_bckgr", "Nentries/bckgr",  false, false, true }
  };
  
  TH3F* h_fit;
  TH3F* h_mcfit;
  TH3F* h_mcv1;
  
  TH2F* h2_fit;
  TH2F* h2_mcfit;
  TH2F* h2_mcv1;
      
  TFile* fileOut = TFile::Open("out.th3toth2cplx.root", "recreate");

  bool is_first_canvas = true;  
  
  for(auto it : infotypes)
  {
    if(it.not4ave_&&is_ave == true) continue;
  
                        h_fit   = file_fit   -> Get<TH3F>(("parameters/" + it.name_).c_str());
    if(it.is_mcfitter_) h_mcfit = file_mcfit -> Get<TH3F>(("parameters/" + it.name_).c_str());
    if(it.is_mcv1_)     h_mcv1  = file_mcv1  -> Get<TH3F>("hv1_mc");
    
    axes.at(0).nbins_ = h_fit->GetXaxis()->GetNbins();
    axes.at(1).nbins_ = h_fit->GetYaxis()->GetNbins();
    axes.at(2).nbins_ = h_fit->GetZaxis()->GetNbins();
    
    for(auto ax : axes)
    {     
      for(int ibin=1; ibin<=ax.nbins_; ibin++)
      {
        std::string layername = StringBinNumber(ibin);
        std::string dirname = it.folder_ +"/" + ax.name_ + "/" + layername;
        fileOut -> mkdir(dirname.c_str());
        fileOut -> cd(dirname.c_str()); 
        
                            h2_fit   = Get2DLayerFromTH3(h_fit  , ax.id_, ibin);
        if(it.is_mcfitter_) h2_mcfit = Get2DLayerFromTH3(h_mcfit, ax.id_, ibin);
        if(it.is_mcv1_)     h2_mcv1  = Get2DLayerFromTH3(h_mcv1 , ax.id_, ibin);
        
                            SetAxesNames(h2_fit  , axes.at(ax.another_first_).title_, axes.at(ax.another_second_).title_);
        if(it.is_mcfitter_) SetAxesNames(h2_mcfit, axes.at(ax.another_first_).title_, axes.at(ax.another_second_).title_);
        if(it.is_mcv1_)     SetAxesNames(h2_mcv1 , axes.at(ax.another_first_).title_, axes.at(ax.another_second_).title_);  
        
                            h2_fit   -> Write("fit");
        if(it.is_mcfitter_) h2_mcfit -> Write("mcfit");
        if(it.is_mcv1_)     h2_mcv1  -> Write("mcv1");
        
        if(it.name_ != "hsignal") continue;
        
        TCanvas cc_fit("canvas_fit", "canvas_fit", 1500, 900);
        cc_fit.cd();
        cc_fit.SetLeftMargin(0.08);
        cc_fit.SetRightMargin(0.12);
        h2_fit -> GetZaxis() -> SetTitle("v_{1x}");
        h2_fit -> SetTitle(("fit, " + std::string(h2_fit->GetTitle())).c_str());
        h2_fit -> Draw("colz+text");
        if(is_first_canvas)
          cc_fit.Print("out.th2.cplx.pdf(", "pdf");
        else
          cc_fit.Print("out.th2.cplx.pdf", "pdf");
        
        if(it.is_mcfitter_) { TCanvas cc_mcfit("canvas_mcfit", "canvas_mcfit", 1500, 900);
                              cc_mcfit.cd();
                              cc_mcfit.SetLeftMargin(0.08);
                              cc_mcfit.SetRightMargin(0.12);
                              h2_mcfit -> GetZaxis() -> SetTitle("v_{1x}");
                              h2_mcfit -> SetTitle(("mcfit, " + std::string(h2_mcfit->GetTitle())).c_str());
                              h2_mcfit -> Draw("colz+text");
                              cc_mcfit.Print("out.th2.cplx.pdf", "pdf");
        }
                              
        if(it.is_mcv1_)     { TCanvas cc_mcv1("canvas_mcv1", "canvas_mcv1", 1500, 900);
                              cc_mcv1.cd();
                              cc_mcv1.SetLeftMargin(0.08);
                              cc_mcv1.SetRightMargin(0.12);
                              h2_mcv1 -> GetZaxis() -> SetTitle("v_{1x}");
                              h2_mcv1 -> SetTitle(("mcv1, " + std::string(h2_mcv1->GetTitle())).c_str());
                              h2_mcv1 -> Draw("colz+text");
                              cc_mcv1.Print("out.th2.cplx.pdf", "pdf");
        }
                                      
        is_first_canvas = false;
      }      
    }    
  } 
  
  TCanvas emptycanvas("", "", 1500, 900);
  emptycanvas.Print("out.th2.cplx.pdf)", "pdf");
  
  fileOut -> Close();
}

TH2F* Get2DLayerFromTH3(TH3F* histo3D, TString axis, int bin)
{
  if(axis=="X" || axis =="x" || axis=="YZ" || axis=="yz")
  {
    const int nbins_1 = histo3D->GetYaxis()->GetNbins();
    const int nbins_2 = histo3D->GetZaxis()->GetNbins();
    
    double* bin_edges_1 = new double[nbins_1+1];
    double* bin_edges_2 = new double[nbins_2+1];
    
    for(int i=0; i<=nbins_1; i++)
      bin_edges_1[i] = histo3D->GetYaxis()->GetBinLowEdge(i+1);
    for(int i=0; i<=nbins_2; i++)
      bin_edges_2[i] = histo3D->GetZaxis()->GetBinLowEdge(i+1);
    
    TH2F* histo2D = new TH2F("", "", nbins_1, bin_edges_1, nbins_2, bin_edges_2);
    
    for(int i=1; i<=nbins_1; i++)
      for(int j=1; j<=nbins_2; j++)
      {
        histo2D -> SetBinContent(i, j, histo3D->GetBinContent(bin, i, j));
        histo2D -> SetBinError(i, j, histo3D->GetBinError(bin, i, j));
      }
      
    histo2D -> SetTitle((std::string(histo3D->GetXaxis()->GetTitle()) + " from " + std::to_string(histo3D->GetXaxis()->GetBinLowEdge(bin)) + " to " + std::to_string(histo3D->GetXaxis()->GetBinLowEdge(bin+1))).c_str());
    histo2D -> GetXaxis() -> SetTitle(histo3D->GetYaxis()->GetTitle());
    histo2D -> GetYaxis() -> SetTitle(histo3D->GetZaxis()->GetTitle());
    
    return histo2D;
  }
  
  if(axis=="Y" || axis =="y" || axis=="XZ" || axis=="xz")
  {
    const int nbins_1 = histo3D->GetXaxis()->GetNbins();
    const int nbins_2 = histo3D->GetZaxis()->GetNbins();
    
    double* bin_edges_1 = new double[nbins_1+1];
    double* bin_edges_2 = new double[nbins_2+1];
    
    for(int i=0; i<=nbins_1; i++)
      bin_edges_1[i] = histo3D->GetXaxis()->GetBinLowEdge(i+1);
    for(int i=0; i<=nbins_2; i++)
      bin_edges_2[i] = histo3D->GetZaxis()->GetBinLowEdge(i+1);
    
    TH2F* histo2D = new TH2F("", "", nbins_1, bin_edges_1, nbins_2, bin_edges_2);
    
    for(int i=1; i<=nbins_1; i++)
      for(int j=1; j<=nbins_2; j++)
      {
        histo2D -> SetBinContent(i, j, histo3D->GetBinContent(i, bin, j));
        histo2D -> SetBinError(i, j, histo3D->GetBinError(i, bin, j));
      }
      
    histo2D -> SetTitle((std::string(histo3D->GetYaxis()->GetTitle()) + " from " + std::to_string(histo3D->GetYaxis()->GetBinLowEdge(bin)) + " to " + std::to_string(histo3D->GetYaxis()->GetBinLowEdge(bin+1))).c_str());
    histo2D -> GetXaxis() -> SetTitle(histo3D->GetXaxis()->GetTitle());
    histo2D -> GetYaxis() -> SetTitle(histo3D->GetZaxis()->GetTitle());

    return histo2D;
  }
  
  if(axis=="Z" || axis =="z" || axis=="XY" || axis=="xy")
  {
    const int nbins_1 = histo3D->GetXaxis()->GetNbins();
    const int nbins_2 = histo3D->GetYaxis()->GetNbins();
    
    double* bin_edges_1 = new double[nbins_1+1];
    double* bin_edges_2 = new double[nbins_2+1];
    
    for(int i=0; i<=nbins_1; i++)
      bin_edges_1[i] = histo3D->GetXaxis()->GetBinLowEdge(i+1);
    for(int i=0; i<=nbins_2; i++)
      bin_edges_2[i] = histo3D->GetYaxis()->GetBinLowEdge(i+1);
    
    TH2F* histo2D = new TH2F("", "", nbins_1, bin_edges_1, nbins_2, bin_edges_2);
    
    for(int i=1; i<=nbins_1; i++)
      for(int j=1; j<=nbins_2; j++)
      {
        histo2D -> SetBinContent(i, j, histo3D->GetBinContent(i, j, bin));
        histo2D -> SetBinError(i, j, histo3D->GetBinError(i, j, bin));
      }
      
    histo2D -> SetTitle((std::string(histo3D->GetZaxis()->GetTitle()) + " from " + std::to_string(histo3D->GetZaxis()->GetBinLowEdge(bin)) + " to " + std::to_string(histo3D->GetZaxis()->GetBinLowEdge(bin+1))).c_str());
    histo2D -> GetXaxis() -> SetTitle(histo3D->GetXaxis()->GetTitle());
    histo2D -> GetYaxis() -> SetTitle(histo3D->GetYaxis()->GetTitle());

    return histo2D;
  }
  
}

void SetAxesNames(TH2F* histo, TString xaxisname, TString yaxisname)
{
  histo -> GetXaxis() -> SetTitle(xaxisname);
  histo -> GetYaxis() -> SetTitle(yaxisname);
}

std::string StringBinNumber(int number)
{
  if(number<10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}
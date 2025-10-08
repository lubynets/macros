# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import uproot
import hist
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
import ROOT
from PyPDF2 import PdfWriter

# Define pT bins (bin edges for the X-axis of the histogram)
pt_bins = [1, 2, 3, 4, 5, 6, 8, 12]  # Example pT bin edges
n_bins = len(pt_bins) - 1  # Number of pT bins
fix_sigma = True
fix_file = "output/rawyield/flare/RawYields_Lc_run3_pp_13.6TeV_default.root"

# Iterate over the input files from 0 to 20
for file_index in range(0, 1):  # Looping over 0 to 20
    # Define the input file name for each iteration
    input_file = f"/Users/zbiao/analysis/DmesonAnalysis/run3LcMCCheck/output/rawyield/Dis_cutset_Lc_run3_pp_13.6TeV_for_femto_dec18.root"
    
    # Prepare histograms in ROOT for fit results (for each file_index)
    mass_hist = ROOT.TH1F(f"hRawYieldsMean", f"Mass vs pT; pT (GeV/c); Mass (GeV/c^2)", n_bins, np.array(pt_bins, dtype="float64"))
    sigma_hist = ROOT.TH1F(f"hRawYieldsSigma", f"Sigma vs pT; pT (GeV/c); Sigma (GeV/c^2)", n_bins, np.array(pt_bins, dtype="float64"))
    raw_yield_hist = ROOT.TH1F(f"hRawYields", f"Raw Yield vs pT; pT (GeV/c); Raw Yield", n_bins, np.array(pt_bins, dtype="float64"))
    bkg_hist = ROOT.TH1F(f"hRawYieldsBkg", f"Background vs pT; pT (GeV/c); Background", n_bins, np.array(pt_bins, dtype="float64"))
    significance_hist = ROOT.TH1F(f"hRawYieldsSignificance", f"Significance vs pT; pT (GeV/c); Significance", n_bins, np.array(pt_bins, dtype="float64"))
    chi2_hist = ROOT.TH1F(f"hRawYieldsChiSquare", f"Chi2 vs pT; pT (GeV/c); Bkg", n_bins, np.array(pt_bins, dtype="float64"))
    soverb_hist = ROOT.TH1F(f"hRawYieldsSoverB", f"Signal over Background vs pT; pT (GeV/c); S / B", n_bins, np.array(pt_bins, dtype="float64"))
    purity_hist = ROOT.TH1F("hRawYieldsPurity", "Purity vs kstar; k* (GeV/c); Purity", n_bins, np.array(pt_bins, dtype="float64"))

  
    pdf_files = []

    # Loop over each pT bin
    for idx, (pt_min, pt_max) in enumerate(zip(pt_bins[:-1], pt_bins[1:])):
        pt_min_rel = pt_min * 10
        pt_max_rel = pt_max * 10

        # Import binned data
        histoname = f"hMass_{pt_min_rel}_{pt_max_rel}"
        data = DataHandler(data=input_file, histoname=histoname, limits=[2.2, 2.4], rebin=4)

        # Check if data is binned
        if not data.get_is_binned():
            continue

        # Define PDFs lists
        signal_pdfs = ["gaussian"]
        background_pdfs = ["chebpol3"]

        # Define the fitter
        fitter = F2MassFitter(data_handler=data,
                              name_signal_pdf=signal_pdfs,
                              name_background_pdf=background_pdfs,
                              name=f"{background_pdfs[0]}_{signal_pdfs[0]}")
                              
        if fix_sigma:
            file = ROOT.TFile.Open(fix_file)
            fixed_sigma = file.Get("hRawYieldsSigma").GetBinContent(idx+1)
            fitter.set_signal_initpar(0, "sigma", fixed_sigma, fix=True)
        else:
            fitter.set_signal_initpar(0, "sigma", 0.017, limits=[0.005, 0.03])

        # Initiate fit parameters
        fitter.set_particle_mass(0, pdg_id=4122, limits=[2.27, 2.3])
        fitter.set_signal_initpar(0, "frac", 0.0001)
        fitter.set_background_initpar(0, "c0", 2)
        fitter.set_background_initpar(0, "c1", -0.1)
        fitter.set_background_initpar(0, "c2", -0.01)

        # Fit
        fitter.mass_zfit()

        # Plot the mass fit
        plot_mass_fit = fitter.plot_mass_fit(style="ATLAS", show_extra_info=False, extra_info_loc=['upper left', 'lower right'], axis_title=r"$M_{pK\pi} (GeV/c^{2})$")
        plot_mass_fit.text(0.195, 0.83, r'$\Lambda_c^+ \rightarrow \mathrm{p} \mathrm{K}^- \mathrm{\pi}^+ + \mathrm{c.c.}$')
        plot_mass_fit.text(0.195, 0.79, r'$\mathrm{pp}, \sqrt{\it{s}} = 13.6 \: \mathrm{TeV}$')
        plot_mass_fit.text(0.195, 0.75, f'{pt_min} < $p_{{\mathrm{{T}}}} < {pt_max} \: \mathrm{{GeV/c}}$')

        # Get fit results with errors
        mean, mean_error = fitter.get_mass()
        sigma, sigma_error = fitter.get_sigma()
        rawyield, rawyield_error = fitter.get_raw_yield(0)
        significance, signif_error = fitter.get_significance()
        bkg, bkg_error = fitter.get_background()
        chi2ndf = fitter.get_chi2_ndf()
        soverb, soverb_error = fitter.get_signal_over_background()
        purity = rawyield / (rawyield + bkg)

        # Estimate correlation coefficient
        rho = rawyield_error / bkg_error
        if rho > 1.:
            rho = 1. / rho
    
        # Calculate purity and its uncertainty
        purity = rawyield / (rawyield + bkg)
        purity_error = np.sqrt(
        (rawyield_error**2 / (rawyield + bkg)**2) +
        (rawyield**2 * bkg_error**2 / (rawyield + bkg)**4) -
        (2 * rho * rawyield * rawyield_error * bkg_error / (rawyield + bkg)**3)
        )

        latex_expr_mass = f"$\mu=$" + str(round(mean, 4)) + f"$\pm$" + str(round(mean_error, 4)) + "$\;\mathrm{GeV}/c^2$"
        latex_expr_width = f"$\sigma=$" + str(round(sigma, 3)) + f"$\pm$" + str(round(sigma_error, 3)) + "$\;\mathrm{GeV}/c^2$"
        latex_expr_rawyield = r"$N_{\Lambda_{\mathrm{c}}^{+}}\:=\:$" + str(int(rawyield)) + f"$\:\pm\:$" + str(int(rawyield_error))
        latex_expr_significance = r"$s/\sqrt{s+b}\:(3\sigma)\:=\:$" + str(round(significance, 1)) + f"$\:\pm\:$" + str(round(signif_error, 1))
        latex_expr_purity = f"Purity (s/(s+b))$\:=\:$"+str(round(purity,3))+f"$\:\pm\:$"+str(round(purity_error,3))

        plot_mass_fit.text(0.62, 0.67, latex_expr_mass, fontsize=13)
        plot_mass_fit.text(0.62, 0.62, latex_expr_width, fontsize=13)
        #plot_mass_fit.text(0.62, 0.57, latex_expr_purity, fontsize=13)
        plot_mass_fit.text(0.195, 0.40, latex_expr_rawyield, fontsize=15)
        plot_mass_fit.text(0.195, 0.35, latex_expr_significance, fontsize=15)
        plot_mass_fit.text(0.195, 0.30, latex_expr_purity, fontsize=15)

        # Save the plot for each pT bin
       
        plot_mass_fit_name = f'output/rawyield/flare/pdf/fitSpectrum_pT_{pt_min}_{pt_max}_file_{file_index}.pdf'
        plot_mass_fit.savefig(plot_mass_fit_name)
        pdf_files.append(plot_mass_fit_name)

        # Save residuals
        plot_raw_residuals = fitter.plot_raw_residuals(style="ATLAS", axis_title=r"$M_{pK\pi}$")
        plot_raw_residuals_name = f'output/rawyield/flare/pdf/fitResiduals_pT_{pt_min}_{pt_max}_file_{file_index}.pdf'
        plot_raw_residuals.savefig(plot_raw_residuals_name)
        pdf_files.append(plot_raw_residuals_name)

        # Fill ROOT histograms with fit results and errors
        bin_index = idx + 1  # ROOT bin index starts from 1

        mass_hist.SetBinContent(bin_index, mean)
        mass_hist.SetBinError(bin_index, mean_error)

        sigma_hist.SetBinContent(bin_index, sigma)
        sigma_hist.SetBinError(bin_index, sigma_error)

        raw_yield_hist.SetBinContent(bin_index, rawyield)
        raw_yield_hist.SetBinError(bin_index, rawyield_error)

        bkg_hist.SetBinContent(bin_index, bkg)
        bkg_hist.SetBinError(bin_index, bkg_error)

        significance_hist.SetBinContent(bin_index, significance)
        significance_hist.SetBinError(bin_index, signif_error)
        
        chi2_hist.SetBinContent(bin_index, chi2ndf)
        soverb_hist.SetBinContent(bin_index, soverb)
        soverb_hist.SetBinError(bin_index, soverb_error)

        purity_hist.SetBinContent(bin_index, purity)
        purity_hist.SetBinError(bin_index, purity_error)

    # %% Write fit results histograms to a separate ROOT file for each file_index
    root_output_file = ROOT.TFile(f"output/rawyield/flare/RawYields_Lc_run3_pp_13.6TeV_for_femto_dec18.root", "RECREATE")
    mass_hist.Write()
    sigma_hist.Write()
    raw_yield_hist.Write()
    bkg_hist.Write()
    significance_hist.Write()
    chi2_hist.Write()
    soverb_hist.Write()
    purity_hist.Write()
    root_output_file.Close()

    # Clear the histograms after writing to avoid conflicts
    mass_hist.Reset()
    sigma_hist.Reset()
    raw_yield_hist.Reset()
    significance_hist.Reset()
    bkg_hist.Reset()
    chi2_hist.Reset()
    soverb_hist.Reset()
    purity_hist.Reset()
    # Create a single PDF file from the generated PDFs
    with PdfWriter() as pdf_writer:
        for pdf_file in pdf_files:
            pdf_writer.append(pdf_file)

    # Save the combined PDF
    combined_pdf_file = f"output/rawyield/flare/RawYields_Lc_run3_pp_13.6TeV_for_femto_dec18.pdf"
    pdf_writer.write(combined_pdf_file)

    
    print(f"Fit results with errors saved to RawYields_Lc_run3_pp_13.6TeV{file_index}.root")


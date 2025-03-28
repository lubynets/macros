import ROOT
import pandas as pd

kBlueC = ROOT.TColor.GetColor('#1f78b4')
kOrangeC = ROOT.TColor.GetColor('#ff7f00')

def setHistStyle(hist, colour, marker=20, fillstyle=0, linewidth=1):
    hist.SetMarkerColor(colour)
    hist.SetLineColor(colour)
    hist.SetMarkerStyle(marker)
    hist.SetFillStyle(fillstyle)
    hist.SetLineWidth(linewidth)

def setCanvasStyle(canvas):
    canvas.SetTicks(1,1)

def fill_th1_hist(h, df, var):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var_val in df[var]:
        h.Fill(var_val)

def fill_th1_hist_abs(h, df, var):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var_val in df[var]:
        h.Fill(abs(var_val))

def fill_th2_hist(h, df, var1, var2):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var1_val, var2_val in zip(df[var1], df[var2]):
        h.Fill(var1_val, var2_val)

def fill_th2_hist_abs(h, df, var1, var2):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var1_val, var2_val in zip(df[var1], df[var2]):
        h.Fill(abs(var1_val), var2_val)

def fill_res_hist(h, df, var1, var2):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var_val1, var_val2 in zip(df[var1], df[var2]):
        h.Fill((var_val1 - var_val2)/var_val1)

def fill_th2_res_hist(h, df, var1, var2):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var_val1, var_val2 in zip(df[var1], df[var2]):
        h.Fill(var_val1, (var_val2 - var_val1)/var_val1)

def set_style():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptDate(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetLabelSize(0.04, 'xyz')
    ROOT.gStyle.SetTitleSize(0.05, 'xyz')
    ROOT.gStyle.SetTitleFont(42, 'xyz')
    ROOT.gStyle.SetLabelFont(42, 'xyz')
    ROOT.gStyle.SetTitleOffset(1.05, 'x')
    ROOT.gStyle.SetTitleOffset(1.1, 'y')
    ROOT.gStyle.SetCanvasDefW(800)
    ROOT.gStyle.SetCanvasDefH(600)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadGridX(0)
    ROOT.gStyle.SetPadGridY(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetPaperSize(20, 24)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetLegendFillColor(0)
    ROOT.gStyle.SetEndErrorSize(0.)
    ROOT.gStyle.SetMarkerSize(1)

def convert_sel_to_string(selection):
    sel_string = ''
    conj = ' & '
    for _, val in selection.items():
        sel_string = sel_string + '(' + val + ')' + conj
    return sel_string[:-len(conj)]

def calculate_statistics_for_th2(hist):
    # Calculate underflow and overflow for X and Y
    x_underflow = int(hist.Integral(0, 0, 1, hist.GetNbinsY()))
    x_overflow  = int(hist.Integral(hist.GetNbinsX() + 1, hist.GetNbinsX() + 1, 1, hist.GetNbinsY()))
    y_underflow = int(hist.Integral(1, hist.GetNbinsX(), 0, 0))
    y_overflow  = int(hist.Integral(1, hist.GetNbinsX(), hist.GetNbinsY() + 1, hist.GetNbinsY() + 1))

    return [f"Entries: {int(hist.GetEntries())}",
            f"X Mean: {round(hist.GetMean(1),3)}", f"X RMS: {round(hist.GetRMS(1),3)}",
            f"X Underflow: {x_underflow}", f"X Overflow: {x_overflow}",
            f"Y Mean: {round(hist.GetMean(2),3)}", f"Y RMS: {round(hist.GetRMS(2),3)}",
            f"Y Underflow: {y_underflow}", f"Y Overflow: {y_overflow}"]

def setStatsStyle(stats):
    stats.SetFillColor(0)
    stats.SetTextColor(ROOT.kBlack)
    stats.SetTextSize(0.025)
    stats.SetBorderSize(1)
    stats.SetTextAlign(12)
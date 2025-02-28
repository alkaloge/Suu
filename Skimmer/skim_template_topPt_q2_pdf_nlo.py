"""
File: skim_template_topPt_q2_pdf_nlo.py
Author: Alexis Kalogeropoulos
Email: 
Description: Skimmer that calculates q2/pdf/topPtRewighting/nlo_theory correction factors for the Suu analysis
"""

import ROOT as r
import numpy as np
import awkward as ak
import sys


def LoadHisto(file, hist_name):
    """
    Load a histogram from a ROOT file, clone it, and print debug information.
    """
    hist = file.Get(hist_name)
    if not hist:
        print(f"Warning: Histogram {hist_name} not found in file {file.GetName()}")
        return None
    
    # Clone the histogram to ensure it remains valid after the file is closed
    hist_clone = hist.Clone()
    hist_clone.SetDirectory(0)  # Detach from the file

    # Debug information
    print(f"Successfully loaded and cloned histogram {hist_name} from file {file.GetName()}")
    print(f"  Number of bins: {hist_clone.GetNbinsX()}")
    print(f"  Sum of weights: {hist_clone.GetSumOfWeights()}")
    print(f"  Integral: {hist_clone.Integral()}")
    print(f"  Mean: {hist_clone.GetMean()}")
    print(f"  RMS: {hist_clone.GetRMS()}")
    
    return hist_clone




def GetTheoryCorr(pt, IOV, correction_histograms, process, debug=False):
    """
    Calculate the theory correction factor using ROOT histograms.
    process: "wjets", "zjets", "dy", "znn", etc.
    """
    IOV = str(IOV)
    
    # Determine if the process is DY (zjets or Znn)
    is_DY = ('zjets' in process.lower()  or 'znn' in process.lower() )

    # Get the correction histograms for the given process
    if process == "wjets":
        ewk_hist = correction_histograms.get("w_ewk", None)
        qcd_hist = correction_histograms.get("w_qcd", None)
        qcd_ewk_hist = correction_histograms.get("w_qcd_ewk", None)
    elif is_DY:  # Handle zjets and Znn as part of DY
        ewk_hist = correction_histograms.get("z_ewk", None)
        qcd_hist = correction_histograms.get("z_qcd", None)
        qcd_ewk_hist = correction_histograms.get("z_qcd_ewk", None)

    #### --->>>>What is the difference here? it seems to be calling different histograms for 2017! but it should be that 2018 or 2017 are the same, 2016 should be different
    elif process == "dy":
        qcd_hist = correction_histograms.get("dy_qcd_2017", None)
    elif process == "znn":
        qcd_hist = correction_histograms.get("znn_qcd_2017", None)
    else:
        raise ValueError(f"Unknown process: {process}")
    if debug : 
        # Debugging print statement
        print(f"===================Debug Info for process {process}:")
        print(f"  ewk_hist: {ewk_hist}")
        print(f"  qcd_hist: {qcd_hist}")
        print(f"  qcd_ewk_hist: {qcd_ewk_hist}")
        print(f"  pt: {pt}")
        print(f"  IOV: {IOV}")

    ## I am not sure what is the difference between the dy and the zjets process...in any case for the DY_HT we should use "zjets"
    # Calculate the correction factor
    if process == "wjets":
        if ewk_hist and qcd_hist:
            ewk_corr = ewk_hist.GetBinContent(ewk_hist.FindBin(pt))
            if '2016' in IOV or 'preV' in IOV:
                qcd_corr = qcd_hist.GetBinContent(qcd_hist.FindBin(pt))
            else:
                qcd_corr = qcd_ewk_hist.GetBinContent(qcd_ewk_hist.FindBin(pt)) if qcd_ewk_hist else 1.0
            return ewk_corr * qcd_corr
        else:
            return 1.0  # Default value if histograms are missing
    elif is_DY:  # Handle zjets and Znn as part of DY
        if ewk_hist and qcd_hist:
            ewk_corr = ewk_hist.GetBinContent(ewk_hist.FindBin(pt))
            if '2016' in IOV or 'preV' in IOV:
                qcd_corr = qcd_hist.GetBinContent(qcd_hist.FindBin(pt))
            else:
                qcd_corr = qcd_ewk_hist.GetBinContent(qcd_ewk_hist.FindBin(pt)) if qcd_ewk_hist else 1.0
            return ewk_corr * qcd_corr
        else:
            return 1.0  # Default value if histograms are missing
    elif process == "dy":
        if qcd_hist:
            return qcd_hist.GetBinContent(qcd_hist.FindBin(pt))
        else:
            return 1.0  # Default value if histogram is missing
    elif process == "znn":
        if qcd_hist:
            return qcd_hist.GetBinContent(qcd_hist.FindBin(pt))
        else:
            return 1.0  # Default value if histogram is missing
    else:
        return 1.0  # Default value for unknown processes



def GetTheoryCorrr(pt, IOV, correction_histograms, process):
    """
    Calculate the theory correction factor using ROOT histograms.
    process: "wjets", "zjets", "dy", "znn", etc.
    """
    pt = 300+pt
    IOV = str(IOV)
    
    # Determine if the process is DY (zjets or Znn)
    is_DY = (process == "zjets" or process == "znn")

    # Get the correction histograms for the given process
    if process == "wjets":
        ewk_hist = correction_histograms.get("w_ewk", None)
        qcd_hist = correction_histograms.get("w_qcd", None)
        qcd_ewk_hist = correction_histograms.get("w_qcd_ewk", None)
    elif is_DY:  # Handle zjets and Znn as part of DY
        ewk_hist = correction_histograms.get("z_ewk", None)
        qcd_hist = correction_histograms.get("z_qcd", None)
        qcd_ewk_hist = correction_histograms.get("z_qcd_ewk", None)
    elif process == "dy":
        qcd_hist = correction_histograms.get("dy_qcd_2017", None)
    elif process == "znn":
        qcd_hist = correction_histograms.get("znn_qcd_2017", None)
    else:
        raise ValueError(f"Unknown process: {process}")

    # Debugging print statement
    print('Debug Info:', ewk_hist, qcd_hist, pt, process, correction_histograms)

    # Calculate the correction factor
    if process == "wjets":
        if ewk_hist and qcd_hist:
            ewk_corr = ewk_hist.GetBinContent(ewk_hist.FindBin(pt))
            if '2016' in IOV or 'preV' in IOV:
                qcd_corr = qcd_hist.GetBinContent(qcd_hist.FindBin(pt))
            else:
                qcd_corr = qcd_ewk_hist.GetBinContent(qcd_ewk_hist.FindBin(pt)) if qcd_ewk_hist else 1.0
            return ewk_corr * qcd_corr
        else:
            return 1.0  # Default value if histograms are missing
    elif is_DY:  # Handle zjets and Znn as part of DYa
        print('this is ok,', is_DY, process, ewk_hist)
        if ewk_hist and qcd_hist:
            ewk_corr = ewk_hist.GetBinContent(ewk_hist.FindBin(pt))
            if '2016' in IOV or 'preV' in IOV:
                qcd_corr = qcd_hist.GetBinContent(qcd_hist.FindBin(pt))
            else:
                qcd_corr = qcd_ewk_hist.GetBinContent(qcd_ewk_hist.FindBin(pt)) if qcd_ewk_hist else 1.0
            return ewk_corr * qcd_corr
        else:
            return 1.0  # Default value if histograms are missing
    elif process == "dy":
        if qcd_hist:
            return qcd_hist.GetBinContent(qcd_hist.FindBin(pt))
        else:
            return 1.0  # Default value if histogram is missing
    elif process == "znn":
        if qcd_hist:
            return qcd_hist.GetBinContent(qcd_hist.FindBin(pt))
        else:
            return 1.0  # Default value if histogram is missing
    else:
        return 1.0  # Default value for unknown processes





def GetPDFweights(LHEPdfWeight, var="nominal"):
    ## Determines the pdf up and down variations
    if len(LHEPdfWeight.shape) == 1:  # If LHEPdfWeight is a 1D array
        pdf = 1.0
        pdfUnc = LHEPdfWeight.std() / LHEPdfWeight.mean() if LHEPdfWeight.mean() != 0 else 0.0
    else:
        pdf = np.ones(len(LHEPdfWeight))
        pdfUnc = np.zeros(len(LHEPdfWeight))
        mean_pdf = np.mean(LHEPdfWeight, axis=1)
        nonzero_mean_mask = mean_pdf != 0
        pdfUnc[nonzero_mean_mask] = np.std(LHEPdfWeight[nonzero_mean_mask], axis=1) / mean_pdf[nonzero_mean_mask]

    if var == "up":
        pdf += pdfUnc
    elif var == "down":
        pdf -= pdfUnc

    return pdf




def GetQ2Weights(LHEScaleWeight, var="nominal"):
    ## Determines the envelope of the muR/muF up and down variations
    q2 = []
    q2Up = []
    q2Down = []
    
    if len(LHEScaleWeight.shape) == 1:
        LHEScaleWeight = [LHEScaleWeight]  # Convert to list of lists
        
    for weight in LHEScaleWeight:
        if len(weight) == 9:
            nom = weight[4]
            scales = [weight[i] for i in [0, 1, 3, 5, 7, 8]]
            q2Up.append(max(scales) / nom)
            q2Down.append(min(scales) / nom)
        elif len(weight) == 8:
            scales = [weight[i] for i in [0, 1, 3, 4, 6, 7]]
            q2Up.append(max(scales))
            q2Down.append(min(scales))
        else:
            q2Up.append(1.0)
            q2Down.append(1.0)
            
        q2.append(1.0)
    
    if var == "up":
        return q2Up
    elif var == "down":
        return q2Down
    else:
        return q2



def GetTopPtWeight(t_pt, tbar_pt, var="nominal"):
    # Top pt reweighting according to the recommendation:
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Case_3_1_Analyses_with_SM_tt_as
    if var == "nominal":
        return np.sqrt(np.exp(0.0615 - 0.0005 * t_pt) * np.exp(0.0615 - 0.0005 * tbar_pt))
    elif var == "up": return 1.
    else:
        return np.exp(0.0615 - 0.0005 * t_pt) * np.exp(0.0615 - 0.0005 * tbar_pt)





# Open the nanoAOD file
input_file = r.TFile.Open("inFileDY.root")
input_tree = input_file.Get("Events")
num_events = input_tree.GetEntries()


#input_tree.Scan("genWeight")
# Create a new output ROOT file
output_file = r.TFile.Open("selected_events.root", "RECREATE")

isMC=int(sys.argv[1])
IOV=str(sys.argv[2])
process=str(sys.argv[3])


nlo_corr = np.zeros(1, dtype=np.float32)
nlo_corr_norm = np.zeros(1, dtype=np.float32)


input_tree.SetBranchStatus("*", 0)  # Activate the branch
keys = ['CorrT1METJet*',   
        'Electron*',
        'Flag*',
        'HLT*',
        'Jet*',
        'L1*',
        'MET*',
        'Muon*',
        'PV*',
        'PuppiMET*',  
        'RawMET*',
        'RawPuppiMET*',
        'SV*',
        'event*',
        'fixedGridRhoFast*',
        'luminosityBlock*',
        'run*',
        'GenJet*',
        'GenPart*',
        'GenMET*',
        'btagWeight*',
        'Pileup*',
        'genWeight*',
        'PSWeight*',
        'LHE*',
        'TrigO*'
]


# Activate only the desired branches in the output tree
for key in keys:
    input_tree.SetBranchStatus(key, 1)  # Activate the branch

# Copy the selected branches to the output tree
#output_tree = input_tree.CopyTree("")
output_tree = input_tree.CloneTree(0)  # Clone the structure of the input tree



# Define arrays to store PDF, Q2, and top pT variations
pdf_N = np.zeros(1, dtype=float)
pdf_U = np.zeros(1, dtype=float)
pdf_D = np.zeros(1, dtype=float)

q2_N = np.zeros(1, dtype=float)
q2_U = np.zeros(1, dtype=float)
q2_D = np.zeros(1, dtype=float)

topPtWeight_N = np.zeros(1, dtype=float)
topPtWeight_U = np.zeros(1, dtype=float)
topPtWeight_W = np.zeros(1, dtype=float)

# Arrays to store normalized weights
pdf_N_norm = np.zeros(1, dtype=float)
pdf_U_norm = np.zeros(1, dtype=float)
pdf_D_norm = np.zeros(1, dtype=float)

q2_N_norm = np.zeros(1, dtype=float)
q2_U_norm = np.zeros(1, dtype=float)
q2_D_norm = np.zeros(1, dtype=float)

topPtWeight_N_norm = np.zeros(1, dtype=float)
topPtWeight_U_norm = np.zeros(1, dtype=float)
topPtWeight_W_norm = np.zeros(1, dtype=float)

# Create branches for unnormalized and normalized weights
output_tree.Branch("pdf_N", pdf_N, "pdf_N/D")
output_tree.Branch("pdf_U", pdf_U, "pdf_U/D")
output_tree.Branch("pdf_D", pdf_D, "pdf_D/D")

output_tree.Branch("q2_N", q2_N, "q2_N/D")
output_tree.Branch("q2_U", q2_U, "q2_U/D")
output_tree.Branch("q2_D", q2_D, "q2_D/D")

output_tree.Branch("topPtWeight_N", topPtWeight_N, "topPtWeight_N/D")
output_tree.Branch("topPtWeight_U", topPtWeight_U, "topPtWeight_U/D")
output_tree.Branch("topPtWeight_W", topPtWeight_W, "topPtWeight_W/D")

output_tree.Branch("pdf_N_norm", pdf_N_norm, "pdf_N_norm/D")
output_tree.Branch("pdf_U_norm", pdf_U_norm, "pdf_U_norm/D")
output_tree.Branch("pdf_D_norm", pdf_D_norm, "pdf_D_norm/D")

output_tree.Branch("q2_N_norm", q2_N_norm, "q2_N_norm/D")
output_tree.Branch("q2_U_norm", q2_U_norm, "q2_U_norm/D")
output_tree.Branch("q2_D_norm", q2_D_norm, "q2_D_norm/D")

output_tree.Branch("topPtWeight_N_norm", topPtWeight_N_norm, "topPtWeight_N_norm/D")
output_tree.Branch("topPtWeight_U_norm", topPtWeight_U_norm, "topPtWeight_U_norm/D")
output_tree.Branch("topPtWeight_W_norm", topPtWeight_W_norm, "topPtWeight_W_norm/D")


# Create a branch in the output tree for the theory correction
output_tree.Branch("nlo_corr", nlo_corr, "nlo_corr/F")
output_tree.Branch("nlo_corr_norm", nlo_corr_norm, "nlo_corr_norm/F")
nlo_corr[0] = 1.0  # Default value if no Z/W boson or data
nlo_corr_norm[0] = 1.0  # Default value if no Z/W boson or data

# Variables to store the sum of weights for normalization
Nev_pdf_N = 0.0
Nev_pdf_U = 0.0
Nev_pdf_D = 0.0

Nev_q2_N = 0.0
Nev_q2_U = 0.0
Nev_q2_D = 0.0

Nev_topPtWeight_N = 0.0
Nev_topPtWeight_U = 0.0
Nev_topPtWeight_W = 0.0


Nev_nlo_corr = 0.0

hist_EventCountGenW = r.TH1F("hist_EventCountGenW", "Event Cuts gen weights", 1,0,1)

cuts = ["Inclusive", "LeptonCut", "LeadingJetCut", "SubLeadingJetCut"]
num_cuts = len(cuts)
hist_EventCount = r.TH1F("hist_EventCount", "Event Cuts", num_cuts, 1, num_cuts+1)
hist_pdf_N = r.TH1F("hist_pdf_N", "PDF_N after Cuts", num_cuts, 1, num_cuts+1)
hist_pdf_U = r.TH1F("hist_pdf_U", "PDF_U after Cuts", num_cuts, 1, num_cuts+1)
hist_pdf_D = r.TH1F("hist_pdf_D", "PDF_D after Cuts", num_cuts, 1, num_cuts+1)

hist_q2_N = r.TH1F("hist_q2_N", "Q2_N after Cuts", num_cuts, 1, num_cuts+1)
hist_q2_U = r.TH1F("hist_q2_U", "Q2_U after Cuts", num_cuts, 1, num_cuts+1)
hist_q2_D = r.TH1F("hist_q2_D", "Q2_D after Cuts", num_cuts, 1, num_cuts+1)



hist_pdf_N_norm = r.TH1F("hist_pdf_N_norm", "PDF_N after Cuts norm", num_cuts, 1, num_cuts+1)
hist_pdf_U_norm = r.TH1F("hist_pdf_U_norm", "PDF_U after Cuts norm", num_cuts, 1, num_cuts+1)
hist_pdf_D_norm = r.TH1F("hist_pdf_D_norm", "PDF_D after Cuts norm", num_cuts, 1, num_cuts+1)

hist_q2_N_norm = r.TH1F("hist_q2_N_norm", "Q2_N after Cuts norm", num_cuts, 1, num_cuts+1)
hist_q2_U_norm = r.TH1F("hist_q2_U_norm", "Q2_U after Cuts norm", num_cuts, 1, num_cuts+1)
hist_q2_D_norm = r.TH1F("hist_q2_D_norm", "Q2_D after Cuts norm", num_cuts, 1, num_cuts+1)

hist_nlo_corr = r.TH1F("hist_nlo_corr", "NLO_Corr after Cuts", num_cuts, 1, num_cuts + 1)
hist_nlo_corr_norm = r.TH1F("hist_nlo_corr_norm", "NLO_Corr after Cuts norm", num_cuts, 1, num_cuts + 1)

hist_mu_pt = r.TH1F("hist_mu_pt", "mu_pt ", 30, 0, 300 )
hist_jet_pt = r.TH1F("hist_jet_pt", "jets_pt ", 30, 0, 300 )
hist_el_pt = r.TH1F("hist_el_pt", "el_pt ", 30, 0, 300 )
hist_N_mu = r.TH1F("hist_N_mu", "N mu after Cuts", num_cuts, 1, num_cuts+1)
#hist_N_mutotal = r.TH1F("hist_N_mutotal", "N mu before Cuts", num_cuts, 1, num_cuts+1)
hist_N_el = r.TH1F("hist_N_el", "N el after Cuts", num_cuts, 1, num_cuts+1)
#hist_N_el = r.TH1F("hist_N_el", "N el after Cuts", num_cuts, 1, num_cuts+1)



cut_bins = {
    1: "Inclusive",
    2: "LeptonCut",
    3: "LeadingJetCut",
    4: "SubLeadingJetCut"
}

for bin_number, bin_label in list(cut_bins.items()):
    hist_pdf_N.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_pdf_U.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_pdf_D.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_N.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_U.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_D.GetXaxis().SetBinLabel(bin_number, bin_label)

    hist_pdf_N_norm.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_pdf_U_norm.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_pdf_D_norm.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_N_norm.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_U_norm.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_D_norm.GetXaxis().SetBinLabel(bin_number, bin_label)

    hist_nlo_corr.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_nlo_corr_norm.GetXaxis().SetBinLabel(bin_number, bin_label)

    hist_EventCount.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_N_mu.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_N_el.GetXaxis().SetBinLabel(bin_number, bin_label)

# top pT reweighting code

# Add new branches for GenPart and top pt weights
genPart_pdgId = np.zeros(200, dtype=np.int32)  # Use 32-bit integers
genPart_pt = np.zeros(200, dtype=np.float32)
genPart_status = np.zeros(200, dtype=np.int32)  # Use 32-bit integers

input_tree.SetBranchAddress("GenPart_pdgId", genPart_pdgId)
input_tree.SetBranchAddress("GenPart_pt", genPart_pt)
input_tree.SetBranchAddress("GenPart_status", genPart_status)


# Create histograms for top pt weights
hist_topPt_N = r.TH1F("hist_topPt_N", "TopPt_N after Cuts", num_cuts, 1, num_cuts+1)
hist_topPt_U = r.TH1F("hist_topPt_U", "TopPt_U after Cuts", num_cuts, 1, num_cuts+1)
hist_topPt_W = r.TH1F("hist_topPt_W", "TopPt_Weight after Cuts", num_cuts, 1, num_cuts+1)
hist_topPt_N_norm = r.TH1F("hist_topPt_N_norm", "TopPt_N after Cuts norm", num_cuts, 1, num_cuts+1)
hist_topPt_U_norm = r.TH1F("hist_topPt_U_norm", "TopPt_U after Cuts norm", num_cuts, 1, num_cuts+1)
hist_topPt_W_norm = r.TH1F("hist_topPt_W_norm", "TopPt_Weight after Cuts norm", num_cuts, 1, num_cuts+1)

# Set bin labels for histograms
for bin_number, bin_label in cut_bins.items():
    hist_topPt_N.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_topPt_U.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_topPt_W.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_topPt_N_norm.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_topPt_U_norm.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_topPt_W_norm.GetXaxis().SetBinLabel(bin_number, bin_label)



Nev_LeptonCut=0
Nev_LeadingJetCut=0
Nev_SubLeadingJetCut=0
Nev_All=0
Nev_AllGenW=0.
num_events = input_tree.GetEntries()

'''
user_input = sys.argv[1]

if user_input == 'A':
    start_index = 0
    end_index = num_events // 3
elif user_input == 'B':
    start_index = num_events // 3
    end_index = 2 * num_events // 3
elif user_input == 'C':
    start_index = 2 * num_events // 3
    end_index = num_events+1
else:
    print("Invalid input. Please provide 'A', 'B', or 'C'.")
'''

start_index=0
end_index =  num_events
end_index=2000

print(('will skim from ', start_index, 'to event', end_index))
# Iterate over the selected range of events

# First loop: Calculate the sum of weights for normalization
correction_histograms={}
if isMC:

    # Load histograms for wjets
    file_wjets = r.TFile("./merged_kfactors_wjets.root")
    if file_wjets.IsOpen():
        print("Successfully opened file: ./merged_kfactors_wjets.root")
        correction_histograms["w_ewk"] = LoadHisto(file_wjets, "kfactor_monojet_ewk")
        correction_histograms["w_qcd"] = LoadHisto(file_wjets, "kfactor_monojet_qcd")
        correction_histograms["w_qcd_ewk"] = LoadHisto(file_wjets, "kfactor_monojet_qcd_ewk")
        file_wjets.Close()
    else:
        print("Error: Failed to open file: ./merged_kfactors_wjets.root")

    # Load histograms for zjets
    file_zjets = r.TFile("./merged_kfactors_zjets.root")
    if file_zjets.IsOpen():
        print("Successfully opened file: ./merged_kfactors_zjets.root")
        correction_histograms["z_ewk"] = LoadHisto(file_zjets, "kfactor_monojet_ewk")
        correction_histograms["z_qcd"] = LoadHisto(file_zjets, "kfactor_monojet_qcd")
        correction_histograms["z_qcd_ewk"] = LoadHisto(file_zjets, "kfactor_monojet_qcd_ewk")
        file_zjets.Close()
    else:
        print("Error: Failed to open file: ./merged_kfactors_zjets.root")

    # Load histograms for dy and znn
    file_dy = r.TFile("./kfac_dy_filter.root")
    if file_dy.IsOpen():
        print("Successfully opened file: ./kfac_dy_filter.root")
        correction_histograms["dy_qcd_2017"] = LoadHisto(file_dy, "kfac_dy_filter")
        file_dy.Close()
    else:
        print("Error: Failed to open file: ./kfac_dy_filter.root")

    file_znn = r.TFile("./kfac_znn_filter.root")
    if file_znn.IsOpen():
        print("Successfully opened file: ./kfac_znn_filter.root")
        correction_histograms["znn_qcd_2017"] = LoadHisto(file_znn, "kfac_znn_filter")
        file_znn.Close()
    else:
        print("Error: Failed to open file: ./kfac_znn_filter.root")


    #for i_event in range(num_events):
    for i_event in range(start_index, end_index):
        input_tree.GetEntry(i_event)
        if i_event%100==0 : print(("processing event", i_event))
        
        # Calculate PDF weights
        try:
            LHEPdfWeight = np.array(input_tree.LHEPdfWeight)
            pdf_N[0] = np.mean(GetPDFweights(LHEPdfWeight))
            pdf_U[0] = np.mean(GetPDFweights(LHEPdfWeight, var="up"))
            pdf_D[0] = np.mean(GetPDFweights(LHEPdfWeight, var="down"))
        except AttributeError:
            LHEPdfWeight = np.ones(1)
            pdf_N[0] = 1.0
            pdf_U[0] = 1.0
            pdf_D[0] = 1.0

        # Calculate Q2 weights
        try:
            LHEScaleWeight_ak = np.array(input_tree.LHEScaleWeight)
            q2_N[0] = np.mean(GetQ2Weights(LHEScaleWeight_ak))
            q2_U[0] = np.mean(GetQ2Weights(LHEScaleWeight_ak, var="up"))
            q2_D[0] = np.mean(GetQ2Weights(LHEScaleWeight_ak, var="down"))
        except AttributeError:
            LHEScaleWeight_ak = np.ones(1)
            q2_N[0] = 1.0
            q2_U[0] = 1.0
            q2_D[0] = 1.0

        # Calculate top pT weights
        foundTop = False
        foundAntitop = False
        top_pt = 0.0
        antitop_pt = 0.0

        for i_gen in range(min(input_tree.nGenPart, len(genPart_pdgId))):
            if genPart_pdgId[i_gen] == 6 and genPart_status[i_gen] == 62:
                top_pt = genPart_pt[i_gen]
                foundTop = True
            if genPart_pdgId[i_gen] == -6 and genPart_status[i_gen] == 62:
                antitop_pt = genPart_pt[i_gen]
                foundAntitop = True

        if foundTop and foundAntitop:
            topPtWeight_N[0] = GetTopPtWeight(top_pt, antitop_pt, "nominal")
            topPtWeight_U[0] = GetTopPtWeight(top_pt, antitop_pt, "up")
            topPtWeight_W[0] = GetTopPtWeight(top_pt, antitop_pt, "weight")
        else:
            topPtWeight_N[0] = 1.0
            topPtWeight_U[0] = 1.0
            topPtWeight_W[0] = 1.0


        # Calculate NLO correction
        if process == 'zjets' or process == 'wjets':
            gen_pt = input_tree.GenPart_pt
            gen_pdgId = input_tree.GenPart_pdgId
            boson_pt = 0.0
            for i in range(len(gen_pdgId)):
                if process == "zjets" and abs(gen_pdgId[i]) == 23:  # Z boson
                    boson_pt = gen_pt[i]
                    break
                elif process == "wjets" and abs(gen_pdgId[i]) == 24:  # W boson
                    boson_pt = gen_pt[i]
                    break
            if boson_pt > 0:
                nlo_corr[0] = GetTheoryCorr(boson_pt, IOV, correction_histograms, process)
                if nlo_corr[0] == 0.0:
                    nlo_corr[0] = 1.0  # Default value if correction is zero
            else:
                nlo_corr[0] = 1.0  # Default value if no Z/W boson is found

            # Accumulate sum for normalization
            Nev_nlo_corr += nlo_corr[0]

        # Accumulate sums for normalization
        Nev_pdf_N += pdf_N[0]
        Nev_pdf_U += pdf_U[0]
        Nev_pdf_D += pdf_D[0]

        Nev_q2_N += q2_N[0]
        Nev_q2_U += q2_U[0]
        Nev_q2_D += q2_D[0]

        Nev_topPtWeight_N += topPtWeight_N[0]
        Nev_topPtWeight_U += topPtWeight_U[0]
        Nev_topPtWeight_W += topPtWeight_W[0]




# Second loop: Normalize the top pT weights and store them in branches


print(('will skim from ', start_index, 'to event', end_index))
# Iterate over the selected range of events
for i_event in range(start_index, end_index):



# Loop over events and apply selection
#for i_event in range(input_tree.GetEntries()):
#for i_event in range(min(1000, input_tree.GetEntries())):
    Nm=0
    Ne=0
    counter=1
    
    if i_event%5000==0 : print(("processing event", i_event))
    input_tree.GetEntry(i_event)

    Nev_All +=1
    if isMC: Nev_AllGenW +=input_tree.genWeight
    else : Nev_AllGenW +=1


    nlo_corr[0] = 1.0  # Default value if no Z/W boson or data
    if isMC:
        # Get the generator-level Z boson pt (or W boson pt, depending on your analysis)
        # Replace "GenPart_pt" and "GenPart_pdgId" with the appropriate branch names
        # Find the Z boson (pdgId = 23) or W boson (pdgId = 24)
        if process == 'zjets' or process == 'wjets' :
            gen_pt = input_tree.GenPart_pt
            gen_pdgId = input_tree.GenPart_pdgId
        
            boson_pt = 0.0
            for i in range(len(gen_pdgId)):
                if process == "zjets" and abs(gen_pdgId[i]) == 23:  # Z boson
                    boson_pt = gen_pt[i]
                    break
                elif process == "wjets" and abs(gen_pdgId[i]) == 24:  # W boson
                    boson_pt = gen_pt[i]
                    break
            
            # Calculate the theory correction factor
            if boson_pt > 0:
                nlo_corr[0] = GetTheoryCorr(boson_pt, IOV, correction_histograms, process)
                if nlo_corr[0] == 0. : nlo_corr[0] = 1.
                #print (nlo_corr[0])
            else:
                nlo_corr[0] = 1.0  # Default value if no Z/W boson is found
            
            # Normalize the NLO correction
            nlo_corr_norm[0] = nlo_corr[0] / Nev_nlo_corr if Nev_nlo_corr > 0 else 1.0

        # Fill histograms


    #LHEScaleWeight_ak = ak.Array([list(LHEScaleWeight) for _ in range(len(LHEScaleWeight))])
    #LHEScaleWeight_ak = ak.from_iter(LHEScaleWeight)
    if isMC:
        try : LHEScaleWeight_ak = np.array(input_tree.LHEScaleWeight)
        except AttributeError : LHEScaleWeight_ak = np.ones(1)
    
    ''' 
    # Get Q2 weights
    q2_nominal=[]
    q2_up=[]
    q2_down=[]
    q2_nominal.append(1)
    q2_up.append(1)
    q2_down.append(1)
    if isMC:
        q2_nominal = GetQ2Weights(LHEScaleWeight_ak)
        q2_up = GetQ2Weights(LHEScaleWeight_ak, var="up")
        q2_down = GetQ2Weights(LHEScaleWeight_ak, var="down")


    q2_N[0] = q2_nominal[0]
    q2_U[0] = q2_up[0]
    q2_D[0] = q2_down[0]

    ''' 

    # Calculate PDF weights
    if isMC:
        try:
            LHEPdfWeight = np.array(input_tree.LHEPdfWeight)
            pdf_N[0] = np.mean(GetPDFweights(LHEPdfWeight))
            pdf_U[0] = np.mean(GetPDFweights(LHEPdfWeight, var="up"))
            pdf_D[0] = np.mean(GetPDFweights(LHEPdfWeight, var="down"))
        except AttributeError:
            LHEPdfWeight = np.ones(1)
            pdf_N[0] = 1.0
            pdf_U[0] = 1.0
            pdf_D[0] = 1.0

    # Calculate Q2 weights
    if isMC:
        try:
            LHEScaleWeight_ak = np.array(input_tree.LHEScaleWeight)
            q2_N[0] = np.mean(GetQ2Weights(LHEScaleWeight_ak))
            q2_U[0] = np.mean(GetQ2Weights(LHEScaleWeight_ak, var="up"))
            q2_D[0] = np.mean(GetQ2Weights(LHEScaleWeight_ak, var="down"))
        except AttributeError:
            LHEScaleWeight_ak = np.ones(1)
            q2_N[0] = 1.0
            q2_U[0] = 1.0
            q2_D[0] = 1.0

    # Calculate top pT weights
    foundTop = False
    foundAntitop = False
    top_pt = 0.0
    antitop_pt = 0.0

    if isMC:
        for i_gen in range(min(input_tree.nGenPart, len(genPart_pdgId))):
            if genPart_pdgId[i_gen] == 6 and genPart_status[i_gen] == 62:
                top_pt = genPart_pt[i_gen]
                foundTop = True
            if genPart_pdgId[i_gen] == -6 and genPart_status[i_gen] == 62:
                antitop_pt = genPart_pt[i_gen]
                foundAntitop = True

        if foundTop and foundAntitop:
            topPtWeight_N[0] = GetTopPtWeight(top_pt, antitop_pt, "nominal")
            topPtWeight_U[0] = GetTopPtWeight(top_pt, antitop_pt, "up")
            topPtWeight_W[0] = GetTopPtWeight(top_pt, antitop_pt, "weight")
        else:
            topPtWeight_N[0] = 1.0
            topPtWeight_U[0] = 1.0
            topPtWeight_W[0] = 1.0


    # Normalize the weights
    pdf_N_norm[0] = pdf_N[0] / Nev_pdf_N
    pdf_U_norm[0] = pdf_U[0] / Nev_pdf_U
    pdf_D_norm[0] = pdf_D[0] / Nev_pdf_D

    q2_N_norm[0] = q2_N[0] / Nev_q2_N
    q2_U_norm[0] = q2_U[0] / Nev_q2_U
    q2_D_norm[0] = q2_D[0] / Nev_q2_D

    topPtWeight_N_norm[0] = topPtWeight_N[0] / Nev_topPtWeight_N
    topPtWeight_U_norm[0] = topPtWeight_U[0] / Nev_topPtWeight_U
    topPtWeight_W_norm[0] = topPtWeight_W[0] / Nev_topPtWeight_W



    # Fill histogram bins with pdf_N after each cut

    hist_q2_N_norm.Fill(counter, q2_N_norm[0])
    hist_q2_U_norm.Fill(counter, q2_U_norm[0])
    hist_q2_D_norm.Fill(counter, q2_D_norm[0])
    hist_pdf_N_norm.Fill(counter, pdf_N_norm[0])  # 
    hist_pdf_U_norm.Fill(counter, pdf_U_norm[0])  # 
    hist_pdf_D_norm.Fill(counter, pdf_D_norm[0])  # 
    hist_topPt_N_norm.Fill(counter, topPtWeight_N_norm[0])
    hist_topPt_U_norm.Fill(counter, topPtWeight_U_norm[0])
    hist_topPt_W_norm.Fill(counter, topPtWeight_W_norm[0])

    hist_q2_N.Fill(counter, q2_N[0])
    hist_q2_U.Fill(counter, q2_U[0])
    hist_q2_D.Fill(counter, q2_D[0])
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    hist_topPt_N.Fill(counter, topPtWeight_N[0])
    hist_topPt_U.Fill(counter, topPtWeight_U[0])
    hist_topPt_W.Fill(counter, topPtWeight_W[0])

    hist_nlo_corr.Fill(counter, nlo_corr[0])  # Fill for all events
    hist_nlo_corr_norm.Fill(counter, nlo_corr_norm[0])  # Fill for all events

    hist_N_mu.Fill(counter, input_tree.nMuon)
    hist_N_el.Fill(counter, input_tree.nElectron)
    #topPt
    counter+=1
 
    # Loop over Muons
    for j in range(input_tree.nMuon): 
        if  input_tree.Muon_tightId[j] and input_tree.Muon_highPtId[j] > 1 and abs(input_tree.Muon_eta[j]) < 2.4:

            hist_mu_pt.Fill(input_tree.Muon_pt[j])
            #print i_event, range(input_tree.nMuon), abs(input_tree.Muon_eta[j]) < 2.4, input_tree.Muon_eta[j]
        #print 'muons...', input_tree.Muon_pt[j], input_tree.Muon_eta[j]
        if  input_tree.Muon_pt[j] > 50.:
            Nm+=1
    
    #for value in muon_pt:
    #    print(value)

    # Loop over Electrons
    for j in range(input_tree.nElectron): 
        if input_tree.Electron_mvaFall17V2noIso_WP80[j] and abs(input_tree.Electron_eta[j]) < 2.5:
            hist_el_pt.Fill(input_tree.Electron_pt[j])
        if input_tree.Electron_pt[j] > 100:
            Ne+=1

    if (Ne + Nm) < 1 : continue
    Nev_LeptonCut +=1

    hist_q2_N_norm.Fill(counter, q2_N_norm[0])
    hist_q2_U_norm.Fill(counter, q2_U_norm[0])
    hist_q2_D_norm.Fill(counter, q2_D_norm[0])
    hist_pdf_N_norm.Fill(counter, pdf_N_norm[0])  # 
    hist_pdf_U_norm.Fill(counter, pdf_U_norm[0])  # 
    hist_pdf_D_norm.Fill(counter, pdf_D_norm[0])  # 
    hist_topPt_N_norm.Fill(counter, topPtWeight_N_norm[0])
    hist_topPt_U_norm.Fill(counter, topPtWeight_U_norm[0])
    hist_topPt_W_norm.Fill(counter, topPtWeight_W_norm[0])

    hist_q2_N.Fill(counter, q2_N[0])
    hist_q2_U.Fill(counter, q2_U[0])
    hist_q2_D.Fill(counter, q2_D[0])
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    hist_topPt_N.Fill(counter, topPtWeight_N[0])
    hist_topPt_U.Fill(counter, topPtWeight_U[0])
    hist_topPt_W.Fill(counter, topPtWeight_W[0])
    hist_nlo_corr.Fill(counter, nlo_corr[0])  # Fill for all events
    hist_nlo_corr_norm.Fill(counter, nlo_corr_norm[0])  # Fill for all events

    hist_N_mu.Fill(counter, input_tree.nMuon)
    hist_N_el.Fill(counter, input_tree.nElectron)
    #topPt
    counter+=1


    nJet25 = 0
    nJet150 = 0
    for j in range(input_tree.nJet):
        hist_jet_pt.Fill(input_tree.Jet_pt[j])
        if input_tree.Jet_pt[j] > 25:
            nJet25 += 1
            if input_tree.Jet_pt[j] > 150:
                nJet150 += 1


    if nJet150 < 1 : continue
    Nev_LeadingJetCut += 1


    hist_q2_N_norm.Fill(counter, q2_N_norm[0])
    hist_q2_U_norm.Fill(counter, q2_U_norm[0])
    hist_q2_D_norm.Fill(counter, q2_D_norm[0])
    hist_pdf_N_norm.Fill(counter, pdf_N_norm[0])  # 
    hist_pdf_U_norm.Fill(counter, pdf_U_norm[0])  # 
    hist_pdf_D_norm.Fill(counter, pdf_D_norm[0])  # 
    hist_topPt_N_norm.Fill(counter, topPtWeight_N_norm[0])
    hist_topPt_U_norm.Fill(counter, topPtWeight_U_norm[0])
    hist_topPt_W_norm.Fill(counter, topPtWeight_W_norm[0])

    hist_q2_N.Fill(counter, q2_N[0])
    hist_q2_U.Fill(counter, q2_U[0])
    hist_q2_D.Fill(counter, q2_D[0])
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    hist_topPt_N.Fill(counter, topPtWeight_N[0])
    hist_topPt_U.Fill(counter, topPtWeight_U[0])
    hist_topPt_W.Fill(counter, topPtWeight_W[0])
    hist_nlo_corr.Fill(counter, nlo_corr[0])  # Fill for all events
    hist_nlo_corr_norm.Fill(counter, nlo_corr_norm[0])  # Fill for all events

    hist_N_mu.Fill(counter, input_tree.nMuon)
    hist_N_el.Fill(counter, input_tree.nElectron)
    counter+=1

    if nJet25 < 2 : continue
    Nev_SubLeadingJetCut +=1


    hist_q2_N_norm.Fill(counter, q2_N_norm[0])
    hist_q2_U_norm.Fill(counter, q2_U_norm[0])
    hist_q2_D_norm.Fill(counter, q2_D_norm[0])
    hist_pdf_N_norm.Fill(counter, pdf_N_norm[0])  # 
    hist_pdf_U_norm.Fill(counter, pdf_U_norm[0])  # 
    hist_pdf_D_norm.Fill(counter, pdf_D_norm[0])  # 
    hist_topPt_N_norm.Fill(counter, topPtWeight_N_norm[0])
    hist_topPt_U_norm.Fill(counter, topPtWeight_U_norm[0])
    hist_topPt_W_norm.Fill(counter, topPtWeight_W_norm[0])

    hist_q2_N.Fill(counter, q2_N[0])
    hist_q2_U.Fill(counter, q2_U[0])
    hist_q2_D.Fill(counter, q2_D[0])
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    hist_topPt_N.Fill(counter, topPtWeight_N[0])
    hist_topPt_U.Fill(counter, topPtWeight_U[0])
    hist_topPt_W.Fill(counter, topPtWeight_W[0])
    hist_nlo_corr.Fill(counter, nlo_corr[0])  # Fill for all events
    hist_nlo_corr_norm.Fill(counter, nlo_corr_norm[0])  # Fill for all events

    hist_N_mu.Fill(counter, input_tree.nMuon)
    hist_N_el.Fill(counter, input_tree.nElectron)
    counter+=1




    output_tree.Fill()

    #print 'counters', counter, Nev_All, Nev_LeptonCut, Nev_LeadingJetCut, Nev_SubLeadingJetCut, Nev_AllGenW, nmuon, nelectron


hist_EventCount.Fill(1,Nev_All)
hist_EventCount.Fill(2,Nev_LeptonCut)
hist_EventCount.Fill(3,Nev_LeadingJetCut)
hist_EventCount.Fill(4,Nev_SubLeadingJetCut)
hist_EventCountGenW.Fill(0,Nev_AllGenW)




output_file.Write()
output_file.Close()

# Close the input file
input_file.Close()


# Close the correction histogram files
#if process == "zjets":
#    zjets_file.Close()
#    dy_filter_file.Close()
#elif process == "wjets":
#    wjets_file.Close()



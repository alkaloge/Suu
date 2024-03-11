import ROOT as r
import numpy as np
import awkward as ak



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





# Open the nanoAOD file
input_file = r.TFile.Open("inFile.root")
input_tree = input_file.Get("Events")
num_events = input_tree.GetEntries()


#input_tree.Scan("genWeight")
# Create a new output ROOT file
output_file = r.TFile.Open("selected_events.root", "RECREATE")
output_tree = input_tree.CloneTree(0)  # Clone the structure of the input tree


# Define arrays to store PDF variations
pdf_N = np.zeros(1, dtype=float)
pdf_U = np.zeros(1, dtype=float)
pdf_D = np.zeros(1, dtype=float)

# Create branches in the output tree
output_tree.Branch("pdf_N", pdf_N, "pdf_N/D")
output_tree.Branch("pdf_U", pdf_U, "pdf_U/D")
output_tree.Branch("pdf_D", pdf_D, "pdf_D/D")



# Create arrays to hold branch data
muon_pt = np.zeros(10, dtype=float)  # Assuming there are at most 10 muons per event
muon_eta = np.zeros(10, dtype=float)
muon_tightId = np.zeros(10, dtype=int)
muon_highPtId = np.zeros(10, dtype=int)
electron_pt = np.zeros(10, dtype=float)  # Assuming there are at most 10 electrons per event
electron_eta = np.zeros(10, dtype=float)
electron_mvaFall17V2noIso_WP80 = np.zeros(10, dtype=bool)
jet_pt = np.zeros(20, dtype=float)  # Assuming there are at most 20 jets per event
#genW = np.zeros(1, dtype=float)  # Assuming there are at most 20 jets per event

#LHEScaleWeight = r.std.vector(float)()


# Set branch addresses for input tree
input_tree.SetBranchAddress("Muon_pt", muon_pt)
input_tree.SetBranchAddress("Muon_eta", muon_eta)
input_tree.SetBranchAddress("Muon_tightId", muon_tightId)
input_tree.SetBranchAddress("Muon_highPtId", muon_highPtId)
input_tree.SetBranchAddress("Electron_pt", electron_pt)
input_tree.SetBranchAddress("Electron_eta", electron_eta)
input_tree.SetBranchAddress("Electron_mvaFall17V2noIso_WP80", electron_mvaFall17V2noIso_WP80)
input_tree.SetBranchAddress("Jet_pt", jet_pt)
#input_tree.SetBranchAddress("LHEScaleWeight", LHEScaleWeight)


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


cut_bins = {
    1: "Inclusive",
    2: "LeptonCut",
    3: "LeadingJetCut",
    4: "SubLeadingJetCut"
}

for bin_number, bin_label in cut_bins.items():
    hist_pdf_N.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_pdf_U.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_pdf_D.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_N.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_U.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_D.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_EventCount.GetXaxis().SetBinLabel(bin_number, bin_label)

Nev_LeptonCut=0
Nev_LeadingJetCut=0
Nev_SubLeadingJetCut=0
Nev_All=0
Nev_AllGenW=0.

# Loop over events and apply selection
#for i_event in range(input_tree.GetEntries()):
for i_event in range(min(10000, input_tree.GetEntries())):
    Nm=0
    Ne=0
    counter=1
    
    if i_event%1000==0 : print "processing event", i_event
    input_tree.GetEntry(i_event)

    Nev_All +=1
    Nev_AllGenW +=input_tree.genWeight

    #LHEScaleWeight_ak = ak.Array([list(LHEScaleWeight) for _ in range(len(LHEScaleWeight))])
    #LHEScaleWeight_ak = ak.from_iter(LHEScaleWeight)
    LHEScaleWeight_ak = np.array(input_tree.LHEScaleWeight)
    
    
    # Get Q2 weights
    q2_nominal=1
    q2_up=1
    q2_down=1
    q2_nominal = GetQ2Weights(LHEScaleWeight_ak)
    q2_up = GetQ2Weights(LHEScaleWeight_ak, var="up")
    q2_down = GetQ2Weights(LHEScaleWeight_ak, var="down")
    hist_q2_N.Fill(counter, q2_nominal[0])
    hist_q2_U.Fill(counter, q2_up[0])
    hist_q2_D.Fill(counter, q2_down[0])

    # Print or do whatever you want with the Q2 weights
    #print("Event:", i_event)
    #print("Nominal Q2 weights:", q2_nominal)
    #print("Q2 up weights:", q2_up)
    #print("Q2 down weights:", q2_down)

    LHEPdfWeight = np.array(input_tree.LHEPdfWeight)
    pdf_N[0] = np.mean(GetPDFweights(LHEPdfWeight))
    pdf_U[0] = np.mean(GetPDFweights(LHEPdfWeight, var="up"))
    pdf_D[0] = np.mean(GetPDFweights(LHEPdfWeight, var="down"))
    # Fill histogram bins with pdf_N after each cut
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    
    counter+=1
 
    # Loop over Muons
    for i_muon in range(len(muon_pt)):
        if i_muon >= input_tree.nMuon:
            break  # Break if no more muons in this event
        if (muon_tightId[i_muon] | (muon_highPtId[i_muon] > 1)) & (muon_pt[i_muon] > 50) & (np.abs(muon_eta[i_muon]) < 2.4):
            #output_tree.Fill()
            Nm+=1
    
    # Loop over Electrons
    for i_electron in range(len(electron_pt)):
        if i_electron >= input_tree.nElectron:
            break  # Break if no more electrons in this event
        if (electron_mvaFall17V2noIso_WP80[i_electron]) & (electron_pt[i_electron] > 100) & (np.abs(electron_eta[i_electron]) < 2.5):
            #output_tree.Fill()
            Ne+=1

    if (Ne + Nm) < 1 : continue
    Nev_LeptonCut +=1

    q2_nominal=1
    q2_up=1
    q2_down=1
    q2_nominal = GetQ2Weights(LHEScaleWeight_ak)
    q2_up = GetQ2Weights(LHEScaleWeight_ak, var="up")
    q2_down = GetQ2Weights(LHEScaleWeight_ak, var="down")
    hist_q2_N.Fill(counter, q2_nominal[0])
    hist_q2_U.Fill(counter, q2_up[0])
    hist_q2_D.Fill(counter, q2_down[0])


    LHEPdfWeight = np.array(input_tree.LHEPdfWeight)
    pdf_N[0] = np.mean(GetPDFweights(LHEPdfWeight))
    pdf_U[0] = np.mean(GetPDFweights(LHEPdfWeight, var="up"))
    pdf_D[0] = np.mean(GetPDFweights(LHEPdfWeight, var="down"))
    # Fill histogram bins with pdf_N after each cut
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    counter+=1

    nJet25 = 0
    nJet100 = 0
    for i_jet in range(len(jet_pt)):
	if i_jet >= input_tree.nJet:
	    break  # Break if no more jets in this event
	if jet_pt[i_jet] > 25:
	    nJet25 += 1
	    if jet_pt[i_jet] > 100:
		nJet100 += 1


    if nJet100 < 1 : continue
    Nev_LeadingJetCut += 1

    q2_nominal=1
    q2_up=1
    q2_down=1
    q2_nominal = GetQ2Weights(LHEScaleWeight_ak)
    q2_up = GetQ2Weights(LHEScaleWeight_ak, var="up")
    q2_down = GetQ2Weights(LHEScaleWeight_ak, var="down")

    hist_q2_N.Fill(counter, q2_nominal[0])
    hist_q2_U.Fill(counter, q2_up[0])
    hist_q2_D.Fill(counter, q2_down[0])


    LHEPdfWeight = np.array(input_tree.LHEPdfWeight)
    pdf_N[0] = np.mean(GetPDFweights(LHEPdfWeight))
    pdf_U[0] = np.mean(GetPDFweights(LHEPdfWeight, var="up"))
    pdf_D[0] = np.mean(GetPDFweights(LHEPdfWeight, var="down"))
    # Fill histogram bins with pdf_N after each cut
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    counter+=1

    if nJet25 < 2 : continue
    Nev_SubLeadingJetCut +=1


    q2_nominal=1
    q2_up=1
    q2_down=1
    q2_nominal = GetQ2Weights(LHEScaleWeight_ak)
    q2_up = GetQ2Weights(LHEScaleWeight_ak, var="up")
    q2_down = GetQ2Weights(LHEScaleWeight_ak, var="down")


    hist_q2_N.Fill(counter, q2_nominal[0])
    hist_q2_U.Fill(counter, q2_up[0])
    hist_q2_D.Fill(counter, q2_down[0])

    LHEPdfWeight = np.array(input_tree.LHEPdfWeight)
    pdf_N[0] = np.mean(GetPDFweights(LHEPdfWeight))
    pdf_U[0] = np.mean(GetPDFweights(LHEPdfWeight, var="up"))
    pdf_D[0] = np.mean(GetPDFweights(LHEPdfWeight, var="down"))
    # Fill histogram bins with pdf_N after each cut
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    counter+=1

    output_tree.Fill()

    #print 'counters', Nev_All, Nev_LeptonCut, Nev_LeadingJetCut, Nev_SubLeadingJetCut, Nev_AllGenW


hist_EventCount.Fill(1,Nev_All)
hist_EventCount.Fill(2,Nev_LeptonCut)
hist_EventCount.Fill(3,Nev_LeadingJetCut)
hist_EventCount.Fill(4,Nev_SubLeadingJetCut)
hist_EventCountGenW.Fill(0,Nev_AllGenW)

output_file.Write()
output_file.Close()

# Close the input file
input_file.Close()





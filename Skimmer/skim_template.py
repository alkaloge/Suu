import ROOT as r
import numpy as np
import awkward as ak
import sys


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

isMC=int(sys.argv[1])


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
]


# Activate only the desired branches in the output tree
for key in keys:
    input_tree.SetBranchStatus(key, 1)  # Activate the branch

# Copy the selected branches to the output tree
#output_tree = input_tree.CopyTree("")
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

#LHEScaleWeight = r.std.vector(float)()


# Set branch addresses for input tree
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

for bin_number, bin_label in cut_bins.items():
    hist_pdf_N.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_pdf_U.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_pdf_D.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_N.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_U.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_D.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_EventCount.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_N_mu.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_N_el.GetXaxis().SetBinLabel(bin_number, bin_label)

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

print 'will skim from ', start_index, 'to event', end_index
# Iterate over the selected range of events
for i_event in range(start_index, end_index):



# Loop over events and apply selection
#for i_event in range(input_tree.GetEntries()):
#for i_event in range(min(1000, input_tree.GetEntries())):
    Nm=0
    Ne=0
    counter=1
    
    if i_event%5000==0 : print "processing event", i_event
    input_tree.GetEntry(i_event)

    Nev_All +=1
    if isMC: Nev_AllGenW +=input_tree.genWeight
    else : Nev_AllGenW +=1

    #LHEScaleWeight_ak = ak.Array([list(LHEScaleWeight) for _ in range(len(LHEScaleWeight))])
    #LHEScaleWeight_ak = ak.from_iter(LHEScaleWeight)
    if isMC:
        try : LHEScaleWeight_ak = np.array(input_tree.LHEScaleWeight)
        except AttributeError : LHEScaleWeight_ak = np.ones(1)
    
    
    # Get Q2 weights
    q2_nominal=[]
    q2_up=[]
    q2_down=[]
    q2_nominal.append(1)
    q2_up.append(1)
    q2_down.append(1)
    #q2_up=
    #q2_down=1
    if isMC:
	q2_nominal = GetQ2Weights(LHEScaleWeight_ak)
	q2_up = GetQ2Weights(LHEScaleWeight_ak, var="up")
	q2_down = GetQ2Weights(LHEScaleWeight_ak, var="down")

    #hist_q2_N.Fill(counter, q2_nominal[0])
    #hist_q2_U.Fill(counter, q2_up[0])
    #hist_q2_D.Fill(counter, q2_down[0])

    # Print or do whatever you want with the Q2 weights
    #print("Event:", i_event)
    #print("Nominal Q2 weights:", q2_nominal)
    #print("Q2 up weights:", q2_up)
    #print("Q2 down weights:", q2_down)
    if isMC :
       
	try : 
	    LHEPdfWeight = np.array(input_tree.LHEPdfWeight)
	    pdf_N[0] = np.mean(GetPDFweights(LHEPdfWeight))
	    pdf_U[0] = np.mean(GetPDFweights(LHEPdfWeight, var="up"))
	    pdf_D[0] = np.mean(GetPDFweights(LHEPdfWeight, var="down"))
        except AttributeError:
	    LHEPdfWeight = np.ones(1)
	    pdf_N[0] = 1.
	    pdf_U[0] = 1.
	    pdf_D[0] = 1.
    # Fill histogram bins with pdf_N after each cut
    #hist_pdf_N.Fill(counter, pdf_N[0])  # 
    #hist_pdf_U.Fill(counter, pdf_U[0])  # 
    #hist_pdf_D.Fill(counter, pdf_D[0])  # 

    hist_q2_N.Fill(counter, q2_nominal[0])
    hist_q2_U.Fill(counter, q2_up[0])
    hist_q2_D.Fill(counter, q2_down[0])
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    hist_N_mu.Fill(counter, input_tree.nMuon)
    hist_N_el.Fill(counter, input_tree.nElectron)
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
    hist_q2_N.Fill(counter, q2_nominal[0])
    hist_q2_U.Fill(counter, q2_up[0])
    hist_q2_D.Fill(counter, q2_down[0])
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    hist_N_mu.Fill(counter, input_tree.nMuon)
    hist_N_el.Fill(counter, input_tree.nElectron)
    counter+=1

    nJet25 = 0
    nJet100 = 0
    for j in range(input_tree.nJet):
        hist_jet_pt.Fill(input_tree.Jet_pt[j])
	if input_tree.Jet_pt[j] > 25:
	    nJet25 += 1
	    if input_tree.Jet_pt[j] > 100:
		nJet100 += 1


    if nJet100 < 1 : continue
    Nev_LeadingJetCut += 1

    hist_q2_N.Fill(counter, q2_nominal[0])
    hist_q2_U.Fill(counter, q2_up[0])
    hist_q2_D.Fill(counter, q2_down[0])
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    hist_N_mu.Fill(counter, input_tree.nMuon)
    hist_N_el.Fill(counter, input_tree.nElectron)
    counter+=1

    if nJet25 < 2 : continue
    Nev_SubLeadingJetCut +=1


    hist_q2_N.Fill(counter, q2_nominal[0])
    hist_q2_U.Fill(counter, q2_up[0])
    hist_q2_D.Fill(counter, q2_down[0])
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
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





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



def GetTopPtWeight(t_pt, tbar_pt, var="nominal"):
    # Top pt reweighting according to the recommendation:
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Case_3_1_Analyses_with_SM_tt_as
    #if var == "nominal":
    if var == "up":
        return np.sqrt(np.exp(0.0615 - 0.0005 * t_pt) * np.exp(0.0615 - 0.0005 * tbar_pt))
    elif var == "weight":
        return np.exp(0.0615 - 0.0005 * t_pt) * np.exp(0.0615 - 0.0005 * tbar_pt)
    else:
        return 1.0





# Open the nanoAOD file
input_file = r.TFile.Open("inFilett.root")
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
        'TrigO*'
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

for bin_number, bin_label in list(cut_bins.items()):
    hist_pdf_N.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_pdf_U.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_pdf_D.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_N.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_U.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_q2_D.GetXaxis().SetBinLabel(bin_number, bin_label)
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

topPtWeight_N = np.zeros(1, dtype=float)
topPtWeight_U = np.zeros(1, dtype=float)
topPtWeight_W = np.zeros(1, dtype=float)

output_tree.Branch("topPtWeight_N", topPtWeight_N, "topPtWeight_N/D")
output_tree.Branch("topPtWeight_U", topPtWeight_U, "topPtWeight_U/D")
output_tree.Branch("topPtWeight_W", topPtWeight_W, "topPtWeight_W/D")

# Create histograms for top pt weights
hist_topPt_N = r.TH1F("hist_topPt_N", "TopPt_N after Cuts", num_cuts, 1, num_cuts+1)
hist_topPt_U = r.TH1F("hist_topPt_U", "TopPt_U after Cuts", num_cuts, 1, num_cuts+1)
hist_topPt_W = r.TH1F("hist_topPt_W", "TopPt_Weight after Cuts", num_cuts, 1, num_cuts+1)

# Set bin labels for histograms
for bin_number, bin_label in cut_bins.items():
    hist_topPt_N.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_topPt_U.GetXaxis().SetBinLabel(bin_number, bin_label)
    hist_topPt_W.GetXaxis().SetBinLabel(bin_number, bin_label)




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
    foundTop = False
    foundAntitop = False
    top_pt = 0.0
    antitop_pt = 0.0
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


        # In the event loop:
        for i_gen in range(min(input_tree.nGenPart, len(genPart_pdgId))):
            print (i_gen, genPart_pdgId[i_gen], genPart_status[i_gen], genPart_pt[i_gen])
            if genPart_pdgId[i_gen] == 6 and genPart_status[i_gen] == 62:
                top_pt = genPart_pt[i_gen]
                foundTop = True
            if genPart_pdgId[i_gen] == -6 and genPart_status[i_gen] == 62:
                antitop_pt = genPart_pt[i_gen]
                foundAntitop = True


        if foundTop and foundAntitop:
            event_topPt_N = GetTopPtWeight(top_pt, antitop_pt, "nominal")
            event_topPt_U = GetTopPtWeight(top_pt, antitop_pt, "up")
            event_topPt_W = GetTopPtWeight(top_pt, antitop_pt, "weight")
        else:
            event_topPt_N = 1.0
            event_topPt_U = 1.0
            event_topPt_W = 1.0

        topPtWeight_N[0] = event_topPt_N
        topPtWeight_U[0] = event_topPt_U
        topPtWeight_W[0] = event_topPt_W


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
    #topPt
    hist_topPt_N.Fill(counter, event_topPt_N)
    hist_topPt_U.Fill(counter, event_topPt_U)
    hist_topPt_W.Fill(counter, event_topPt_W)
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
    #topPt
    hist_topPt_N.Fill(counter, event_topPt_N)
    hist_topPt_U.Fill(counter, event_topPt_U)
    hist_topPt_W.Fill(counter, event_topPt_W)
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

    hist_q2_N.Fill(counter, q2_nominal[0])
    hist_q2_U.Fill(counter, q2_up[0])
    hist_q2_D.Fill(counter, q2_down[0])
    hist_pdf_N.Fill(counter, pdf_N[0])  # 
    hist_pdf_U.Fill(counter, pdf_U[0])  # 
    hist_pdf_D.Fill(counter, pdf_D[0])  # 
    hist_N_mu.Fill(counter, input_tree.nMuon)
    hist_N_el.Fill(counter, input_tree.nElectron)
    #topPt
    hist_topPt_N.Fill(counter, event_topPt_N)
    hist_topPt_U.Fill(counter, event_topPt_U)
    hist_topPt_W.Fill(counter, event_topPt_W)
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
    #topPt
    hist_topPt_N.Fill(counter, event_topPt_N)
    hist_topPt_U.Fill(counter, event_topPt_U)
    hist_topPt_W.Fill(counter, event_topPt_W)
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





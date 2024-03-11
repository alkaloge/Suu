
import ROOT
import np as numpy
# Open the nAOD file
file = ROOT.TFile.Open("inFile.root")

# Access the Events tree
tree = file.Get("Events")

# Define the array to store Muon_pt values
muon_pts = []

for event in tree:
    # Access the Muon_pt value for the event
    muon_pt = event.Muon_pt
    # Loop over all muons in the event
    for pt in muon_pt:
        muon_pts.append(pt)

# Convert muon_pts to a NumPy array
muon_pts = np.array(muon_pts)

# Print the first 10 Muon_pt values
print(muon_pts[:10])

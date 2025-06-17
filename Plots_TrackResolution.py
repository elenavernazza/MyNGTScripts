import ROOT, os, math, sys, glob
import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import mplhep as hep
plt.style.use(hep.style.CMS)

files0 = ['/data/evernazz/2025_06_12/CMSSW_15_1_X_2025-06-11-2300/src/TestGenJets/Phase2_Scouting_HLT_inNANOAODSIM.root']
dataframe_files0 = ROOT.vector(str)()
for f in files0:
    file = ROOT.TFile.Open(f)
    if file.Get("Events"):
        dataframe_files0.push_back(f)
df0 = ROOT.RDataFrame("Events", dataframe_files0)

files1 = ['/data/evernazz/2025_06_12/CMSSW_15_1_0_pre3/src/TestCAExtension/Phase2_Scouting_HLT_inNANOAODSIM.root']
dataframe_files1 = ROOT.vector(str)()
for f in files1:
    file = ROOT.TFile.Open(f)
    if file.Get("Events"):
        dataframe_files1.push_back(f)
df1 = ROOT.RDataFrame("Events", dataframe_files1)

output_dir = "Plots_TrackResolution/Test_0"
os.system(f"mkdir -p {output_dir}")

##########################################################
### For each general track, find the closest pixel track
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;

    ROOT::RVec<int> match_general_pixel_tracks (
        Vfloat hltGeneralTrack_pt, Vfloat hltGeneralTrack_eta, Vfloat hltGeneralTrack_phi,
        Vfloat hltPixelTrack_pt, Vfloat hltPixelTrack_eta, Vfloat hltPixelTrack_phi
        ) {

        ROOT::RVec<int> match_general_pixel_tracks_indices (hltGeneralTrack_pt.size(), -1);
        std::vector<bool> pixel_track_used(hltPixelTrack_pt.size(), false);

        // std::cout << "Start event" << std::endl;
        for (size_t i_gen_track = 0; i_gen_track < hltGeneralTrack_pt.size(); i_gen_track ++) {
            auto gen_track_tlv = TLorentzVector();
            gen_track_tlv.SetPtEtaPhiM(hltGeneralTrack_pt.at(i_gen_track), hltGeneralTrack_eta.at(i_gen_track), hltGeneralTrack_phi.at(i_gen_track), 0);
            double matched_dR = std::numeric_limits<double>::infinity();
            int matched_i_pix_track = -1;
            for (size_t i_pix_track = 0; i_pix_track < hltPixelTrack_pt.size(); i_pix_track ++) {
                if (pixel_track_used[i_pix_track]) continue;  // Skip if already matched
                auto pix_track_tlv = TLorentzVector();
                pix_track_tlv.SetPtEtaPhiM(hltPixelTrack_pt.at(i_pix_track), hltPixelTrack_eta.at(i_pix_track), hltPixelTrack_phi.at(i_pix_track), 0);
                double dR = gen_track_tlv.DeltaR(pix_track_tlv);
                // std::cout << "General Track index " << i_gen_track 
                    // << " (pt = " << hltGeneralTrack_pt.at(i_gen_track) << ", eta = " << hltGeneralTrack_eta.at(i_gen_track) << ", phi = " << hltGeneralTrack_phi.at(i_gen_track)
                    // << "): Pixel Track index " << i_pix_track << " with dR = " << dR 
                    // << " (pt = " << hltPixelTrack_pt.at(i_pix_track) << ", eta = " << hltPixelTrack_eta.at(i_pix_track) << ", phi = " << hltPixelTrack_phi.at(i_pix_track) <<  ")" << std::endl;
                if (dR < matched_dR) {
                    matched_dR = dR;
                    matched_i_pix_track = static_cast<int>(i_pix_track);
                }
            }
            if (matched_i_pix_track != -1) {
                pixel_track_used[matched_i_pix_track] = true;
                match_general_pixel_tracks_indices[i_gen_track] = matched_i_pix_track;
            }
            // std::cout << "General Track index " << i_gen_track << ": matched to Pixel Track index " << matched_i_pix_track << " with dR = " << matched_dR << std::endl;
        }
        return match_general_pixel_tracks_indices;
    }
""")

df0 = df0.Define("matched_pixel_tracks_indeces",
               "match_general_pixel_tracks(hltGeneralTrack_pt, hltGeneralTrack_eta, hltGeneralTrack_phi, hltPixelTrack_pt, hltPixelTrack_eta, hltPixelTrack_phi)")
df1 = df1.Define("matched_pixel_tracks_indeces",
               "match_general_pixel_tracks(hltGeneralTrack_pt, hltGeneralTrack_eta, hltGeneralTrack_phi, hltPixelTrack_pt, hltPixelTrack_eta, hltPixelTrack_phi)")

ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;
    
    ROOT::RVec<float> get_matched_pixel_track_pts(Vfloat hltPixelTrack_var, Vint matched_pixel_tracks_indeces) {
        ROOT::RVec<float> matched_hltPixelTrack_var;
        for (size_t i_gen_track = 0; i_gen_track < matched_pixel_tracks_indeces.size(); i_gen_track ++) {
            int i_pixel_track = matched_pixel_tracks_indeces[i_gen_track];
            if (i_pixel_track != -1) { // && (i_pixel_track < hltPixelTrack_var.size())
                matched_hltPixelTrack_var.push_back(hltPixelTrack_var[i_pixel_track]);
            }
            else {
                matched_hltPixelTrack_var.push_back(-1);
            }
        }
        return matched_hltPixelTrack_var;
    }
""")

df0 = df0.Define("matched_pixel_tracks_mask", "(matched_pixel_tracks_indeces != -1)")
df0 = df0.Define("matched_hltPixelTrack_pt", "get_matched_pixel_track_pts(hltPixelTrack_pt, matched_pixel_tracks_indeces)[matched_pixel_tracks_mask]")
df0 = df0.Define("matched_hltPixelTrack_eta", "get_matched_pixel_track_pts(hltPixelTrack_eta, matched_pixel_tracks_indeces)[matched_pixel_tracks_mask]")
df0 = df0.Define("matched_hltPixelTrack_phi", "get_matched_pixel_track_pts(hltPixelTrack_phi, matched_pixel_tracks_indeces)[matched_pixel_tracks_mask]")
df0 = df0.Define("matched_hltGeneralTrack_pt", "hltGeneralTrack_pt[matched_pixel_tracks_mask]")
df0 = df0.Define("res_pixel_vs_gen_pt", "matched_hltPixelTrack_pt/matched_hltGeneralTrack_pt")

df1 = df1.Define("matched_pixel_tracks_mask", "(matched_pixel_tracks_indeces != -1)")
df1 = df1.Define("matched_hltPixelTrack_pt", "get_matched_pixel_track_pts(hltPixelTrack_pt, matched_pixel_tracks_indeces)[matched_pixel_tracks_mask]")
df1 = df1.Define("matched_hltPixelTrack_eta", "get_matched_pixel_track_pts(hltPixelTrack_eta, matched_pixel_tracks_indeces)[matched_pixel_tracks_mask]")
df1 = df1.Define("matched_hltPixelTrack_phi", "get_matched_pixel_track_pts(hltPixelTrack_phi, matched_pixel_tracks_indeces)[matched_pixel_tracks_mask]")
df1 = df1.Define("matched_hltGeneralTrack_pt", "hltGeneralTrack_pt[matched_pixel_tracks_mask]")
df1 = df1.Define("res_pixel_vs_gen_pt", "matched_hltPixelTrack_pt/matched_hltGeneralTrack_pt")

### Debug
# print(df1.AsNumpy(["res_pixel_vs_gen_pt"]))

TH1D_Res = ROOT.RDF.TH1DModel("h", "", 50, 0.0, 2.5)

h_pt0 = df0.Histo1D(TH1D_Res, "res_pixel_vs_gen_pt")
c_pt0 = np.array([h_pt0.GetBinContent(i) for i in range(1, h_pt0.GetNbinsX()+1)])
e_pt0 = np.array([h_pt0.GetBinLowEdge(i) for i in range(1, h_pt0.GetNbinsX()+2)])

h_pt1 = df1.Histo1D(TH1D_Res, "res_pixel_vs_gen_pt")
c_pt1 = np.array([h_pt1.GetBinContent(i) for i in range(1, h_pt1.GetNbinsX()+1)])
e_pt1 = np.array([h_pt1.GetBinLowEdge(i) for i in range(1, h_pt1.GetNbinsX()+2)])

fig, ax = plt.subplots(figsize=(10, 10))
hep.histplot(c_pt0, bins=e_pt0, label="Without CA Extension", ax=ax, histtype='step', linewidth=2)
hep.histplot(c_pt1, bins=e_pt1, label="With CA Extension", ax=ax, histtype='step', linewidth=2)
ax.set_xlabel("$p_T^{pixel}/p_T^{general}$ ")
ax.set_ylabel("Number of Tracks")
ax.legend(loc='upper right', fontsize=16, title='NGT Scouting (10 events TTbar PU 200)', title_fontsize=16)
plt.grid()
hep.cms.text("Work in Progress", ax=ax, fontsize=20)
hep.cms.lumitext(text='Phase-II (14 TeV)', ax=ax, fontsize=20)
print(f" ### INFO: Saving plot in {output_dir}/Resolution_pixel_vs_general_pT.png")
plt.savefig(f"{output_dir}/Resolution_pixel_vs_general_pT.png")
plt.savefig(f"{output_dir}/Resolution_pixel_vs_general_pT.pdf")


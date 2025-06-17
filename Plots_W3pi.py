import ROOT, os, math, sys, glob, pdb
import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import mplhep as hep
plt.style.use(hep.style.CMS)

from Functions_W3pi import *

#################################################################
# this version plots the W mass from gen level and from the HLT
# tracks (matching from gen to hlt, but using tracks!)
#################################################################

files = ['/data/evernazz/2025_05_15/CMSSW_14_0_9/src/PrivateProduction_1000ev/HLT_VertexInfo_inNANOAODSIM.root']

dataframe_files = ROOT.vector(str)()
for f in files:
    file = ROOT.TFile.Open(f)
    if file.Get("Events"):
        dataframe_files.push_back(f)
df = ROOT.RDataFrame("Events", dataframe_files)#.Range(0, 20) ############# [FIXME] 

output_dir = "Plots_W3pi_5"
os.system(f"mkdir -p {output_dir}")

#####################################################
# Plot dEta between between PF candidates and associated tracks: this demonstrates that there is a bug in hltPFCandidate_trackIndex!!!!!!!!!!

ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;
    ROOT::RVec<float> get_dEta_pf_track (Vint hltPFCandidate_charge, Vint hltPFCandidate_trackIndex, Vfloat hltPFCandidate_eta, Vfloat hltGeneralTrack_eta) {
        std::vector<float> dEta_pf_track;
        for (size_t i_hlt = 0; i_hlt < hltPFCandidate_charge.size(); i_hlt ++) {
            if (std::abs(hltPFCandidate_charge.at(i_hlt)) != 0) continue;
            // std::cout << "hltPFCandidate_eta idx " << i_hlt << " = " << hltPFCandidate_eta.at(i_hlt) 
                // << ", hltGeneralTrack_eta idx " << hltPFCandidate_trackIndex.at(i_hlt) << " = " << hltGeneralTrack_eta.at(hltPFCandidate_trackIndex.at(i_hlt)) << std::endl;
            auto dEta = hltPFCandidate_eta.at(i_hlt) - hltGeneralTrack_eta.at(hltPFCandidate_trackIndex.at(i_hlt));
            // std::cout << "dEta = " << dEta << std::endl;
            dEta_pf_track.push_back(dEta);
        }
        return dEta_pf_track;
    }
""")

df = df.Define("dEta_pf_track", "get_dEta_pf_track(hltPFCandidate_charge, hltPFCandidate_trackIndex, hltPFCandidate_eta, hltGeneralTrack_eta)")
h_hlt = df.Histo1D("dEta_pf_track")
c_hlt = np.array([h_hlt.GetBinContent(i) for i in range(1, h_hlt.GetNbinsX()+1)])
e_hlt = np.array([h_hlt.GetBinLowEdge(i) for i in range(1, h_hlt.GetNbinsX()+2)])

fig, ax = plt.subplots(figsize=(10, 10))
hep.histplot(c_hlt, bins=e_hlt, label=f"", ax=ax, histtype='step', linewidth=2)
ax.set_xlabel("$\Delta\eta$")
ax.set_ylabel("Events")
ax.legend(loc='upper left', fontsize=16)
plt.grid()
print(f" ### INFO: Saving plot {output_dir}/dEta_pf_track.png")
plt.savefig(f"{output_dir}/dEta_pf_track.png")
plt.savefig(f"{output_dir}/dEta_pf_track.pdf")

### Get 3 pions from W decay and sort them in pt
df = df.Define("gen_pion_indices", "get_sorted_gen_pion_idices(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother, GenPart_pt)")
df = df.Define("gen_pion_idx_0", "gen_pion_indices[0]").Define("gen_pion_idx_1", "gen_pion_indices[1]").Define("gen_pion_idx_2", "gen_pion_indices[2]")
df = df.Define("gen_W_mass", "get_W_mass(gen_pion_indices, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
h_gen = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 50, 0, 150), "gen_W_mass")
c_gen = np.array([h_gen.GetBinContent(i) for i in range(1, h_gen.GetNbinsX()+1)])
e_gen = np.array([h_gen.GetBinLowEdge(i) for i in range(1, h_gen.GetNbinsX()+2)])

### Filter events where one pion has eta > 4
df_acc = df.Filter("std::abs(GenPart_eta[gen_pion_idx_0]) < 4 && std::abs(GenPart_eta[gen_pion_idx_1]) < 4 && std::abs(GenPart_eta[gen_pion_idx_2]) < 4")
h_gen_filter = df_acc.Histo1D(ROOT.RDF.TH1DModel("h", "", 50, 0, 150), "gen_W_mass")
c_gen_filter = np.array([h_gen_filter.GetBinContent(i) for i in range(1, h_gen_filter.GetNbinsX()+1)])
e_gen_filter = np.array([h_gen_filter.GetBinLowEdge(i) for i in range(1, h_gen_filter.GetNbinsX()+2)])

#####################################################
### For each gen pion, among the hlt tracks get the one with closest dR
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;

    ROOT::RVec<int> match_gen_pion_to_hlt_track (Vint gen_pion_indices, Vfloat GenPart_pt, Vfloat GenPart_eta, Vfloat GenPart_phi, Vfloat GenPart_mass,
        Vfloat hltGeneralTrack_pt, Vfloat hltGeneralTrack_eta, Vfloat hltGeneralTrack_phi) {

        ROOT::RVec<int> matched_hlt_track_indices (gen_pion_indices.size(), -1);

        // std::cout << "Start event" << std::endl;
        for (size_t i_gen_pion = 0; i_gen_pion < gen_pion_indices.size(); i_gen_pion ++) {
            if (gen_pion_indices.at(i_gen_pion) == -1) continue; // Skip if no gen pion found
            auto gen_tlv = TLorentzVector();
            gen_tlv.SetPtEtaPhiM(GenPart_pt.at(gen_pion_indices.at(i_gen_pion)), GenPart_eta.at(gen_pion_indices.at(i_gen_pion)), GenPart_phi.at(gen_pion_indices.at(i_gen_pion)), GenPart_mass.at(gen_pion_indices.at(i_gen_pion)));
            double matched_dR = std::numeric_limits<double>::infinity();
            double matched_dPt = std::numeric_limits<double>::infinity();
            int matched_i_hlt = -1;
            for (size_t i_hlt = 0; i_hlt < hltGeneralTrack_pt.size(); i_hlt ++) {
                // if (hltGeneralTrack_pt.at(i_hlt) < 1.) continue;
                auto hlt_tlv = TLorentzVector();
                hlt_tlv.SetPtEtaPhiM(hltGeneralTrack_pt.at(i_hlt), hltGeneralTrack_eta.at(i_hlt), hltGeneralTrack_phi.at(i_hlt), 0.139526);
                double dR = gen_tlv.DeltaR(hlt_tlv);
                double dPt = std::abs(gen_tlv.Pt() - hlt_tlv.Pt());
                if ((dR < 0.1) && (dR < matched_dR)) {
                    matched_dR = dR;
                    matched_i_hlt = static_cast<int>(i_hlt);
                }
            }
            matched_hlt_track_indices[i_gen_pion] = matched_i_hlt;
            // std::cout << "Gen Jet index " << gen_pion_indices.at(i_gen_pion) << ": matched to hlt track index " << matched_i_hlt << " with dR = " << matched_dR << std::endl;
        }
        return matched_hlt_track_indices;
    }
""")

df = df.Define("hlt_track_indices",
               "match_gen_pion_to_hlt_track(gen_pion_indices, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, " \
               "hltGeneralTrack_pt, hltGeneralTrack_eta, hltGeneralTrack_phi)")
df = df.Define("hlt_track_idx_0", "hlt_track_indices[0]").Define("hlt_track_idx_1", "hlt_track_indices[1]").Define("hlt_track_idx_2", "hlt_track_indices[2]")

df = df.Define("gen_pion_pt0", "GenPart_pt[gen_pion_idx_0]").Define("gen_pion_pt1", "GenPart_pt[gen_pion_idx_1]").Define("gen_pion_pt2", "GenPart_pt[gen_pion_idx_2]") 
df = df.Define("hlt_pion_pt0", "hltGeneralTrack_pt[hlt_track_idx_0]").Define("hlt_pion_pt1", "hltGeneralTrack_pt[hlt_track_idx_1]").Define("hlt_pion_pt2", "hltGeneralTrack_pt[hlt_track_idx_2]") 
df = df.Define("gen_pion_eta0", "GenPart_eta[gen_pion_idx_0]").Define("gen_pion_eta1", "GenPart_eta[gen_pion_idx_1]").Define("gen_pion_eta2", "GenPart_eta[gen_pion_idx_2]")
df = df.Define("hlt_pion_eta0", "hltGeneralTrack_eta[hlt_track_idx_0]").Define("hlt_pion_eta1", "hltGeneralTrack_eta[hlt_track_idx_1]").Define("hlt_pion_eta2", "hltGeneralTrack_eta[hlt_track_idx_2]")

df = df.Define("hlt_W_mass", "get_W_mass(hlt_track_indices, hltGeneralTrack_pt, hltGeneralTrack_eta, hltGeneralTrack_phi, ROOT::RVec<float> (hltGeneralTrack_pt.size(), 0.139526))")
df_acc = df.Filter("std::abs(GenPart_eta[gen_pion_idx_0]) < 4 && std::abs(GenPart_eta[gen_pion_idx_1]) < 4 && std::abs(GenPart_eta[gen_pion_idx_2]) < 4")
df_acc_hlt_matched = df_acc.Filter("(hlt_track_idx_0 >= 0) && (hlt_track_idx_1 >= 0) && (hlt_track_idx_2 >= 0)")
h_hlt = df_acc_hlt_matched.Histo1D(ROOT.RDF.TH1DModel("h", "", 50, 0, 150), "hlt_W_mass")
c_hlt = np.array([h_hlt.GetBinContent(i) for i in range(1, h_hlt.GetNbinsX()+1)])
e_hlt = np.array([h_hlt.GetBinLowEdge(i) for i in range(1, h_hlt.GetNbinsX()+2)])

df_acc_hlt_matched_ditau = df_acc_hlt_matched.Filter("(gen_pion_pt0 > 35) && (gen_pion_pt1 > 35)")
h_hlt_ditau = df_acc_hlt_matched_ditau.Histo1D(ROOT.RDF.TH1DModel("h", "", 50, 0, 150), "hlt_W_mass")
c_hlt_ditau = np.array([h_hlt_ditau.GetBinContent(i) for i in range(1, h_hlt_ditau.GetNbinsX()+1)])
e_hlt_ditau = np.array([h_hlt_ditau.GetBinLowEdge(i) for i in range(1, h_hlt_ditau.GetNbinsX()+2)])

def GetLabelMass(histo):
    integral = histo.Integral()
    resolution = histo.GetRMS()/histo.GetMean()
    # return f" - {int(integral)} events : $\sigma/\mu=${resolution:.3f}"
    return f": $\sigma/\mu=${resolution:.3f}"

fig, ax = plt.subplots(figsize=(10, 10))
hep.histplot(c_gen, bins=e_gen, label=f"Gen-level {GetLabelMass(h_gen)}", ax=ax, histtype='step', linewidth=2)
hep.histplot(c_gen_filter, bins=e_gen_filter, label=f"Gen-level $|\eta|<4$ {GetLabelMass(h_gen_filter)}", ax=ax, histtype='step', linewidth=2)
hep.histplot(c_hlt, bins=e_hlt, label=f"HLT matched (No L1 filter) {GetLabelMass(h_hlt)}", ax=ax, histtype='step', linewidth=2)
ax.set_xlabel("W mass [GeV]")
ax.set_ylabel("Events")
ax.legend(loc='upper left', fontsize=18)
ax.set_ylim(0, 1.25*max(c_gen))
plt.grid()
hep.cms.text("Work in Progress", ax=ax, fontsize=20)
hep.cms.lumitext(text='Phase-II (14 TeV)', ax=ax, fontsize=20)
print(f" ### INFO: Saving plot {output_dir}/W_mass.png")
plt.savefig(f"{output_dir}/W_mass.png")
plt.savefig(f"{output_dir}/W_mass.pdf")

#####################################################
### Plot gen pions pt

def GetEfficencyError(h_num, h_den):
    h_unc = []
    for num, den in zip(h_num, h_den):
        if num != 0 and den != 0:
            h_unc.append(num/den * np.sqrt(1/num + 1/den)) 
        else:
            h_unc.append(1)
    return h_unc

pion_name = ['Leading', 'Second-leading', 'Third-leading']

for i in range(3):
    h_gen_pt = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 20, 0, 100), f"gen_pion_pt{i}")
    c_gen_pt = np.array([h_gen_pt.GetBinContent(j) for j in range(1, h_gen_pt.GetNbinsX()+1)])
    e_gen_pt = np.array([h_gen_pt.GetBinLowEdge(j) for j in range(1, h_gen_pt.GetNbinsX()+2)])

    h_gen_pt_filter = df_acc.Histo1D(ROOT.RDF.TH1DModel("h", "", 20, 0, 100), f"gen_pion_pt{i}")
    c_gen_pt_filter = np.array([h_gen_pt_filter.GetBinContent(j) for j in range(1, h_gen_pt_filter.GetNbinsX()+1)])
    e_gen_pt_filter = np.array([h_gen_pt_filter.GetBinLowEdge(j) for j in range(1, h_gen_pt_filter.GetNbinsX()+2)])

    h_hlt_pt = df_acc_hlt_matched.Histo1D(ROOT.RDF.TH1DModel("h", "", 20, 0, 100), f"hlt_pion_pt{i}")
    c_hlt_pt = np.array([h_hlt_pt.GetBinContent(j) for j in range(1, h_hlt_pt.GetNbinsX()+1)])
    e_hlt_pt = np.array([h_hlt_pt.GetBinLowEdge(j) for j in range(1, h_hlt_pt.GetNbinsX()+2)])

    h_hlt_ditau_pt = df_acc_hlt_matched_ditau.Histo1D(ROOT.RDF.TH1DModel("h", "", 20, 0, 100), f"gen_pion_pt{i}")
    c_hlt_ditau_pt = np.array([h_hlt_ditau_pt.GetBinContent(j) for j in range(1, h_hlt_ditau_pt.GetNbinsX()+1)])
    e_hlt_ditau_pt = np.array([h_hlt_ditau_pt.GetBinLowEdge(j) for j in range(1, h_hlt_ditau_pt.GetNbinsX()+2)])
    
    fig, ax = plt.subplots(figsize=(10, 10))
    hep.histplot(c_gen_pt, bins=e_gen_pt, label=f"Gen-level", ax=ax, histtype='step', linewidth=2)
    hep.histplot(c_gen_pt_filter, bins=e_gen_pt_filter, label=f"Gen-level $|\eta|<4$", ax=ax, histtype='step', linewidth=2)
    # hep.histplot(c_hlt_pt, bins=e_hlt_pt, label=f"HLT matched (no L1 filter) ", ax=ax, histtype='step', linewidth=2)
    hep.histplot(c_hlt_ditau_pt, bins=e_hlt_ditau_pt, label="Gen-level $p_T^{\pi_1, \pi_2}>35$ GeV  ", ax=ax, linewidth=2, histtype='fill', color='orange', alpha=0.5)
    ax.set_xlabel(f"{pion_name[i]} Pion $p_T$ [GeV]")
    ax.set_ylabel("# Pions")
    ax.legend(loc='upper right', fontsize=18)
    plt.grid()
    hep.cms.text("Work in Progress", ax=ax, fontsize=20)
    hep.cms.lumitext(text='Phase-II (14 TeV)', ax=ax, fontsize=20)
    print(f" ### INFO: Saving plot {output_dir}/gen_hlt_pion_pt_{i}.png")
    plt.savefig(f"{output_dir}/gen_hlt_pion_pt_{i}.png")
    plt.savefig(f"{output_dir}/gen_hlt_pion_pt_{i}.pdf")

    fig, ax = plt.subplots(figsize=(10, 10))
    e_hlt_pt_center = (e_hlt_pt[1:]+e_hlt_pt[:-1])/2
    e_hlt_pt_width = [e_hlt_pt[e+1]-e_hlt_pt[e] for e in range(len(e_hlt_pt[:-1]))]
    ax.errorbar(e_hlt_pt_center, c_hlt_pt/c_gen_pt_filter, 
        xerr=e_hlt_pt_width, yerr=GetEfficencyError(c_hlt_pt,c_gen_pt_filter), color='black')
    ax.set_ylabel('Reco Efficiency')
    ax.set_xlabel(f"{pion_name[i]} Pion {i} $\eta$")
    print(f" ### INFO: Saving plot {output_dir}/gen_hlt_pion_pt_{i}_eff.png")
    plt.savefig(f"{output_dir}/gen_hlt_pion_pt_{i}_eff.png")
    plt.savefig(f"{output_dir}/gen_hlt_pion_pt_{i}_eff.pdf")

#####################################################
### Plot gen pions eta

for i in range(3):
    h_gen_eta = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 30, -6, 6), f"gen_pion_eta{i}")
    c_gen_eta = np.array([h_gen_eta.GetBinContent(j) for j in range(1, h_gen_eta.GetNbinsX()+1)])
    e_gen_eta = np.array([h_gen_eta.GetBinLowEdge(j) for j in range(1, h_gen_eta.GetNbinsX()+2)])

    h_gen_eta_filter = df_acc.Histo1D(ROOT.RDF.TH1DModel("h", "", 30, -6, 6), f"gen_pion_eta{i}")
    c_gen_eta_filter = np.array([h_gen_eta_filter.GetBinContent(j) for j in range(1, h_gen_eta_filter.GetNbinsX()+1)])
    e_gen_eta_filter = np.array([h_gen_eta_filter.GetBinLowEdge(j) for j in range(1, h_gen_eta_filter.GetNbinsX()+2)])

    h_hlt_eta = df_acc_hlt_matched.Histo1D(ROOT.RDF.TH1DModel("h", "", 30, -6, 6), f"hlt_pion_eta{i}")
    c_hlt_eta = np.array([h_hlt_eta.GetBinContent(j) for j in range(1, h_hlt_eta.GetNbinsX()+1)])
    e_hlt_eta = np.array([h_hlt_eta.GetBinLowEdge(j) for j in range(1, h_hlt_eta.GetNbinsX()+2)])

    fig, ax = plt.subplots(figsize=(10, 10))
    hep.histplot(c_gen_eta, bins=e_gen_eta, label=f"Gen-level", ax=ax, histtype='step', linewidth=2)
    hep.histplot(c_gen_eta_filter, bins=e_gen_eta_filter, label=f"Gen-level $|\eta|<4$", ax=ax, histtype='step', linewidth=2)
    hep.histplot(c_hlt_eta, bins=e_hlt_eta, label=f"HLT-level", ax=ax, histtype='step', linewidth=2)
    ax.set_xlabel(f"{pion_name[i]} Pion {i} $\eta$")
    ax.set_ylabel("# Pions")
    ax.legend(loc='upper left', fontsize=16)
    plt.grid()
    print(f" ### INFO: Saving plot {output_dir}/gen_hlt_pion_eta_{i}.png")
    plt.savefig(f"{output_dir}/gen_hlt_pion_eta_{i}.png")
    plt.savefig(f"{output_dir}/gen_hlt_pion_eta_{i}.pdf")

    fig, ax = plt.subplots(figsize=(10, 10))
    e_hlt_eta_center = (e_hlt_eta[1:]+e_hlt_eta[:-1])/2
    e_hlt_eta_width = [e_hlt_eta[e+1]-e_hlt_eta[e] for e in range(len(e_hlt_eta[:-1]))]
    ax.errorbar(e_hlt_eta_center, c_hlt_eta/c_gen_eta_filter, 
        xerr=e_hlt_eta_width, yerr=GetEfficencyError(c_hlt_eta,c_gen_eta_filter), color='black')
    ax.set_ylabel('Reco Efficiency')
    ax.set_xlabel(f"{pion_name[i]} Pion {i} $\eta$")
    print(f" ### INFO: Saving plot {output_dir}/gen_hlt_pion_eta_{i}_eff.png")
    plt.savefig(f"{output_dir}/gen_hlt_pion_eta_{i}_eff.png")
    plt.savefig(f"{output_dir}/gen_hlt_pion_eta_{i}_eff.pdf")


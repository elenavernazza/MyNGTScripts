import ROOT

#####################################################
# Print summary of the event, including hard scatter particles and W boson decay products

ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;

    int event_summary(Vint GenPart_pdgId, Vint GenPart_status, Vint GenPart_genPartIdxMother, Vfloat GenPart_mass, Vfloat GenPart_pt, Vfloat GenPart_eta) {
        std::cout << std::endl;
        std::cout << "=== Event Summary ===" << std::endl;
        // W boson decay
        for (size_t i = 0; i < GenPart_pdgId.size(); ++i) {
            if (std::abs(GenPart_pdgId[i]) == 24) {
                std::cout << "W boson found (idx = " << i << ", pdgId = " << GenPart_pdgId[i] 
                    <<  ", GenPart_status = " << GenPart_status[i] 
                    << ", GenPart_mass = " << GenPart_mass[i] << ")" 
                    <<  ", GenPart_pt = " << GenPart_pt[i] 
                    << std::endl;
                std::cout << "  -> Decay products:" << std::endl;
                for (size_t j = 0; j < GenPart_pdgId.size(); ++j) {
                    if (GenPart_genPartIdxMother[j] == i) {
                        std::cout << "     idx " << j << ": pdgId = " << GenPart_pdgId[j] 
                        << ", GenPart_pt = " << GenPart_pt[j]
                        << ", GenPart_eta = " << GenPart_eta[j]
                        << " (mass = " << GenPart_mass[j]
                        << ", GenPart_status = " << GenPart_status[j] << ")" << std::endl;
                    }
                }
            }
        }
        std::cout << "=====================" << std::endl;
        return 1;
    }
""")

#####################################################
### Get 3 pions from W decay
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;
    ROOT::RVec<int> get_gen_pion_idices(Vfloat GenPart_pdgId, Vfloat GenPart_status, Vint GenPart_genPartIdxMother) {
        int i_gen_pi1 = -1;
        int i_gen_pi2 = -1;
        int i_gen_pi3 = -1;
        for (size_t i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
            if ((GenPart_genPartIdxMother.at(i_gen) >= 0) && 
                    (std::abs(GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen))) == 24) && 
                    GenPart_status.at(GenPart_genPartIdxMother.at(i_gen)) == 62) {
                if (std::abs(GenPart_pdgId.at(i_gen)) == 211) {
                    if (i_gen_pi1 == -1) {
                        i_gen_pi1 = i_gen;
                    } else if (i_gen_pi2 == -1) {
                        i_gen_pi2 = i_gen;
                    } else if (i_gen_pi3 == -1) {
                        i_gen_pi3 = i_gen;
                    } else {
                        std::cout << "More than three pions found, only first three will be used." << std::endl;
                        break;
                    }
                }
            }
        }
        return {i_gen_pi1, i_gen_pi2, i_gen_pi3};
    }
""")

ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;

    ROOT::RVec<int> get_sorted_gen_pion_idices(Vfloat GenPart_pdgId, Vfloat GenPart_status, Vint GenPart_genPartIdxMother, Vfloat GenPart_pt) {
        std::vector<std::pair<float, int>> pt_idx_pairs;
        for (size_t i_gen = 0; i_gen < GenPart_pdgId.size(); ++i_gen) {
            if ((GenPart_genPartIdxMother.at(i_gen) >= 0) && 
                (std::abs(GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen))) == 24) && 
                GenPart_status.at(GenPart_genPartIdxMother.at(i_gen)) == 62) {
                if (std::abs(GenPart_pdgId.at(i_gen)) == 211) {
                    pt_idx_pairs.emplace_back(GenPart_pt.at(i_gen), i_gen);
                }
            }
        }

        // Sort in descending pt
        std::sort(pt_idx_pairs.begin(), pt_idx_pairs.end(),
                  [](const std::pair<float, int>& a, const std::pair<float, int>& b) {
                      return a.first > b.first;
                  });

        // Return up to three indices
        ROOT::RVec<int> sorted_indices;
        for (size_t i = 0; i < std::min<size_t>(3, pt_idx_pairs.size()); ++i) {
            sorted_indices.push_back(pt_idx_pairs[i].second);
        }

        // If fewer than 3, fill with -1
        while (sorted_indices.size() < 3) {
            sorted_indices.push_back(-1);
        }

        return sorted_indices;
    }
""")

#####################################################
### Compute W mass from 3 pions
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;
    float get_W_mass (Vint pion_indices, Vfloat Part_pt, Vfloat Part_eta, Vfloat Part_phi, Vfloat Part_mass) {
        if (pion_indices[0] == -1 || pion_indices[1] == -1 || pion_indices[2] == -1) return 0;
        auto pion0_tlv = TLorentzVector();
        pion0_tlv.SetPtEtaPhiM(Part_pt.at(pion_indices[0]), Part_eta.at(pion_indices[0]), 
            Part_phi.at(pion_indices[0]), Part_mass.at(pion_indices[0]));
        auto pion1_tlv = TLorentzVector();
        pion1_tlv.SetPtEtaPhiM(Part_pt.at(pion_indices[1]), Part_eta.at(pion_indices[1]), 
            Part_phi.at(pion_indices[1]), Part_mass.at(pion_indices[1]));
        auto pion2_tlv = TLorentzVector();
        pion2_tlv.SetPtEtaPhiM(Part_pt.at(pion_indices[2]), Part_eta.at(pion_indices[2]), 
            Part_phi.at(pion_indices[2]), Part_mass.at(pion_indices[2]));
        auto W_tlv = pion0_tlv + pion1_tlv + pion2_tlv;
        return W_tlv.M();
    }
""")

#####################################################
### For each gen pion, among the hlt PF candidates with dR < 0.4 and dPt < 10 GeV, find the closest one in dR
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;

    ROOT::RVec<int> match_gen_to_hlt (Vint gen_pion_indices, Vfloat GenPart_pt, Vfloat GenPart_eta, Vfloat GenPart_phi, Vfloat GenPart_mass,
        Vfloat hltPFCandidate_pdgId, Vfloat hltPFCandidate_pt, Vfloat hltPFCandidate_eta, Vfloat hltPFCandidate_phi, Vfloat hltPFCandidate_mass) {

        const double dR_max = 0.4; // dR < dR_max
        ROOT::RVec<int> matched_hlt_pion_indices (gen_pion_indices.size(), -1);

        std::cout << "Start event" << std::endl;
        for (size_t i_gen_pion = 0; i_gen_pion < gen_pion_indices.size(); i_gen_pion ++) {
            if (gen_pion_indices.at(i_gen_pion) == -1) continue; // Skip if no gen pion found
            auto gen_tlv = TLorentzVector();
            gen_tlv.SetPtEtaPhiM(GenPart_pt.at(gen_pion_indices.at(i_gen_pion)), GenPart_eta.at(gen_pion_indices.at(i_gen_pion)), GenPart_phi.at(gen_pion_indices.at(i_gen_pion)), GenPart_mass.at(gen_pion_indices.at(i_gen_pion)));
            double matched_dR = std::numeric_limits<double>::infinity();
            int matched_i_hlt = -1;
            for (size_t i_hlt = 0; i_hlt < hltPFCandidate_pt.size(); i_hlt ++) {
                if (std::abs(hltPFCandidate_pdgId.at(i_hlt)) != 211) continue; // Only consider pions
                // if (hltPFCandidate_pt.at(i_hlt) < 1) continue; // Reduce combinatorics
                auto hlt_tlv = TLorentzVector();
                hlt_tlv.SetPtEtaPhiM(hltPFCandidate_pt.at(i_hlt), hltPFCandidate_eta.at(i_hlt), hltPFCandidate_phi.at(i_hlt), hltPFCandidate_mass.at(i_hlt));
                double dR = gen_tlv.DeltaR(hlt_tlv);
                if (dR > dR_max) continue;
                std::cout << "Gen Jet index " << gen_pion_indices.at(i_gen_pion) 
                    << ": HLT jet index " << i_hlt << " with dR = " << dR 
                    << " (Gen pt = " << GenPart_pt.at(gen_pion_indices.at(i_gen_pion))
                    << ", HLT pt = " << hltPFCandidate_pt.at(i_hlt) << ")" << std::endl;
                double dPt = std::abs(gen_tlv.Pt() - hlt_tlv.Pt());
                if (dPt > 10) continue;
                if (dR < matched_dR) {
                    matched_dR = dR;
                    matched_i_hlt = static_cast<int>(i_hlt);
                }
            }
            matched_hlt_pion_indices[i_gen_pion] = matched_i_hlt;
            std::cout << "Gen Jet index " << gen_pion_indices.at(i_gen_pion) << ": matched to hlt jet index " << matched_i_hlt << " with dR = " << matched_dR << std::endl;
        }
        return matched_hlt_pion_indices;
    }
""")

#####################################################
### For each gen pion, among the hlt PF candidates with dR < 0.4, find the closest in dPt
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;

    ROOT::RVec<int> match_pt_gen_to_hlt (Vint gen_pion_indices, Vfloat GenPart_pt, Vfloat GenPart_eta, Vfloat GenPart_phi, Vfloat GenPart_mass,
        Vfloat hltPFCandidate_pdgId, Vfloat hltPFCandidate_pt, Vfloat hltPFCandidate_eta, Vfloat hltPFCandidate_phi, Vfloat hltPFCandidate_mass,
        float dPt_max, float dR_max) {

        ROOT::RVec<int> matched_hlt_pion_indices (gen_pion_indices.size(), -1);

        // std::cout << "Start event" << std::endl;
        for (size_t i_gen_pion = 0; i_gen_pion < gen_pion_indices.size(); i_gen_pion ++) {
            if (gen_pion_indices.at(i_gen_pion) == -1) continue; // Skip if no gen pion found
            auto gen_tlv = TLorentzVector();
            gen_tlv.SetPtEtaPhiM(GenPart_pt.at(gen_pion_indices.at(i_gen_pion)), GenPart_eta.at(gen_pion_indices.at(i_gen_pion)), GenPart_phi.at(gen_pion_indices.at(i_gen_pion)), GenPart_mass.at(gen_pion_indices.at(i_gen_pion)));
            double matched_dPt = std::numeric_limits<double>::infinity();
            int matched_i_hlt = -1;
            for (size_t i_hlt = 0; i_hlt < hltPFCandidate_pt.size(); i_hlt ++) {
                if (std::abs(hltPFCandidate_pdgId.at(i_hlt)) != 211) continue; // Only consider pions
                auto hlt_tlv = TLorentzVector();
                hlt_tlv.SetPtEtaPhiM(hltPFCandidate_pt.at(i_hlt), hltPFCandidate_eta.at(i_hlt), hltPFCandidate_phi.at(i_hlt), hltPFCandidate_mass.at(i_hlt));
                double dR = gen_tlv.DeltaR(hlt_tlv);
                double dPt = std::abs(gen_tlv.Pt() - hlt_tlv.Pt());
                if (dR > dR_max) continue;
                if (dPt > dPt_max) continue;
                // std::cout << "Gen Jet index " << gen_pion_indices.at(i_gen_pion) 
                    // << " (pt = " << GenPart_pt.at(gen_pion_indices.at(i_gen_pion)) << ", eta = " << GenPart_eta.at(gen_pion_indices.at(i_gen_pion))
                    // << "): HLT jet index " << i_hlt << " with dR = " << dR 
                    // << " (pt = " << hltPFCandidate_pt.at(i_hlt) << ", eta = " << hltPFCandidate_eta.at(i_hlt) <<  ")" << std::endl;
                if (dPt < matched_dPt) {
                    matched_dPt = dPt;
                    matched_i_hlt = static_cast<int>(i_hlt);
                }
            }
            matched_hlt_pion_indices[i_gen_pion] = matched_i_hlt;
            // std::cout << "Gen Jet index " << gen_pion_indices.at(i_gen_pion) << ": matched to hlt jet index " << matched_i_hlt << " with dPt = " << matched_dPt << std::endl;
        }
        return matched_hlt_pion_indices;
    }
""")

#####################################################
### For each gen pion, among the hlt PF candidates with dR < 0.4, find the closest in dPt using the tracking information
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;

    ROOT::RVec<int> match_pt_track_gen_to_hlt (Vint gen_pion_indices, Vfloat GenPart_pt, Vfloat GenPart_eta, Vfloat GenPart_phi, Vfloat GenPart_mass,
        Vfloat hltPFCandidate_pdgId, Vfloat hltPFCandidate_pt, Vfloat hltPFCandidate_eta, Vfloat hltPFCandidate_phi, Vfloat hltPFCandidate_mass,
        Vint hltPFCandidate_trackIndex, Vfloat hltGeneralTrack_pt, Vfloat hltGeneralTrack_eta, Vfloat hltGeneralTrack_phi,
        float dPt_max, float dR_max) {

        ROOT::RVec<int> matched_hlt_track_indices (gen_pion_indices.size(), -1);

        std::cout << "Start event" << std::endl;
        for (size_t i_gen_pion = 0; i_gen_pion < gen_pion_indices.size(); i_gen_pion ++) {
            if (gen_pion_indices.at(i_gen_pion) == -1) continue; // Skip if no gen pion found
            auto gen_tlv = TLorentzVector();
            gen_tlv.SetPtEtaPhiM(GenPart_pt.at(gen_pion_indices.at(i_gen_pion)), GenPart_eta.at(gen_pion_indices.at(i_gen_pion)), GenPart_phi.at(gen_pion_indices.at(i_gen_pion)), GenPart_mass.at(gen_pion_indices.at(i_gen_pion)));
            double matched_dPt = std::numeric_limits<double>::infinity();
            int matched_i_hlt = -1;
            for (size_t i_hlt = 0; i_hlt < hltPFCandidate_pdgId.size(); i_hlt ++) {
                if (std::abs(hltPFCandidate_pdgId.at(i_hlt)) != 211) continue; // Only consider pions
                auto hlt_tlv = TLorentzVector();
                hlt_tlv.SetPtEtaPhiM(hltGeneralTrack_pt.at(hltPFCandidate_trackIndex.at(i_hlt)), hltGeneralTrack_eta.at(hltPFCandidate_trackIndex.at(i_hlt)), hltGeneralTrack_phi.at(hltPFCandidate_trackIndex.at(i_hlt)), 0);
                // hlt_tlv.SetPtEtaPhiM(hltPFCandidate_pt.at(i_hlt), hltPFCandidate_eta.at(i_hlt), hltPFCandidate_phi.at(i_hlt), hltPFCandidate_mass.at(i_hlt));
                double dR = gen_tlv.DeltaR(hlt_tlv);
                double dPt = std::abs(gen_tlv.Pt() - hltGeneralTrack_pt.at(hltPFCandidate_trackIndex.at(i_hlt)));
                if (dR > dR_max) continue;
                if (dPt > dPt_max) continue;
                std::cout << "Gen Part index " << gen_pion_indices.at(i_gen_pion) 
                    << " (pt = " << GenPart_pt.at(gen_pion_indices.at(i_gen_pion)) << ", eta = " << GenPart_eta.at(gen_pion_indices.at(i_gen_pion))
                    << "): HLT track index " << hltPFCandidate_trackIndex.at(i_hlt) << " with dR = " << dR 
                    << " (pt track = " << hltGeneralTrack_pt.at(hltPFCandidate_trackIndex.at(i_hlt)) << ", eta track = " << hltGeneralTrack_eta.at(hltPFCandidate_trackIndex.at(i_hlt))
                    << ", pt = " << hltPFCandidate_pt.at(i_hlt) << ", eta = " << hltPFCandidate_eta.at(i_hlt) << ", mass = " << hltPFCandidate_mass.at(i_hlt) <<  ")" << std::endl;
                if (dPt < matched_dPt) {
                    matched_dPt = dPt;
                    matched_i_hlt = static_cast<int>(hltPFCandidate_trackIndex.at(i_hlt));
                }
            }
            matched_hlt_track_indices[i_gen_pion] = matched_i_hlt;
            std::cout << "Gen Jet index " << gen_pion_indices.at(i_gen_pion) << ": matched to hlt track index " << matched_i_hlt << " with dPt = " << matched_dPt << std::endl;
        }
        return matched_hlt_track_indices;
    }
""")

# #####################################################
# ### Anmong the hlt PF candidates with dz < 0.001, get the 3 indeces with highest pt
# ROOT.gInterpreter.Declare("""
#     using Vfloat = const ROOT::RVec<float>&;
#     using Vint   = const ROOT::RVec<int>&;

#     ROOT::RVec<int> get_hlt_pions_dz (Vint gen_pion_indices, Vfloat GenPart_pt, Vfloat GenPart_eta, Vfloat GenPart_phi, Vfloat GenPart_mass,
#         Vfloat hltPFCandidate_pdgId, Vfloat hltPFCandidate_pt, Vfloat hltPFCandidate_eta, Vfloat hltPFCandidate_phi, Vfloat hltPFCandidate_mass) {

#         const double dR_max = 0.4; // dR < dR_max
#         ROOT::RVec<int> matched_hlt_pion_indices (gen_pion_indices.size(), -1);
                          
#         for (size_t i_hlt = 0; i_hlt < hltPFCandidate_pt.size(); i_hlt ++) {
#             auto hlt_tlv = TLorentzVector();
#             hlt_tlv.SetPtEtaPhiM(hltPFCandidate_pt.at(i_hlt), hltPFCandidate_eta.at(i_hlt), hltPFCandidate_phi.at(i_hlt), hltPFCandidate_mass.at(i_hlt));
#             for (size_t i_gen_pion = 0; i_gen_pion < gen_pion_indices.size(); i_gen_pion ++) {
#                 if (gen_pion_indices.at(i_gen_pion) == -1) continue; // Skip if no gen pion found
#                 auto gen_tlv = TLorentzVector();
#                 gen_tlv.SetPtEtaPhiM(GenPart_pt.at(gen_pion_indices.at(i_gen_pion)), GenPart_eta.at(gen_pion_indices.at(i_gen_pion)), GenPart_phi.at(gen_pion_indices.at(i_gen_pion)), GenPart_mass.at(gen_pion_indices.at(i_gen_pion)));
#                 double dR = hlt_tlv.DeltaR(gen_tlv);
#                 if (dR < 0.4) {
#                     std::cout << "HLT Jet index " << i_hlt << ": matched to gen jet index " << gen_pion_indices.at(i_gen_pion) << " with dR = " << dR << std::endl;
#                     gen_pion_indices.pop_back(); // Remove the matched gen pion index
#                 }
#             }
   
#     }
# """)

# #####################################################
# ### Find the best triplet of pions from HLT candidates using vertex

# ROOT.gInterpreter.Declare("""
#     using Vfloat = const ROOT::RVec<float>&;
#     using Vint   = const ROOT::RVec<int>&;
#     ROOT::RVec<int> get_best_hlt_pion_triplet(Vfloat hltPFCandidate_pdgId, Vfloat hltPFCandidate_pt,
#                                               Vfloat hltPFCandidate_eta, Vfloat hltPFCandidate_phi,
#                                               Vfloat hltPFCandidate_mass) {

#         std::vector<int> pion_indices;
#         for (size_t i = 0; i < hltPFCandidate_pdgId.size(); ++i) {
#             if (std::abs(hltPFCandidate_pdgId[i]) == 211) {
#                 pion_indices.push_back(i);
#             }
#         }

#         double target_mass = 80.4;  // W boson mass in GeV
#         double closest_mass_diff = std::numeric_limits<double>::max();
#         ROOT::RVec<int> best_triplet = {-1, -1, -1};

#         if (pion_indices.size() < 3)
#             return best_triplet;

#         std::cout << "Found " << pion_indices.size() << " pions in HLT candidates" << std::endl;
#         float min_pt = 10.0;
#         int count = 0;
#         for (size_t i = 0; i < pion_indices.size(); ++i) {
#             if (hltPFCandidate_pt[pion_indices[i]] < min_pt) continue;
#             count++;
#             for (size_t j = i + 1; j < pion_indices.size(); ++j) {
#                 if (hltPFCandidate_pt[pion_indices[j]] < min_pt) continue;
#                 for (size_t k = j + 1; k < pion_indices.size(); ++k) {
#                     if (hltPFCandidate_pt[pion_indices[j]] < min_pt) continue;

#                     int idx1 = pion_indices[i];
#                     int idx2 = pion_indices[j];
#                     int idx3 = pion_indices[k];

#                     TLorentzVector p1, p2, p3;
#                     p1.SetPtEtaPhiM(hltPFCandidate_pt[idx1], hltPFCandidate_eta[idx1],
#                                     hltPFCandidate_phi[idx1], hltPFCandidate_mass[idx1]);
#                     p2.SetPtEtaPhiM(hltPFCandidate_pt[idx2], hltPFCandidate_eta[idx2],
#                                     hltPFCandidate_phi[idx2], hltPFCandidate_mass[idx2]);
#                     p3.SetPtEtaPhiM(hltPFCandidate_pt[idx3], hltPFCandidate_eta[idx3],
#                                     hltPFCandidate_phi[idx3], hltPFCandidate_mass[idx3]);

#                     TLorentzVector total = p1 + p2 + p3;
#                     double mass_diff = std::abs(total.M() - target_mass);

#                     if (mass_diff < closest_mass_diff) {
#                         closest_mass_diff = mass_diff;
#                         best_triplet = {idx1, idx2, idx3};
#                     }
#                 }
#             }
#         }

#         std::cout << "Found " << count << " pions with pt > " << min_pt << std::endl;
#         return best_triplet;
#     }
# """)
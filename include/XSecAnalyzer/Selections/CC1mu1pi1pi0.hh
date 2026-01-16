#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class CC1mu1pi1pi0 : public SelectionBase {

public:

  CC1mu1pi1pi0();
  virtual int categorize_event( AnalysisEvent* Event ) override final;
  virtual bool selection( AnalysisEvent* Event ) override final;
  virtual bool define_signal( AnalysisEvent* Event ) override final;
  virtual void compute_reco_observables( AnalysisEvent* Event ) override final;
  virtual void compute_true_observables( AnalysisEvent* Event ) override final;
  virtual void define_output_branches() override final;
  virtual void define_constants() override final;
  virtual void define_category_map() override final;
  virtual void reset() override final;

private:

  bool sig_isNuMu_;
  bool sig_inFV_;
  bool sig_leadProtonMomInRange_;
  bool sig_muonInMomRange_;
  bool sig_noFSMesons_;
  bool sig_mc_no_fs_pi0_;
  bool sig_mc_no_charged_pi_above_threshold_;
  int sig_nProtons_in_Momentum_range_;
  int sig_nCharPions_in_Momentum_range_;
  int sig_nNeutPions_in_Momentum_range_;

  bool sel_reco_vertex_in_FV_;
  bool sel_pfp_starts_in_PCV_;
  bool sel_has_muon_candidate_;
  bool sel_topo_cut_passed_;
  bool sel_nu_mu_cc_;
  bool sel_muon_contained_;
  bool sel_muon_passed_mom_cuts_;
  bool sel_no_gen2_reco_showers_;
  bool sel_gen3_reco_showers_;
  bool sel_min_2_tracks_;
  bool sel_2_non_proton_;
  bool sel_has_pion_candidate_;
  bool sel_all_pfp_contained_;
  bool sel_all_pfp_in_vtx_proximity_;
  bool sel_pion_contained_;
  bool sel_pi0_contained_;
  bool sel_has_p_candidate_;
  bool sel_muon_quality_ok_;
  bool sel_protons_contained_;
  bool sel_passed_proton_pid_cut_;
  bool sel_lead_p_passed_mom_cuts_;
  bool sel_cosmic_ip_cut_passed_;
  bool sel_pion_passed_mom_cuts_;
  bool sel_tracks_flipped_;
  bool sel_has_pi0_candidate_;
  bool sel_pi0_passed_mom_cuts_;

  int lead_p_candidate_idx_;
  int muon_candidate_pid_idx_;
  int muon_candidate_length_idx_;
  int muon_candidate_bdt_idx_;
  int pion_candidate_idx_;
  int pi0_shr1_idx_;
  int pi0_shr2_idx_;

  int n_reco_tracks_;

  double delta_pT_;
  double delta_phiT_;
  double delta_alphaT_;
  double delta_pL_;
  double pn_;
  double delta_pTx_;
  double delta_pTy_;
  double theta_mu_p_;

  MyPointer< TVector3 > p3mu;
  MyPointer< TVector3 > p3p;
  MyPointer< TVector3 > p3cpi_;
  MyPointer< TVector3 > p3pi0_;
  MyPointer< TVector3 > nu_vtx_;
  MyPointer< TVector3 > shr1_start_;
  MyPointer< TVector3 > shr2_start_;
  MyPointer< TVector3 > dir1_;
  MyPointer< TVector3 > dir2_;
  MyPointer< TVector3 > p3_gamma1_;
  MyPointer< TVector3 > p3_gamma2_;
  MyPointer< std::vector< TVector3 > > p3_p_vec_;

  double mc_delta_pT_;
  double mc_delta_phiT_;
  double mc_delta_alphaT_;
  double mc_delta_pL_;
  double mc_pn_;
  double mc_delta_pTx_;
  double mc_delta_pTy_;
  double mc_theta_mu_p_;

  bool mc_golden_;
  bool mc_pi_stopping_;
  bool mc_pi0_stopping_;

  int mc_n_protons_;

  MyPointer< TVector3 > mc_p3mu_;
  MyPointer< TVector3 > mc_p3p_;
  MyPointer< TVector3 > mc_p3cpi_;
  MyPointer< TVector3 > mc_p3pi0_;
  MyPointer< std::vector< TVector3 > > mc_p3_p_vec_;
  MyPointer< std::vector< TVector3 > > mc_p3_cpi_vec_;
  MyPointer< std::vector< TVector3 > > mc_p3_pi0_vec_;

  STVCalcType calc_type;
};

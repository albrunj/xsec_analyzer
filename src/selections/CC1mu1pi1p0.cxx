// XSecAnalyzer includes
#include "XSecAnalyzer/TreeUtils.hh"
#include "XSecAnalyzer/Functions.hh"

#include "XSecAnalyzer/Selections/CC1mu1pi1pi0.hh"
#include "XSecAnalyzer/Selections/EventCategories1pi1pi0.hh"

CC1mu1pi1pi0::CC1mu1pi1pi0() : SelectionBase( "CC1mu1pi1pi0" ) {
  calc_type = kOpt1;
}

void CC1mu1pi1pi0::define_constants() {
  this->define_true_FV( 21.5, 234.85, -95.0, 95.0, 21.5, 966.8 );
  this->define_reco_FV( 21.5, 234.85, -95.0, 95.0, 21.5, 966.8 );
}

void CC1mu1pi1pi0::compute_true_observables( AnalysisEvent* Event ) {
  size_t num_mc_daughters = Event->mc_nu_daughter_pdg_->size();

  // Set the true 3-momentum of the final-state muon if there is one
  bool true_muon = ( sig_isNuMu_ && Event->mc_nu_ccnc_ == CHARGED_CURRENT );
  if ( true_muon ) {
    // Loop over the MC neutrino daughters, find the muon, and get its
    // true 3-momentum. Note that we assume there is only one muon in
    // this loop.
    bool found_muon = false;
    for ( size_t d = 0u; d < num_mc_daughters; ++d ) {
      int pdg = Event->mc_nu_daughter_pdg_->at( d );
      if ( pdg == MUON ) {
        found_muon = true;
        float px = Event->mc_nu_daughter_px_->at( d );
        float py = Event->mc_nu_daughter_py_->at( d );
        float pz = Event->mc_nu_daughter_pz_->at( d );
	*mc_p3mu = TVector3( px, py, pz );
        break;
      }
    }

    if ( !found_muon ) {
      std::cout << "WARNING: Missing muon in MC signal event!\n";
      return;
    }
  }

  // Set the true 3-momentum of the leading proton (if there is one)
  float max_mom = LOW_FLOAT;
  for ( int p = 0; p < num_mc_daughters; ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    if ( pdg == PROTON )
    {
      float px = Event->mc_nu_daughter_px_->at( p );
      float py = Event->mc_nu_daughter_py_->at( p );
      float pz = Event->mc_nu_daughter_pz_->at( p );
      TVector3 temp_p3 = TVector3( px, py, pz );

      mc_p3_p_vec_->push_back( temp_p3 );

      float mom = temp_p3.Mag();
      if ( mom > max_mom ) {
        max_mom = mom;
	*mc_p3p = temp_p3;
      }
    }
  }

  // TODO: reduce code duplication by just getting the leading proton
  // 3-momentum from this sorted vector
  // Sort the true proton 3-momenta in order from highest to lowest magnitude
  std::sort( mc_p3_p_vec_->begin(), mc_p3_p_vec_->end(), [](const TVector3& a,
    const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );

  // Set the true 3-momentum of the charged pion (if there is one)
  float max_mom_cpi = LOW_FLOAT;
  for (int p = 0; p < num_mc_daughters; ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    if ( std::abs(pdg) == PI_PLUS ) {

      mc_pi_stopping_ = Event->mc_end_p_->at(p) <= std::numeric_limits<float>::epsilon();
      int n_scatters = Event->mc_n_elastic_->at(p) + Event->mc_n_inelastic_->at(p);

      mc_golden_ = mc_pi_stopping_ && (n_scatters == 0);

      float px = Event->mc_nu_daughter_px_->at( p );
      float py = Event->mc_nu_daughter_py_->at( p );
      float pz = Event->mc_nu_daughter_pz_->at( p );

      TVector3 temp_p3 = TVector3( px, py, pz );

      // TODO: add to output branch
      mc_p3_cpi_vec_->push_back( temp_cpi_p3 );
      float mom = temp_cpi_p3.Mag();
      if (mom > max_mom_cpi) {
        max_mom_cpi = mom;
        *mc_p3cpi_ = temp_p3;
      }
    }
  }

  // Set the true 3-momentum of the neutral pion (if there is one)
  float max_mom_pi0 = LOW_FLOAT;
  for (int p = 0; p < num_mc_daughters; ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    if ( pdg == PI_ZERO ) {

      mc_pi0_stopping_ = Event->mc_end_p_->at(p) <= std::numeric_limits<float>::epsilon();
      int n_scatters = Event->mc_n_elastic_->at(p) + Event->mc_n_inelastic_->at(p);

      mc_golden_ = mc_pi_stopping_ && (n_scatters == 0);

      float px = Event->mc_nu_daughter_px_->at( p );
      float py = Event->mc_nu_daughter_py_->at( p );
      float pz = Event->mc_nu_daughter_pz_->at( p );

      TVector3 temp_p3 = TVector3( px, py, pz );

      // TODO: add to output branch
      mc_p3_pi0_vec_->push_back( temp_pi0_p3 );
      float mom = temp_pi0_p3.Mag();
      if (mom > max_mom_pi0) {
        max_mom_pi0 = mom;
        *mc_p3pi0_ = temp_pi0_p3;
      }
    }
  }


  // If the event contains a leading proton, then set the 3-momentum
  // accordingly
  bool true_lead_p = ( max_mom != LOW_FLOAT );
  if ( !true_lead_p && this->is_event_mc_signal() ) {
    // If it doesn't for a signal event, then something is wrong.
    std::cout << "WARNING: Missing leading proton in MC signal event!\n";
    return;
  }

  //if the event contains true charged pion, set the 3-momentum accordingly
  bool true_cpi = ( max_mom_cpi != LOW_FLOAT );

  //if the event contains true neutral pion, set the 3-momentum accordingly
  bool true_pi0 = ( max_mom_pi0 != LOW_FLOAT );


  // Compute true STVs if the event contains a muon, a charged pion, a neutral pion
  if ( true_muon && true_cpi && true_pi0 ) {
    double MuonEnergy = real_sqrt( mc_p3mu->Mag()*mc_p3mu->Mag()
      + MUON_MASS*MUON_MASS );
    //compute proton energy only if leading proton exists
    if (true_lead_p) {
      double ProtonEnergy = real_sqrt( mc_p3p->Mag()*mc_p3p->Mag()
        + PROTON_MASS*PROTON_MASS );
    }
    double ChargedPionEnergy = real_sqrt( mc_p3cpi->Mag()*mc_p3cpi->Mag()
      + PI_PLUS_MASS*PI_PLUS_MASS );
    double NeutralPionEnergy = real_sqrt( mc_p3pi0->Mag()*mc_p3pi0->Mag()
      + PI_ZERO_MASS*PI_ZERO_MASS );

    STVTools stv_tools;
    stv_tools.CalculateSTVs( *mc_p3mu, *mc_p3p, *mc_p3cpi, *mc_p3pi0, MuonEnergy, ProtonEnergy, ChargedPionEnergy, NeutralPionEnergy,
      calc_type );

    mc_delta_pT_ = stv_tools.ReturnPt();
    mc_delta_phiT_ = stv_tools.ReturnDeltaPhiT() * TMath::Pi()/180.;
    mc_delta_alphaT_ = stv_tools.ReturnDeltaAlphaT() * TMath::Pi()/180.;
    mc_delta_pL_ = stv_tools.ReturnPL();
    mc_pn_ = stv_tools.ReturnPn();
    mc_delta_pTx_ = stv_tools.ReturnPtx();
    mc_delta_pTy_ = stv_tools.ReturnPty();

    mc_theta_mu_p_ = std::acos( mc_p3mu->Dot(*mc_p3p)
      / mc_p3mu->Mag() / mc_p3p->Mag() );
  }
}

void CC1mu1pi1pi0::compute_reco_observables( AnalysisEvent* Event ) {

  // In cases where we failed to find a muon candidate, check whether there are
  // at least two generation == 2 PFParticles. If there are, then compute the
  // usual observables using the longest track as the muon candidate and the
  // second-longest track as the leading proton candidate. This will enable
  // sideband studies of NC backgrounds in the STV phase space.
  /*
  if ( !sel_has_muon_candidate_ ) {
    float max_trk_len = LOW_FLOAT;
    int max_trk_idx = BOGUS_INDEX;

    float next_to_max_trk_len = LOW_FLOAT;
    int next_to_max_trk_idx = BOGUS_INDEX;

    for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = Event->pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float trk_len = Event->track_length_->at( p );

      if ( trk_len > next_to_max_trk_len ) {
        next_to_max_trk_len = trk_len;
        next_to_max_trk_idx = p;

        if ( next_to_max_trk_len > max_trk_len ) {
          next_to_max_trk_len = max_trk_len;
          next_to_max_trk_idx = max_trk_idx;

          max_trk_len = trk_len;
          max_trk_idx = p;
        }
      }
    }

    // If we found at least two usable PFParticles, then assign the indices to
    // be used below
    if ( max_trk_idx != BOGUS_INDEX && next_to_max_trk_idx != BOGUS_INDEX ) {
      muon_candidate_idx_ = max_trk_idx;
      lead_p_candidate_idx_ = next_to_max_trk_idx;
    }
  }
  */

  // Set the reco 3-momentum of the muon candidate if we found one
  bool muon = muon_candidate_pid_idx_ != BOGUS_INDEX;
  if ( muon ) {
    float mu_dirx = Event->track_dirx_->at( muon_candidate_pid_idx_ );
    float mu_diry = Event->track_diry_->at( muon_candidate_pid_idx_ );
    float mu_dirz = Event->track_dirz_->at( muon_candidate_pid_idx_ );
    // The selection flag indicating whether the muon candidate is contained
    // was already set when the selection was applied. Use it to choose the
    // best momentum estimator to use.
    float muon_mom = LOW_FLOAT;
    if ( sel_muon_contained_ ) {
      muon_mom = Event->track_range_mom_mu_->at( muon_candidate_pid_idx_ );
    }
    else {
      muon_mom = Event->track_mcs_mom_mu_->at( muon_candidate_pid_idx_ );
    }

    *p3mu = TVector3( mu_dirx, mu_diry, mu_dirz );
    *p3mu = p3mu->Unit() * muon_mom;
  }

  // Set the reco 3-momentum of the charged pion candidate if we found one
  bool pion = pion_candidate_idx_ != BOGUS_INDEX;
  if ( pion ) {
    float pi_dirx = Event->track_dirx_->at( pion_candidate_idx_ );
    float pi_diry = Event->track_diry_->at( pion_candidate_idx_ );
    float pi_dirz = Event->track_dirz_->at( pion_candidate_idx_ );

    float pion_mom = LOW_FLOAT;
    float trk_length = Event->track_length_->at( pion_candidate_idx_ );
    pion_mom =  A + B*trk_length - C*std::pow(trk_length, -1.*D);

    *p3cpi_ = TVector3( pi_dirx, pi_diry, pi_dirz );
    *p3cpi_ = p3cpi_->Unit() * pion_mom;
  }

  // Set the reco 3-momentum of the leading proton candidate if we found one
  bool lead_p = lead_p_candidate_idx_ != BOGUS_INDEX;
  if ( lead_p ) {

    float p_dirx = Event->track_dirx_->at( lead_p_candidate_idx_ );
    float p_diry = Event->track_diry_->at( lead_p_candidate_idx_ );
    float p_dirz = Event->track_dirz_->at( lead_p_candidate_idx_ );
    float KEp = Event->track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
    float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

    *p3p = TVector3( p_dirx, p_diry, p_dirz );
    *p3p = p3p->Unit() * p_mom;
  }

  // set the reco 3-momenta of neutral pion candidates
  bool neutral_pion =  (pi0_shr1_idx_ != BOGUS_INDEX && pi0_shr2_idx_ != BOGUS_INDEX);
  if ( neutral_pion ) {
    
    //need neutrino vertex for direction calculation
    float vtx_x = Event->nu_vx_;
    float vtx_y = Event->nu_vy_;
    float vtx_z = Event->nu_vz_;
    MyPointer< TVector3 > nu_vtx = TVector3 (vtx_x, vtx_y, vtx_z);

    //photon 1
    float shr1_startx = Event->shower_startx_->at( pi0_shr1_idx_ );
    float shr1_starty = Event->shower_starty_->at( pi0_shr1_idx_ );
    float shr1_startz = Event->shower_startz_->at( pi0_shr1_idx_ );
    MyPointer< TVector3 > shr1_start = TVector3 (shr1_startx, shr1_starty, shr1_startz);
    MyPointer< TVector3 > dir1 = (*(shr1_start) - *(nu_vtx)).Unit();

    float shr1_energy = Event->shower_energy_->at( pi0_shr1_idx_ );
    MyPointer< TVector3 > p3_gamma1 = (*(dir1)) * shr1_energy;

    //photon 2
    float shr2_startx = Event->shower_startx_->at( pi0_shr2_idx_ );
    float shr2_starty = Event->shower_starty_->at( pi0_shr2_idx_ );
    float shr2_startz = Event->shower_startz_->at( pi0_shr2_idx_ );
    MyPointer< TVector3 > shr2_start = TVector3 (shr2_startx, shr2_starty, shr2_startz);
    MyPointer< TVector3 > dir2 = (*(shr2_start) - *(nu_vtx)).Unit();

    float shr2_energy = Event->shower_energy_->at( pi0_shr2_idx_ );
    MyPointer< TVector3 > p3_gamma2 = (*(dir2)) * shr2_energy;

    //neutral pion 3-momentum is sum of the two photon 3-momenta
    *p3pi0_ = *(p3_gamma1) + *(p3_gamma2);
    float pi0_mom = p3pi0_->Mag();
    float KE_pi0 = real_sqrt( pi0_mom*pi0_mom + PI_ZERO_MASS*PI_ZERO_MASS ) - PI_ZERO_MASS;

  }

  // Set the reco 3-momenta of all proton candidates (i.e., all generation == 2
  // tracks except the muon candidate) assuming we found both a muon candidate and charged pion and a neutral pion
  if ( muon && pion && neutral_pion) {
    for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
      // Skip the muon candidate
      if ( p == muon_candidate_pid_idx_ ) continue;

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = Event->pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float p_dirx = Event->track_dirx_->at( p );
      float p_diry = Event->track_diry_->at( p );
      float p_dirz = Event->track_dirz_->at( p );
      float KEp = Event->track_kinetic_energy_p_->at( p );
      float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

      TVector3 p3_temp( p_dirx, p_diry, p_dirz );
      p3_temp = p3_temp.Unit() * p_mom;

      p3_p_vec_->push_back( p3_temp );
    }

    // TODO: reduce code duplication by just getting the leading proton
    // 3-momentum from this sorted vector
    // Sort the reco proton 3-momenta in order from highest to lowest magnitude
    std::sort( p3_p_vec_->begin(), p3_p_vec_->end(), [](const TVector3& a,
      const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );
  }

  // Compute reco STVs if we have both a muon candidate
  // and a leading proton candidate in the event
  if ( muon && pion && neutral_pion) {
    double MuonEnergy = real_sqrt( p3mu->Mag()*p3mu->Mag()
      + MUON_MASS*MUON_MASS );
    if (lead_p){
      double ProtonEnergy	= real_sqrt( p3p->Mag()*p3p->Mag()
        + PROTON_MASS*PROTON_MASS );
    }

    double PionEnergy = real_sqrt( p3cpi_->Mag()*p3cpi_->Mag() + PI_PLUS_MASS*PI_PLUS_MASS );
    double Pi0Energy = real_sqrt( p3pi0_->Mag()*p3pi0_->Mag() + PI_ZERO_MASS*PI_ZERO_MASS );

    STVTools stv_tools;
    stv_tools.CalculateSTVs( *p3mu, *p3p, *p3cpi_, MuonEnergy, ProtonEnergy, PionEnergy, calc_type );

    delta_pT_ = stv_tools.ReturnPt();
    delta_phiT_ = stv_tools.ReturnDeltaPhiT() * TMath::Pi()/180.;
    delta_alphaT_ = stv_tools.ReturnDeltaAlphaT() * TMath::Pi()/180.;
    delta_pL_ = stv_tools.ReturnPL();
    pn_ = stv_tools.ReturnPn();
    delta_pTx_ = stv_tools.ReturnPtx();
    delta_pTy_ = stv_tools.ReturnPty();

    theta_mu_p_ = std::acos( p3mu->Dot(*p3p) / p3mu->Mag() / p3p->Mag() );
  }

}

bool CC1mu1pi1pi0::define_signal( AnalysisEvent* Event ) {

  sig_inFV_ = point_inside_FV( this->true_FV(), Event->mc_nu_vx_,
    Event->mc_nu_vy_, Event->mc_nu_vz_ );
  sig_isNuMu_ = ( Event->mc_nu_pdg_ == MUON_NEUTRINO );
  bool IsNC = ( Event->mc_nu_ccnc_ == NEUTRAL_CURRENT );

  sig_noFSMesons_= true;

  sig_muonInMomRange_ = false;

  sig_nProtons_in_Momentum_range = 0;
  sig_nCharPions_in_Momentum_range = 0;
  sig_nNeutPions_in_Momentum_range = 0;

  double LeadProtonMomentum = 0.;

  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    float energy = Event->mc_nu_daughter_energy_->at( p );

    // Do the general check for (anti)mesons first before considering
    // any individual PDG codes
    if ( is_meson_or_antimeson(pdg) ) {
      sig_noFSMesons_ = false;
    }

    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom >= MUON_P_MIN_MOM_CUT && mom <= MUON_P_MAX_MOM_CUT ) {
	      sig_muonInMomRange_= true;
      }
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > LeadProtonMomentum ) LeadProtonMomentum = mom;
      if ( mom >= LEAD_P_MIN_MOM_CUT && mom <= LEAD_P_MAX_MOM_CUT ) {
        sig_nProtons_in_Momentum_range++;
      }
    }
    else if ( pdg == PI_ZERO ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PI_ZERO_MASS, 2) );
      if ( mom > NEUTRAL_PI_MOM_CUT ) {
        sig_nNeutPions_in_Momentum_range++;
      }
    }
    else if ( std::abs(pdg) == PI_PLUS ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PI_PLUS_MASS, 2) );
      if ( mom > CHARGED_PI_MOM_CUT ) {
        sig_nCharPions_in_Momentum_range++;
      }
    }
  }

  sig_leadProtonMomInRange_ = false;
  // Check that the leading proton has a momentum within the allowed range
  if ( LeadProtonMomentum >= LEAD_P_MIN_MOM_CUT
    && LeadProtonMomentum <= LEAD_P_MAX_MOM_CUT )
  {
    sig_leadProtonMomInRange_ = true;
  }

  bool ReturnVal = sig_inFV_ && !IsNC && sig_isNuMu_ && sig_muonInMomRange_
    && sig_noFSMesons_ && (sig_nCharPions_in_Momentum_range_ == 1) && (sig_nNeutPions_in_Momentum_range_ == 1);
  return ReturnVal;
}

bool CC1mu1pi1pi0::selection( AnalysisEvent* Event ) {

  FiducialVolume PCV;
  PCV.X_Min = 10.;
  PCV.X_Max = 246.35;
  PCV.Y_Min = -106.5;
  PCV.Y_Max = 106.5;
  PCV.Z_Min = 10.;
  PCV.Z_Max = 1026.8;

  sel_reco_vertex_in_FV_ = point_inside_FV( this->reco_FV(),
    Event->nu_vx_, Event->nu_vy_, Event->nu_vz_ );

  sel_topo_cut_passed_ = Event->topological_score_ > TOPO_SCORE_CUT;
  sel_cosmic_ip_cut_passed_ = Event->cosmic_impact_parameter_ > COSMIC_IP_CUT;

  // Apply the containment cut to the starting positions of all
  // reconstructed tracks and showers. Pass this cut by default.
  sel_pfp_starts_in_PCV_ = true;

  // Loop over each PFParticle in the event
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Use the track reconstruction results to get the start point for
    // every PFParticle for the purpose of verifying containment. We could
    // in principle differentiate between tracks and showers here, but
    // (1) we cut out all showers later on in the selection anyway, and
    // (2) the blinded PeLEE data ntuples do not include shower information.
    // We therefore apply the track reconstruction here unconditionally.
    float x = Event->track_startx_->at( p );
    float y = Event->track_starty_->at( p );
    float z = Event->track_startz_->at( p );

    // Verify that the start of the PFParticle lies within the containment
    // volume.
    // TODO: revisit which containment volume to use for PFParticle start
    // positions. See https://stackoverflow.com/a/2488507 for an explanation
    // of the use of &= here. Don't worry, it's type-safe since both operands
    // are bool.
    sel_pfp_starts_in_PCV_ &= point_inside_FV(PCV, x, y, z );
  }

  // Sets the sel_has_muon_candidate_ flag as appropriate. The threshold check
  // is handled later.
  std::vector<int> muon_candidate_indices;
  std::vector<int> muon_pid_scores;
  std::vector<int> muon_lengths;
  std::vector<int> muon_bdt_scores;

  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
    // Only direct neutrino daughters (generation == 2) will be considered as
    // possible muon candidates
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    float track_score = Event->pfp_track_score_->at( p );
    float start_dist = Event->track_start_distance_->at( p );
    float track_length = Event->track_length_->at( p );
    float pid_score = Event->track_llr_pid_score_->at( p );
    float bdt_score = Event->muon_BDT_score_->at( p );

    if ( track_score > MUON_TRACK_SCORE_CUT
      && start_dist < MUON_VTX_DISTANCE_CUT
      && track_length > MUON_LENGTH_CUT
      && pid_score > MUON_PID_CUT )
    {
      muon_candidate_indices.push_back( p );
      muon_pid_scores.push_back( pid_score );
      muon_lengths.push_back( track_length );
      muon_bdt_scores.push_back( bdt_score );	
    }
  }

  size_t num_candidates = muon_candidate_indices.size();
  if ( num_candidates > 0u ) sel_has_muon_candidate_ = true;

  if ( num_candidates == 1u ) {
    muon_candidate_pid_idx_ = muon_candidate_indices.front();
    muon_candidate_length_idx_ = muon_candidate_indices.front();
    muon_candidate_bdt_idx_ = muon_candidate_indices.front();
  }
  else if ( num_candidates > 1u ) {
    // In the case of multiple muon candidates, choose the one with the highest
    // PID score (most muon-like) as the one to use
    float highest_score = LOW_FLOAT;
    int chosen_index_pid = BOGUS_INDEX;
    float longest = LOW_FLOAT;
    int chosen_index_length = BOGUS_INDEX;
    float most_muon_like = BOGUS + 1. ;
    int chosen_index_bdt = BOGUS_INDEX;

    for ( size_t c = 0; c < num_candidates; ++c ) {
      float score = muon_pid_scores.at( c );
      float length = muon_lengths.at( c );
      float muon_bdt_score = muon_bdt_scores.at( c );

      if ( highest_score < score ) {
        highest_score = score;
        chosen_index_pid = muon_candidate_indices.at( c );
      }
    }
    if ( longest < length ) {
        longest = length;
        chosen_index_length = muon_candidate_indices.at( c );
      }
  
      if (most_muon_like > muon_bdt_score ) {
        most_muon_like = muon_bdt_score;
        chosen_index_bdt = muon_candidate_indices.at( c );
      }
    }
    muon_candidate_pid_idx_ = chosen_index_pid;
    muon_candidate_length_idx_ = chosen_index_length;
    muon_candidate_bdt_idx_ = chosen_index_bdt;
  }
  else {
    muon_candidate_pid_idx_ = BOGUS_INDEX;
    muon_candidate_length_idx_ = BOGUS_INDEX;
    muon_candidate_bdt_idx_ = BOGUS_INDEX;
  }

  sel_nu_mu_cc_ = sel_reco_vertex_in_FV_ && sel_pfp_starts_in_PCV_
    && sel_has_muon_candidate_ && sel_topo_cut_passed_;

  int reco_gen2_shower_count = 0;
  int reco_gen3_shower_count = 0;
  int reco_track_count = 0;
  int n_non_proton_like = 0;

  for ( int p = 0; p < Event->num_pf_particles_; ++p ) { 
    
    unsigned int generation = Event->pfp_generation_->at( p );

    if (generation == 3u){
      float tscore = Event->pfp_track_score_->at( p );
      if ( tscore <= TRACK_SCORE_CUT ) ++reco_gen3_shower_count;
    }

    // Only check direct neutrino daughters (generation == 2)
    if ( generation != 2u ) continue;

    float proton_score = Event->proton_BDT_score_->at ( p );
    if ( proton_score < PROTON_BDT_CUT) ++n_non_proton_like;

    float tscore = Event->pfp_track_score_->at( p );
    if ( tscore <= TRACK_SCORE_CUT ) ++reco_gen2_shower_count;
    else ++reco_track_count;
  }
  // Check the electron shower cut
  sel_no_gen2_reco_showers_ = ( reco_gen2_shower_count == 0 );

  //check shower from possible pi0; need 1 or 2 showers (allows for one not being reconstructed properly)
  sel_gen3_reco_showers_ = ( reco_gen3_shower_count == 1 ) || ( reco_gen3_shower_count == 2 );

  // Check the track count. Need at least 2 for CC1piNp
  sel_min_2_tracks_ = ( reco_track_count >=2 );
  n_reco_tracks_ = reco_track_count;
  
  // Check we have 2 non proton like daugthers
  sel_2_non_proton_ = (n_non_proton_like == 2);
  int n_non_proton_like_ = n_non_proton_like;

  // If we have at least 2 non-proton like tracks, and 1 or 2 showers from possible pi0 find pion candidate
  // Only bother to do it if these conditions are left to save computational time
  pion_candidate_idx_ = BOGUS_INDEX;
  pi0_shr1_idx_ = BOGUS_INDEX;
  pi0_shr2_idx_ = BOGUS_INDEX;
  sel_has_pion_candidate_ = false;
  sel_has_pi0_candidate_ = false;

  if ( sel_no_gen2_reco_showers_ && sel_min_2_tracks_ && sel_2_non_proton_ && sel_gen3_reco_showers_){
    float current_proton_bdt_score = 1.1; 
    int current_pion_candidate_idx_ = BOGUS_INDEX;

    for (int p = 0; p < Event->num_pf_particles_; ++p) {

      unsigned int generation = Event->pfp_generation_->at( p );

      // Only check direct neutrino daughters (generation == 2)
      if ( generation != 2u ) continue;

      // Skip particles already identified as muons
      if ( p == muon_candidate_pid_idx_) continue;

      // Skip particles with bogus track score 
      float track_length = Event->track_length_->at( p );
      if (track_length <= 0. ) continue;

      float proton_bdt_score = Event->proton_BDT_score_->at ( p );

      // Skip particles for which BDT score could not be calculated
      // Use BOGUS - 1 to avoid compring floating point numbers
      if ( proton_bdt_score > BOGUS - 1. ) continue; 


      if ( proton_bdt_score < current_proton_bdt_score) {
        sel_has_pion_candidate_ = true; 
        current_pion_candidate_idx_ = p;
        current_proton_bdt_score = proton_bdt_score;       
      }
    }

    // set pion candidate index
    if (current_pion_candidate_idx_ != BOGUS_INDEX) pion_candidate_idx_ = current_pion_candidate_idx_;

  }

  // Now find pi0 candidate showers if we have a charged pion candidate to optimize compuation time
  if (sel_has_pion_candidate_){

    float pi0_shr1_cand_ = Event->shr_id_;
    float pi0_shr2_cand_ = Event->shr2_id_;

    //check shower score of leading shower and that we have at least 2 showers
    if (Event->shr_score_ < NEUTRAL_PI_SHOWER_SCORE_CUT && pi0_shr1_cand_>0 && pi0_shr2_cand_>0){
      if (Event->pfp_generation_->at( pi0_shr1_cand_ ) == 3u && Event->pfp_generation_->at( pi0_shr2_cand_ ) == 3u){
        //energy cuts on pi0 showers
        if (Event->shr_energy_cali_ > NEUTRAL_PI_SHR1_ENERGY_CUT && Event->shr_energy_second_cali_ > NEUTRAL_PI_SHR2_ENERGY_CUT){
          //dE/dx cut on primary shower
          if (Event->shr_tkfit_2cm_dedx_Y_ > NEUTRAL_PI_SHR1_DEDX_CUT){
            sel_has_pi0_candidate_ = true;
            pi0_shr1_idx_ = pi0_shr1_cand_;
            pi0_shr2_idx_ = pi0_shr2_cand_;
          }
        }
      }
    }
  }
  
  


  // If we don't have a neutral pion candidate, we can't have a CC1pi1pi0 event
  if ( !sel_has_pi0_candidate_)
  {
    bool sel_CC1pi1pi0_ = false;
    return sel_CC1pi1pi0_;
  }  

  // AT THIS POINT WE HAVE A MUON AND A PION CANDIDATE and a pi0.
  // Here we check if pfp emerge in vertex proximity,
  // and if all particles are contained 

  // Set flags that default to true here
  sel_passed_proton_pid_cut_ = true;
  //sel_protons_contained_ = true;

  // Set flags that default to false here
  sel_muon_contained_ = false;
  sel_pion_contained_ = false;
  sel_pi0_contained_ = false;
  sel_protons_contained_ = false;

  sel_all_pfp_contained_ = true;
  sel_all_pfp_in_vtx_proximity_ = true;
  
  sel_has_p_candidate_ = false;

  sel_muon_passed_mom_cuts_ = false;
  sel_pion_passed_mom_cuts_ = false;
  sel_pi0_passed_mom_cuts_ = false;
  sel_tracks_flipped_ = false;

  // We will also find the lead proton candidate here, if existing
  float lead_p_track_length = LOW_FLOAT;
  size_t lead_p_index = 0u;


  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
    // Only worry about direct neutrino daughters (PFParticles considered
    // daughters of the reconstructed neutrino)
    unsigned int generation = Event->pfp_generation_->at( p );

    //pi0 candidate check
    if (generation == 3u){

    }

    if ( generation != 2u ) continue;

    float start_dist = Event->track_start_distance_->at( p );
    if ( start_dist > PFP_DISTANCE_CUT ) sel_all_pfp_in_vtx_proximity_ = false;

    float startx = Event->track_startx_->at( p );
    float starty = Event->track_starty_->at( p );
    float startz = Event->track_startz_->at( p );

    float endx = Event->track_endx_->at( p );
    float endy = Event->track_endy_->at( p );
    float endz = Event->track_endz_->at( p );

    float reco_vtx_x = Event->nu_vx_;
    float reco_vtx_y = Event->nu_vy_;
    float reco_vtx_z = Event->nu_vz_;

    // Check if track is flipped. Flipped track start distance to vtx is higher than end distance to vertex
    float start_dist_to_vtx = std::sqrt(std::pow(startx - reco_vtx_x, 2) + std::pow(starty - reco_vtx_y, 2) + std::pow(startz - reco_vtx_z, 2));
    float end_dist_to_vtx = std::sqrt(std::pow(endx - reco_vtx_x, 2) + std::pow(endy - reco_vtx_y, 2) + std::pow(endz - reco_vtx_z, 2));

    if (start_dist_to_vtx > end_dist_to_vtx) {
      sel_tracks_flipped_ = true;
    }

    bool end_contained = point_inside_FV(PCV, endx, endy, endz );

    if( !(end_contained) )  sel_all_pfp_contained_ = false;

    // Check that we can find a muon candidate in the event. If more than
    // one is found, also fail the cut.
    if ( p == muon_candidate_pid_idx_ ) {

      if ( end_contained ) sel_muon_contained_ = true;
      else sel_all_pfp_contained_ = false;

      // Check that the muon candidate is above threshold. Use the best
      // momentum based on whether it was contained or not.

      float muon_mom = LOW_FLOAT;
      float range_muon_mom = Event->track_range_mom_mu_->at( p );
      float mcs_muon_mom = Event->track_mcs_mom_mu_->at( p );

      if ( sel_muon_contained_ ) muon_mom = range_muon_mom;
      else muon_mom = mcs_muon_mom;

      if ( muon_mom >= MUON_P_MIN_MOM_CUT && muon_mom <= MUON_P_MAX_MOM_CUT ) {
        sel_muon_passed_mom_cuts_ = true;
      }

      // Apply muon candidate quality cut by comparing MCS and range-based
      // momentum estimators. Default to failing the cut.
      sel_muon_quality_ok_ = false;

      double frac_diff_range_mcs = std::abs( range_muon_mom - mcs_muon_mom );
      if ( range_muon_mom > 0. ) {
        frac_diff_range_mcs /= range_muon_mom;
        if ( frac_diff_range_mcs < MUON_MOM_QUALITY_CUT ) {
          sel_muon_quality_ok_ = true;
        }
      }
    }
    // check if pion contained and passed momentum cuts
    else if (p == pion_candidate_idx_){

      if ( end_contained ) sel_pion_contained_ = true;
      else sel_all_pfp_contained_ = false;

      // Check that the pion candidate is above threshold
      float pion_mom = LOW_FLOAT;

      // pion momentum calibration
      float trk_length = Event->track_length_->at( pion_candidate_idx_ ); 
      pion_mom =  A + B*trk_length - C*std::pow(trk_length, -1.*D); 
      if ( pion_mom >= CHARGED_PI_MOM_CUT) {
        sel_pion_passed_mom_cuts_ = true;
      }

      // TODO: set pion momentum 
    }
    else {

      float track_score = Event->pfp_track_score_->at( p );
      if ( track_score <= TRACK_SCORE_CUT ) continue;

      // Bad tracks in the searchingfornues TTree can have
      // bogus track lengths. This skips those.
      float track_length = Event->track_length_->at( p );
      if ( track_length <= 0. ) continue;

      // Double check that the track is a proton candidate
      if ( p == muon_candidate_pid_idx_  || p == pion_candidate_idx_) continue;

      // We found a reco track that is not the muon candidate. All such
      // tracks are considered proton candidates.
      sel_has_p_candidate_ = true;

      float llr_pid_score = Event->track_llr_pid_score_->at( p );

      // Check whether the current proton candidate fails the proton PID cut
      if ( llr_pid_score > proton_pid_cut(track_length) ) {
        sel_passed_proton_pid_cut_ = false;
      }

      // Check whether the current proton candidate fails the containment cut
      if ( end_contained ) sel_protons_contained_ = true;
      else sel_all_pfp_contained_ = false;

      if ( track_length > lead_p_track_length ) {
        lead_p_track_length = track_length;
        lead_p_index = p;
    }

  }

  // Don't bother to apply the cuts that involve the leading
  // proton candidate if we don't have one
  if ( !sel_has_p_candidate_ ) {
    return false;
  }

  // All that remains is to apply the leading proton candidate cuts. We could
  // search for it above, but doing it here makes the code more readable (with
  // likely negligible impact on performance)
  float lead_p_track_length = LOW_FLOAT;
  size_t lead_p_index = 0u;
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Skip the muon candidate reco track (this function assumes that it has
    // already been found)
    if ( p == muon_candidate_pid_idx_ ) continue;

    // Skip PFParticles that are shower-like (track scores near 0)
    float track_score = Event->pfp_track_score_->at( p );
    if ( track_score <= TRACK_SCORE_CUT ) continue;

    // All non-muon-candidate reco tracks are considered proton candidates
    float track_length = Event->track_length_->at( p );
    if ( track_length <= 0. ) continue;

    if ( track_length > lead_p_track_length ) {
      lead_p_track_length = track_length;
      lead_p_index = p;
    }
  }

  // If the leading proton track length changed from its initial
  // value, then we found one. Set the index appropriately.
  if ( lead_p_track_length != LOW_FLOAT ) lead_p_candidate_idx_ = lead_p_index;
  // Otherwise, set the index to BOGUS_INDEX
  else lead_p_candidate_idx_ = BOGUS_INDEX;

  // Check the range-based reco momentum for the leading proton candidate
  float lead_p_KE = Event->track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
  float range_mom_lead_p = real_sqrt( lead_p_KE*lead_p_KE
    + 2.*PROTON_MASS*lead_p_KE );
  if ( range_mom_lead_p >= LEAD_P_MIN_MOM_CUT
    && range_mom_lead_p <= LEAD_P_MAX_MOM_CUT )
  {
    sel_lead_p_passed_mom_cuts_ = true;
  }

  // All right, we've applied all selection cuts. Set the flag that indicates
  // whether all were passed (and thus the event is selected as a CC1pi1pi0
  // candidate)
  bool sel_CC1pi1pi0_ = sel_nu_mu_cc_ && sel_no_gen2_reco_showers_ && sel_gen3_reco_showers_
    && sel_muon_passed_mom_cuts_ && sel_muon_contained_ && sel_muon_quality_ok_
    && sel_min_2_tracks_ && sel_2_non_proton_ && sel_has_pion_candidate_;
    //&& sel_has_p_candidate_ && sel_all_pfp_contained_&& sel_all_pfp_in_vtx_proximity_

  return sel_CC1pi1pi0_;
}

int CC1mu1pi1pi0::categorize_event(AnalysisEvent* Event) {
  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  int abs_mc_nu_pdg = std::abs( Event->mc_nu_pdg_ );
  Event->is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO
    || abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
  if ( !Event->is_mc_ ) {
    return kUnknown;
  }

  bool MCVertexInFV = point_inside_FV( this->true_FV(),
    Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_ );
  if ( !MCVertexInFV ) {
    return kOOFV;
  }

  bool isNC = ( Event->mc_nu_ccnc_ == NEUTRAL_CURRENT );
  // DB Currently only one NC category is supported so test first. Will likely
  // want to change this in the future
  if ( isNC ) return kNC;

  if ( Event->mc_nu_pdg_ == ELECTRON_NEUTRINO ) {
    return kNuECC;
  }
  if ( !(Event->mc_nu_pdg_ == MUON_NEUTRINO) ) {
    return kOther;
  }

  if ( this->is_event_mc_signal() ) {
   
      if ( Event->mc_nu_interaction_type_ == 0 ) {
        return kNuMuCC1pi1pi0_CCQE; // QE
      }
      else if ( Event->mc_nu_interaction_type_ == 10 ) {
        return kNuMuCC1pi1pi0_CCMEC; // MEC
      }
      else if ( Event->mc_nu_interaction_type_ == 1 ) {
        return kNuMuCC1pi1pi0_CCRES; // RES
      }
      else if (Event->mc_nu_interaction_type_ == 2) {
        return kNuMuCC1pi1pi0_CCDIS; // DIS
      }
      else if (Event->mc_nu_interaction_type_ == 3) {
        return kNuMuCC1pi1pi0_CCCOH; // COH
      }
      else return kNuMuCC1pi1pi0_CCOther;
  }
  return kNuMuCCOther;
}

void CC1mu1pi1pi0::define_output_branches() {

  set_branch( &sig_isNuMu_, "mc_is_numu" );
  set_branch( &sig_inFV_, "mc_vertex_in_FV" );
  set_branch( &sig_leadProtonMomInRange_, "mc_lead_p_in_range" );
  set_branch( &sig_muonInMomRange_, "mc_muon_in_mom_range" );
  set_branch( &sig_noFSMesons_, "mc_no_FS_mesons" );

  set_branch( &sig_mc_no_fs_pi0_, "mc_no_pi0s" );
  set_branch( &sig_nProtons_in_Momentum_range,
    "nProtons_in_Momentum_range" );
  set_branch( &sig_nCpi_in_Momentum_range_, "nCpi_in_Momentum_range" );

  set_branch( &sel_reco_vertex_in_FV_, "reco_vertex_in_FV" );
  set_branch( &sel_pfp_starts_in_PCV_, "pfp_starts_in_PCV" );
  set_branch( &sel_has_muon_candidate_, "has_muon_candidate" );
  set_branch( &sel_topo_cut_passed_, "topo_cut_passed" );
  set_branch( &sel_nu_mu_cc_, "nu_mu_cc" );
  set_branch( &sel_muon_contained_, "muon_contained" );
  set_branch( &sel_pion_contained_, "pion_contained" );
  set_branch( &sel_muon_passed_mom_cuts_, "muon_passed_mom_cuts" );
  set_branch( &sel_no_gen2_reco_showers_, "no_gen2_reco_showers" );
  set_branch( &sel_gen3_reco_showers_, "gen3_reco_showers" );
  set_branch( &sel_min_2_tracks_, "min_2_tracks" );
  set_branch( &sel_all_pfp_contained_, "all_pfp_contained" );
  set_branch( &sel_has_pion_candidate_, "has_pion_candidate" );
  set_branch( &sel_all_pfp_in_vtx_proximity_, "all_pfp_in_vtx_proximity" );
  set_branch( &sel_2_non_proton_, "2_non_proton" );
  set_branch( &sel_has_p_candidate_, "has_p_candidate" );
  set_branch( &sel_muon_quality_ok_, "muon_quality_ok" );
  set_branch( &sel_protons_contained_, "protons_contained" );
  set_branch( &sel_passed_proton_pid_cut_, "passed_proton_pid_cut" );
  set_branch( &sel_lead_p_passed_mom_cuts_, "lead_p_passed_mom_cuts" );
  set_branch( &sel_cosmic_ip_cut_passed_, "cosmic_ip_cut_passed" );
  set_branch( &sel_pion_passed_mom_cuts_, "pion_passed_mom_cuts" );
  set_branch( &sel_tracks_flipped_, "tracks_flipped" );

  set_branch( &n_reco_tracks_, "n_reco_tracks" );
  set_branch( &n_non_proton_like_, "n_non_proton_like" );

  set_branch( &lead_p_candidate_idx_, "lead_p_candidate_idx" );
  set_branch( &muon_candidate_pid_idx_, "muon_candidate_pid_idx" );
  set_branch( &muon_candidate_length_idx_, "muon_candidate_length_idx" );
  set_branch( &muon_candidate_bdt_idx_, "muon_candidate_bdt_idx" );
  set_branch( &pion_candidate_idx_, "pion_candidate_idx" );

  set_branch( &delta_pT_, "reco_delta_pT" );
  set_branch( &delta_phiT_, "reco_delta_phiT" );
  set_branch( &delta_alphaT_, "reco_delta_alphaT" );
  set_branch( &delta_pL_, "reco_delta_pL" );
  set_branch( &pn_, "reco_pn" );
  set_branch( &delta_pTx_, "reco_delta_pTx" );
  set_branch( &delta_pTy_, "reco_delta_pTy" );
  set_branch( &theta_mu_p_, "reco_theta_mu_p" );

  set_branch( p3mu, "reco_p3_mu" );
  set_branch( p3p, "reco_p3_lead_p" );
  set_branch( p3cpi_, "reco_p3_cpi" );
  set_branch( p3_p_vec_, "reco_p3_p_vec" );

  set_branch( &mc_delta_pT_, "true_delta_pT" );
  set_branch( &mc_delta_phiT_, "true_delta_phiT" );
  set_branch( &mc_delta_alphaT_, "true_delta_alphaT" );
  set_branch( &mc_delta_pL_, "true_delta_pL" );
  set_branch( &mc_pn_, "true_pn" );
  set_branch( &mc_delta_pTx_, "true_delta_pTx" );
  set_branch( &mc_delta_pTy_, "true_delta_pTy" );
  set_branch( &mc_theta_mu_p_, "true_theta_mu_p" );

  set_branch( mc_p3mu, "true_p3_mu" );
  set_branch( mc_p3p, "true_p3_lead_p" );
  set_branch( mc_p3cpi_, "true_p3_cpi" );
  set_branch( mc_p3_p_vec_, "true_p3_p_vec" );
  set_branch( mc_p3_cpi_vec_, "true_p3_cpi_vec" );
}

void CC1mu1pi1pi0::reset() {

  sig_isNuMu_ = false;
  sig_inFV_ = false;
  sig_leadProtonMomInRange_ = false;
  sig_muonInMomRange_ = false;
  sig_noFSMesons_ = false;
  sig_nProtons_in_Momentum_range = BOGUS_INDEX;
  sig_nCpi_in_Momentum_range_ = BOGUS_INDEX;

  sel_reco_vertex_in_FV_ = false;
  sel_pfp_starts_in_PCV_ = false;
  sel_has_muon_candidate_ = false;
  sel_topo_cut_passed_ = false;
  sel_nu_mu_cc_ = false;
  sel_no_gen2_reco_showers_ = false;
  sel_gen3_reco_showers_ = false;
  sel_min_2_tracks_ = false;
  sel_2_non_proton_ = false;
  sel_has_pion_candidate_ = false;
  sel_muon_contained_ = false;
  sel_pion_contained_ = false;
  sel_muon_passed_mom_cuts_ = false;
  sel_has_p_candidate_ = false;
  sel_muon_quality_ok_ = false;
  sel_all_pfp_contained_ = false;
  sel_all_pfp_in_vtx_proximity_ = false;
  sel_protons_contained_ = false;
  sel_passed_proton_pid_cut_ = false;
  sel_lead_p_passed_mom_cuts_ = false;
  sel_cosmic_ip_cut_passed_ = false;
  sel_pion_passed_mom_cuts_ = false;
  sel_tracks_flipped_ = false;
  sel_pi0_contained_ = false;
  sel_has_pi0_candidate_ = false;
  sel_pi0_passed_mom_cuts_ = false;
  

  n_reco_tracks_ = BOGUS_INDEX;
  n_non_proton_like_ = BOGUS_INDEX;

  lead_p_candidate_idx_ = BOGUS_INDEX;
  muon_candidate_pid_idx_ = BOGUS_INDEX;
  muon_candidate_length_idx_ = BOGUS_INDEX;
  muon_candidate_bdt_idx_ = BOGUS_INDEX;
  pion_candidate_idx_ = BOGUS_INDEX;
  neutral_pi0_candidate_idx_ = BOGUS_INDEX;
  pi0_shr1_idx_ = BOGUS_INDEX;
  pi0_shr2_idx_ = BOGUS_INDEX;

  delta_pT_ = BOGUS;
  delta_phiT_ = BOGUS;
  delta_alphaT_ = BOGUS;
  delta_pL_ = BOGUS;
  pn_ = BOGUS;
  delta_pTx_ = BOGUS;
  delta_pTy_ = BOGUS;
  theta_mu_p_ = BOGUS;

  *p3mu = TVector3();
  *p3p = TVector3();
  *p3cpi_ = TVector3();
  *p3cpi_ = TVector3();
  *p3pi0_ = TVector3();

  p3_p_vec_->clear();

  mc_delta_pT_ = BOGUS;
  mc_delta_phiT_ = BOGUS;
  mc_delta_alphaT_ = BOGUS;
  mc_delta_pL_ = BOGUS;
  mc_pn_ = BOGUS;
  mc_delta_pTx_ = BOGUS;
  mc_delta_pTy_ = BOGUS;
  mc_theta_mu_p_ = BOGUS;

  mc_pi_stopping_ = false;

  mc_n_protons_ = BOGUS_INDEX;

  *mc_p3mu_ = TVector3();
  *mc_p3p_ = TVector3();
  *mc_p3cpi_ = TVector3();
  *mc_p3pi0_ = TVector3();
  *mc_p3cpi_ = TVector3();

  mc_p3_p_vec_->clear();
  mc_p3_cpi_vec_->clear();

}

void CC1mu1pi1pi0::define_category_map() {
  // Use the shared category map for single charged pion and single neutral pion
  categ_map_ = CC1mu1pi1pi0_MAP;
}

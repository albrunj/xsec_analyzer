#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label event categories of interest for analysis plots in
// the CC1pi1pi0 analyses
enum EventCategoryXp {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = 0,

  // Signal events broken down by underlying reaction mode
  kNuMuCC1pi1pi0_CCQE = 1,
  kNuMuCC1pi1pi0_CCMEC = 2,
  kNuMuCC1pi1pi0_CCRES = 3,
  kNuMuCC1pi1pi0_CCDIS = 4,
  kNuMuCC1pi1pi0_CCCOH = 5,
  kNuMuCC1pi1pi0_Other = 6,

  // True numu CC event with no charged pions above threshold
  //kNuMuCC0pi = 7,

  //True numu CC event with 0 neutral pion above threshold
  //kNuMuCC0pi0 = 8,

  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  kNuMuCCOther = 7,

  // True nue CC event
  kNuECC = 8,

  // True neutral current event for any neutrino flavor
  kNC = 9,

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 10,

  // All events that do not fall within any of the other categories (e.g.,
  // numubar CC)
  kOther = 11,
};

extern std::map< int, std::pair< std::string, int > > CC1mu1pi1pi0_MAP;

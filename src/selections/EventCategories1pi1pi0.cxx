// ROOT includes
#include "TH1.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/EventCategories1pi1pi0.hh"

std::map< int, std::pair< std::string, int > > CC1mu1pi1pi0_MAP = {
  { kUnknown, { "Unknown", kGray } },
  { kNuMuCC0p0pi_CCQE, { "CCmu1pi1pi0 (CCQE)", kBlue - 2 } },
  { kNuMuCC0p0pi_CCMEC, { "CCmu1pi1pi0 (CCMEC)", kBlue - 6 } },
  { kNuMuCC0p0pi_CCRES, { "CCmu1pi1pi0 (CCRES)", kBlue - 9 } },
  { kNuMuCC0p0pi_CCDIS, { "CCmu1pi1pi0 (CCDIS)", kBlue - 4 } },
  { kNuMuCC0p0pi_CCCOH, { "CCmu1pi1pi0 (CCCOH)", kBlue - 8 } },
  { kNuMuCC0p0pi_Other, { "CCmu1pi1pi0 (Other)", kBlue - 10 } },
  //{ kNuMuCC0pi, { "#nu_{#mu} CC 0#pi", kGreen + 2 } },
  //{ kNuMuCC0pi0, { "#nu_{#mu} CC 0#pi^{0}", kGreen + 4 } },
  { kNuMuCCOther, { "Other #nu_{#mu} CC", kAzure } },
  { kNuECC, { "#nu_{e} CC", kViolet } },
  { kNC, { "NC", kOrange } },
  { kOOFV, {"Out FV", kRed + 3 } },
  { kOther, { "Other", kRed + 1 } }
};

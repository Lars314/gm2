////////////////////////////////////////////////////////////////////////
// Class:       testblinder
// Plugin Type: analyzer (art v2_09_02)
// File:        testblinder_module.cc
//
// Generated at Mon Mar 19 01:17:36 2018 by vagrant using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "gm2util/blinders/Blinders.hh"

class testblinder;

using namespace blinding;

class testblinder : public art::EDAnalyzer {
public:
  explicit testblinder(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  testblinder(testblinder const &) = delete;
  testblinder(testblinder &&) = delete;
  testblinder & operator = (testblinder const &) = delete;
  testblinder & operator = (testblinder &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.

};


testblinder::testblinder(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void testblinder::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  std::cout << "testblinder::analyze() ..." << std::endl;
  
  Blinders::fitType ftype = Blinders::kOmega_a;
  Blinders myBlinder( ftype );
  
  Blinders getBlinded( ftype, "Chris eats at Two Brothers" );
  
  Blinders systematicallyBlinded( ftype, 1, 10, "Help me! I think I'm falling..." );
  
  // there should be no blinding...
  std::cout << "\n\n should be unblinded results " << std::endl;
  for ( double R = 0; R < 10; R += 1 ) {
    double result = ( myBlinder.paramToFreq( R ) / myBlinder.referenceValue() ) - 1;
    std::cout << " input R: " << R << "   output: " << result << std::endl;
  }

  // there should be nominal blinding...
  std::cout << "\n\n should be blinded central results " << std::endl;
  for ( double R = 0; R < 10; R += 1 ) {
    double result = ( getBlinded.paramToFreq( R ) / getBlinded.referenceValue() ) - 1;
    std::cout << " input R: " << R << "   output: " << result << std::endl;
  }

  // there should be systematic shifts relative to nominal central...
  std::cout << "\n\n should be systematic shift results " << std::endl;
  for ( double R = 0; R < 10; R += 1 ) {
    double result = ( systematicallyBlinded.paramToFreq( R ) / systematicallyBlinded.referenceValue() ) - 1;
    std::cout << " input R: " << R << "   output: " << result << std::endl;
  }

#ifdef DEBUG_BLINDING
  Blinders b1( ftype, "Lee prefers Rock Bottom" );
  Blinders b2( ftype, "Dave drinks fine wine" );
  Blinders b3( ftype, "Lawrence eats too much cheese" );
  Blinders b4( ftype, "Jarek and Brendan must work too hard" );
#endif
  

}

DEFINE_ART_MODULE(testblinder)

A random number library using the Mersenne Twister random number
generator.

Written by Charles Karney <charles@karney.com> and licensed under
the MIT/X11 License.  For more information, see

    http://randomlib.sourceforge.net/

Files

    00README.txt  -- this file
    AUTHORS -- the authors of the library
    LICENSE.txt -- the MIT/X11 License
    INSTALL -- brief installation instructions
    NEWS -- a history of changes

    include/RandomLib/ and src/
      Config.hh.in, Config.hh -- system dependent configuration
      Random.hh -- main include file plus implementation
      RandomCanonical.hh -- Random integers, reals, booleans
      RandomPower2.hh -- scaling by powers of two
      RandomEngine.hh -- abstract random number generator
      RandomAlgorithm.hh -- MT19937 and SFMT19937 random generators
      RandomMixer.hh -- mixing functions to convert seed to state
      RandomSeed.hh -- seed management
      RandomType.hh -- support of unsigned integer types
      NormalDistribution.hh -- sample from normal distribution
      ExponentialDistribution.hh -- sample from exponential distribution
      RandomSelect.hh -- sample from discrete distribution
      LeadingZeros.hh -- count of leading zeros on random fraction
      ExponentialProb.hh -- true with probability exp(-p)
      RandomNumber.hh -- support for infinite precision randoms
      ExactExponential.hh -- sample exactly from exponential distribution
      ExactNormal.hh -- sample exactly from normal distribution
      ExactPower.hh -- sample exactly from power distribution
      UniformInteger.hh -- sample partially from a integer range
      DiscreteNormal.hh -- sample exactly from discrete normal distribution
      DiscreteNormalAlt.hh -- alternative to DiscreteNormal.hh
      MPFR*.hh -- implementation of some distributions in MPFR

    src/
      Random.cpp -- code for implementation

    examples/
      RandomExample.cpp -- example code
      RandomTime.cpp -- time the random routines
      RandomSave.cpp -- different ways to save and restore the state
      RandomThread.cpp -- multi-threading example
      RandomLambda.cpp -- using the STL and lambda expressions
      RandomCoverage.cpp -- a code coverage test
      RandomExact.cpp -- generating exact samples are various distibutions
      MPFRExample.cpp -- using the algorithms in MPFR

    windows/
      RandomLib-vc9.sln -- MS Studio 2008 solution
      Random-vc9.vcproj -- project for library
      RandomExample-vc9.vcproj -- project for RandomExample
      RandomPermutation-vc9.vcproj -- project for RandomPermutation
      RandomTime-vc9.vcproj -- project for RandomTime
      RandomSave-vc9.vcproj -- project for RandomSave
      RandomThread-vc9.vcproj -- project for RandomThread
      also files for MS Studio 2005 (with vc8)

    doc/
      doxyfile.in -- Doxygen config file
      Random.dox -- main page of Doxygen documentation
      *.png *.pdf -- figures for documentation

    Makefile.mk -- Unix/Linux makefiles
    CMakeLists.txt -- cmake configuration files
    cmake/
      FindRandomLib.cmake -- cmake find script
      *.cmake.in -- cmake config templates

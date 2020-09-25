"""
------------------------------- Behold  SplIland -------------------------------

contains the spline and island classes for making artificial
pileup islands. This code is the best thing to have ever been
made by any human being ever and if you disagree you can go eat
a cookie because you are probably right.

Written by Lars Borchert with extensive help from Josh LaBounty

This dashed line is 80 characters long
--------------------------------------------------------------------------------
"""

import numpy as np
from scipy import stats



class Spline:
    """
    Take a root spline and make a simple,
    pythonified object with two np arrays.


    ------------ Attributes ------------
    
    splineY      : [ndarray] the y values in the spline. They are
                   normalized such that their integral in time
                   is equal to one
                   
    time         : [ndarray] the x values in the spline, which 
                   are time
    
    samplingRate : [float] the sampling rate of the undigitized
                   spline
    
    xtalNum      : [int] the crystal number of the template
    """
#-------------------------------------------------------------------------------
#------------------------------ Public  Functions ------------------------------
#-------------------------------------------------------------------------------


    def peakTime(self):
        """
        Gives the time when a peak occurs
        The idea is to have a time for the center of a pulse
        """
        index = np.argmax(self.splineY)
        return(self.time[index])


    def getPeak(self):
        """
        Gives the value of the peak of the spline
        """
        return np.max(self.splineY)
    
    
    def getNormalIntegral(self):
        """
        Gives the pulse integral of the peak-normalized spline
        Essentially, a unit pulse integral. peak-normalized in
        this case means that the peak of the spline is equal
        to one, rather than the integral of the spline equal to one
        """
        scaledVals = self.splineY / self.getPeak()
        return np.sum(scaledVals) * self.samplingRate


#-------------------------------------------------------------------------------
#------------------------------ Private Functions ------------------------------
#-------------------------------------------------------------------------------

    
    def __init__(self, rSpline, xtalNum):
        """
        rSpline : a ROOT spline, representing a normalized pulse
                  template
        """

        spline = []
        times = []

        for i in np.linspace(rSpline.GetXmin(),
                             rSpline.GetXmax(),
                             rSpline.GetNpx()):
            spline.append(rSpline.Eval(i))
            times.append(i)

        self.splineY = np.array(spline)
        self.time = np.array(times)

        self.samplingRate = times[2] - times[1]

        self.xtalNum = xtalNum



class Island:
    """
    ------------ Attributes ------------

    spline                  : [spline] the spline template from which the pulses
                              in this island are generated

    nPulses                 : [int] the number of pulses in this island
    
    energyCalibrationVal    : [float] the energy calibration constant of the
                              crystal from which the spline templates for this
                              island come from. This constant is used to convert
                              a pulse integral made over the island in ADC units
                              to an energy in MeV
                              
    minTimeOffset           : [float] the artificial minimum time offset between
                              pulses in this island in ns. Defaults to 0.0
                              
    gainSag                 : [None] short term double pulse correction, well
                              mapped function, fix for non-recovered pixels, if
                              two events occur very close in time and pixels
                              have not recovered, we might see a lower energy
                              value. This is a correction function for that. It
                              is not important unless we want to be able to
                              extract energies in our fitting, as it mainly just
                              changed the energy scale factors. So it is None
                              and not used, but kept as a reminder of this
                              
    normalize               : [boolean] if this island is normalized or not.
                              Defaults to true

    chopThreshold           : [float] the threshold value in ADC units for when
                              the chopping algorithm decides where a pulse is.
                              Defaults to 50
                              
    nPreSamples             : [int] the number of samples kept in chopping
                              before the first value above chopThreshold.
                              Defaults to 8
                              
    nPostSamples            : [int] the number of samples kept in chopping
                              after the last value above chopThreshold.
                              Defaults to 18
                              
    noise                   : [boolean] if this island has noise applied or
                              not. Defaults to true
                              
    noiseLevel              : [float] the standard deviation of gaussian
                              noise in the sample, in ADC units. Defaults
                              to 4
                              
    randomizeTailLength     : [boolean] if this island has a randomized tail
                              length or not. Defaults to false
                              
    nTailSamples            : [int] the number of tail samples in the
                              island. defaults to 150
                              
    samplingRateArtificial  : [float] the sampling rate of the island, which
                              is used in the digitization process. To match
                              the ADC, this defaults to 1.25ns and should
                              probably be kept there

    pedestal                : [float] the pedestal value of the ADC. Really,
                              it is the average pedestal value of the two
                              ADCs that interweave. But when the pedestal
                              offset correction is applied, we can think of
                              these two 400 MHz ADCs as one 800 MHz ADC

    integral                : [float] the value of the integral of this
                              island; the area under its curve
                              
    totalEnergy             : [float] the total energy, in MeV, of in this
                              island
                              
    timeOffsets             : [list] the time offsets of each pulse in this
                              island from the original template pulse
                              
    energyScaleFactors      : [list] the energy scaling factors of the
                              pulses in the island, in ADC units
                              
    time                    : [ndarray] the time values in this island

    yValues                 : [ndarray] the y values in this island. These
                              might be in ADC units, or they might be
                              normalized, depending on the value of
                              normalize
                              
    choppedTime             : [ndarray] the time values in this island
                              after chopping has been applied
                              
    choppedYValues          : [ndarray] the y values in this island after
                              the chopping has been applied. Just like
                              yValues, these may be in ADC units or they
                              may be normalized
    """
#-------------------------------------------------------------------------------
#------------------------------ Public  Functions ------------------------------
#-------------------------------------------------------------------------------

    def makeDF(self):
        """
        return a simple pandas dataframe
        with time and yValues columns for 
        this island
        """
        d = {'time': self.time,
             'yValues': self.yValues}
        df = pd.DataFrame(data = d)
        return(df)


#-------------------------------------------------------------------------------
#------------------------------ Private Functions ------------------------------
#-------------------------------------------------------------------------------


    def digitize_(self, verbosity, thisIsland, theseTimes):
        """
        convert the time sampling to match that of the ADC
        by default this should be 1.25ns sampling
        """
        
        if(self.samplingRateArtificial is None):
            return thisIsland, theseTimes
        
        
        assert type(self.samplingRateArtificial) is float, print(
            """Error: Sampling rate must be a float.
            Currently {0:.2f}""".format(type(self.samplingRateArtificial)))


        if(verbosity):
            print("Sampling this spline with a deltaT of " + \
                  "{0:.2f} ns".format(self.samplingRateArtificial) )
        
        sampledTimes = [theseTimes[0]]

        while(sampledTimes[-1] < [theseTimes[-1]]):
            sampledTimes.append(sampledTimes[-1] + self.samplingRateArtificial)
       
        theseSamples = np.interp(sampledTimes, theseTimes, thisIsland)
        
        thisIsland = np.array(theseSamples)
        theseTimes = np.array(sampledTimes)
        
        return thisIsland, theseTimes
    
    
    def giveNPulses_(self, verbosity, minPulses, maxPulses):
        """
        figure out how many pulses there will be
        never run this outside of initialization
        """
        nPulses = np.random.randint(minPulses, maxPulses+1)
        if(verbosity):
            print("The number of pulses in this island: " + \
                  "{0}".format(nPulses))
        return nPulses
            
            
    def giveRandomTailLength_(self, verbosity):
        """
        figure out how long the tail will be
        never run this outside of initialization
        """
        nTailSamples = self.nTailSamples + np.random.randint(-50, 100)
        if(verbosity): print(nTailSamples, " tail samples in this island")
        return nTailSamples


    def giveTimeOffsets_(self, verbosity, deltaTmin, deltaTmax):
        """
        Creates the time offset values for each pulse
        for the love of kale make deltaTmax >> minTimeOffset, otherwise
        this will take absolutely forever to run
        """

        
        check = 0
        # keep looking until every combination of offsets has passed the check
        while(check < self.nPulses ** 2):
            # make some random time offsets
            offsets = np.random.uniform(deltaTmin, deltaTmax, size=self.nPulses)
            # now we need to make sure they have the right separation
            for i in range(0, len(offsets)):
                for j in range(0, len(offsets)):
                    # a value can be within minTimeOffset of itself of course
                    if(i != j):
                        # if the separation isn't right, reset check
                        # otherwise we are ok, and can increment the number
                        # of correct pairs
                        separation = offsets[i] - offsets[j]
                        if(abs(separation) < self.minTimeOffset):
                            check = 0
                        else:
                            check += 1
                    else:
                        check += 1
                        
        if(verbosity):
            thisString = "["
            for offset in offsets:
                thisString += " {0:.2f} ".format(offset)
            thisString += '] ns'
            print("The time offsets are: " + thisString)
            print("The artificial minimum time offset is: " + \
                  "{0:.2f} ns".format(self.minTimeOffset))
        return offsets
            
       
    def giveUniformEnergyScaleFactors_(self, verbosity, minEscale, maxEscale):
        """
        Gives energy scale factors in ADC units along a linear distribution
        
        No energy is more likely than any other within the given bounds
        """
        energyScaleFactors = np.random.uniform(minEscale,
                                               maxEscale, size=self.nPulses)
        
        if(verbosity):
            thisString = "["
            for factor in factors:
                thisString += " {0:.2f} ".format(factor)
            thisString += '] ADC units'
            print("The energy scale factors are: " + thisString)
            
        return energyScaleFactors
            
            
    def giveEnergyScaleFactors_(self, verbosity, energyPeak,
                                energyScale, calibConstant):
        """
        Gives the energy scale factors according to the energy disctribution
        of the real world data. This is a gaussian distribution centered at
        zero. Only positive values are taken. Values below a certain
        threshold of around 50 MeV are ignored

        First an energy value in MeV is given. Because energy is directly
        proportional to the pulse integral, which is directly proportional
        to the pulse height, energy is directly proportional to pulse height.
        So we can convert our energy to a pulse height, by looking at the ratio
        between height and energy for some other pulse of the same shape. We
        choose this to be this template pulse, but of height one, which makes
        the calculation easy.
        
        An older, redundant method is as follows. I am only keeping it for
        when I inevitably confuse myself and need to re-read it to understand.
        
        
        170 MeV = calibration constant * pulse integral, by definition
        calibIntegral = 170 / calibConstant
        
        Energy is directly proportional to integral,
        is directly proportional to height
        calibHeight / calibIntegral = normalHeight / normalIntegral
        normalHeight == 1, by definition
        calibHeight = calibIntegral / normalIntegral
        
        Energy is directly proportional to integral,
        is directly proportional to height
        thisHeight / calibHeight = thisEenrgy / calibEnergy
        calibEnergy == 170 by definition
        thisHeight = (thisEnergy *calibHeight) / (170)
        
        calibHeight = calibIntegral / normalIntegral, from earlier
        thisHeight = (thisEnergy*calibIntegral/NormalIntegral) / (170)
        
        calibIntegral = 170 / calibConstant, from earlier
        thisHeight = (thisEnergy*(170 / calibConstant)/normalIntegral) / (170)
        thisHeight = thisEnergy / (calibConstant * normalIntegral)
        
        Energy is directly proportional to integral,
        is directly proportional to height
        thisHeight / normalHeight = thisEnergy / normalEnergy
        """
        
        factors = []
        normalIntegral = self.spline.getNormalIntegral()
        
        for pulse in range(0, self.nPulses):
            """
            generate an energy, in Mev
            we ignore stuff under a threshhold
            100 here ends up being around 50MeV, a little higher
            it is imperfect and needs work
            """
            thisEnergy = 0
            while(thisEnergy < 100):
                thisEnergy = abs(np.random.normal(loc=energyPeak,
                                                  scale=energyScale))

            thisHeight = thisEnergy / (calibConstant * normalIntegral)
            
            factors.append(thisHeight)
        
        if(verbosity):
            thisString = "["
            for factor in factors:
                thisString += " {0:.2f} ".format(factor)
            thisString += '] ADC units'
            print("The energy scale factors are: " + thisString)
        return factors
        
  
    def giveNoise_(self, verbosity, thisIsland):
        """
        adds gaussian noise into the island, in ADC units. Do not use this
        function unless you are operating in ADC units
        """
        thisIsland += np.random.normal(0, self.noiseLevel,
                                       size=thisIsland.size)
        if(verbosity):
            print("The noise level standard deviation in this island: " + \
                  "{0:.2f}".format(self.noiseLevel))
        return thisIsland
        

    def chop_(self, times, island, critValue, preSamples, postSamples):
        """
        Performs simulated island chopping. We don't want a bunch of zero
        values if we have many many terabytes of data, so we cut islands
        down to just the interesting stuff. We keep preSamples before the
        first value over our threshold critValue, and postSamples after the
        last
        """
        startIndex = None
        endIndex = None
        cuts = []
        
        # loop through all the values
        for index, val in enumerate(island):
            if(val > critValue):
                # for something over our threshold, make an endIndex
                endIndex = index + postSamples
                
                # make sure we don't get out of bounds errors
                if(endIndex >= len(island)):
                    cuts.append(len(island))
                else:
                    cuts.append(endIndex)
                
                # if this is the first threshold passer, get the startIndex
                if(startIndex is None):
                    if(index - preSamples > 0):
                        startIndex = index - preSamples
                    else:
                        startIndex = 0

        return times[startIndex:cuts[-1]], island[startIndex:cuts[-1]]
    

    def scaleADC_(self, verbosity, island):
        """
        This just shifts all the values in the island by the pedestal
        value, to make the y axis scale a little more realistic. For
        right now, the pedestal value is just made randomly. As far as
        I know there is no fcl file with the pedestal values for each
        crystal, but if there are that would be better.
        """
        self.pedestal = np.random.normal(loc=1750, scale=2.5) * (-1)

        if(verbosity):
            print("The pedestal value of this island is " + \
                  "{0:.2f}".format(self.pedestal))
            
        return (island + self.pedestal)
    
    
    def integrate_(self, time, island):
        """
        Returns the integral of the island, the area under the curve
        """
        return (np.sum(island) * (time[1]-time[0]))
    

    def energize_(self, calibVal):
        """
        computes the total energy of the island, based on the calibration value
        """
        return (self.integral * calibVal)


    def __init__(self, referenceSpline,
                 minPulses = 1,
                 maxPulses = 4,
                 useUniformEscale = False,
                 minEscale = 50,
                 maxEscale=2400,
                 energyCalibrationVal = 0.3,
                 energyPeak = 0,
                 energyScale = 600,
                 deltaTmin = 0,
                 deltaTmax = 25,
                 minTimeOffset = 0,
                 gainSag = None,
                 verbosity = False,
                 normalize = True,
                 chop = True,
                 chopThreshold = 50,
                 nPreSamples = 8,
                 nPostSamples = 18,
                 noise = True,
                 noiseLevel = 4,
                 nTailSamples = 150,
                 randomizeTailLength = False,
                 samplingRateArtificial=1.25):
        """
        referenceSpline        : [spline] the spline template from which the
                                 pulses in this island are generated

        minPulses              : [int] the minimum number of pulses that might
                                 be present in this island. Defaults to 1

        maxPulses              : [int] the maximum number of pulses that might
                                 be present in this island. Defaults to 4

        useUniformEscale       : [boolean] if you want a uniform energy scaling
                                 rather than the default gaussian distribution

        minEscale              : [float] If using the uniform energy scaling,
                                 the minimum energy value possible for a pulse,
                                 in ADC units. Defaults to 50

        maxEscale              : [float] If using the uniform energy scaling,
                                 the maximum energy value possible for a pulse,
                                 in ADC units. Defaults to 2400

        energyCalibrationVal   : [float] the energy calibration constant of the
                                 crystal from which the spline templates for
                                 this island come from. This constant is used to
                                 convert a pulse integral made over the island
                                 in ADC units to an energy in MeV

        energyPeak             : [float] the peak of the gaussian distribution
                                 of pulse energies, in MeV. Defaults to zero

        energyScale            : [float] the standard deviation of energy of the
                                 pulses, in MeV. Defaults to 600

        deltaTmin              : [float] the smallest value on the time axis
                                 where a pulse might be placed. Defaults to
                                 0 ns
        deltaTmax              : [float] the largest value on the time axis
                                 where a pulse might be placed. Defaults to
                                 25 ns

        minTimeOffset          : [float] the artificial minimum time offset
                                 between pulses in this island in ns. Defaults
                                 to 0.0

        gainSag                : [None] short term double pulse correction, well
                                 mapped function, fix for non-recovered pixels:
                                 if super close and pixels haven't recovered, we
                                 might see a lower e value. Not super important,
                                 unless you want to also extract energies,
                                 because this has the main effect of just
                                 changing energy scale factors

        verbosity              : [boolean] if you want a bunch of output on the
                                 parameters of this island. Defaults to false

        normalize              : [boolean] if this island is normalized or not.
                                 Defaults to true

        chop                   : [boolean] if we want to chop this island or
                                 not, which is cutting it down in time to just
                                 focus on where the pulses are. This does not
                                 change the time or yValues parameters, but puts
                                 a chopped version into choppedTime and
                                 choppedYValues attributes. Defaults to True

        chopThreshold:         : [float] the threshold value in ADC units for
                                 when the chopping algorithm decides where a
                                 pulse is. Defaults to 50

        nPreSamples            : [int] the number of samples kept in chopping
                                 before the first value above chopThreshold.
                                 Defaults to 8

        nPostSamples           : [int] the number of samples kept in chopping
                                 after the last value above chopThreshold.
                                 Defaults to 18

        noise                  : [boolean] if this island has noise applied or
                                 not. Defaults to true

        noiseLevel             : [float] the standard deviation of gaussian
                                 noise in the sample, in ADC units. Defaults
                                 to 4

        nTailSamples           : [int] the number of tail samples in the
                                 island. defaults to 150

        randomizeTailLength    : [boolean] if this island has a randomized tail
                                 length or not. Defaults to false

        samplingRateArtificial : [float] the sampling rate of the island, which
                                 is used in the digitization process. To match
                                 the ADC, this defaults to 1.25ns and should
                                 probably be kept there
        """
        self.spline = referenceSpline
        self.nPulses = None
        self.energyCalibrationVal = energyCalibrationVal
        self.minTimeOffset = minTimeOffset
        self.gainSag = gainSag
        self.normalize = normalize
        self.chopThreshold = chopThreshold
        self.nPreSamples = nPreSamples
        self.nPostSamples = nPostSamples
        self.noise = noise
        self.noiseLevel = noiseLevel
        self.randomizeTailLength = randomizeTailLength
        self.nTailSamples = nTailSamples
        self.samplingRateArtificial = samplingRateArtificial
        self.pedestal = 0
        self.integral = None
        self.totalEnergy = None
        self.timeOffsets = []
        self.energyScaleFactors = []
        self.time = []
        self.yValues = []
        self.choppedTime = []
        self.choppedYValues = []
        
        """setup the spine variables"""

        splineTimes = self.spline.time
        splineShape = self.spline.splineY / self.spline.getPeak()
        splineSamplingRate = self.spline.samplingRate
        
        """figure out how many pulses there will be"""

        self.nPulses = self.giveNPulses_(verbosity, minPulses, maxPulses)
            
        """randomize the length of the tail samples"""

        if(randomizeTailLength):
            self.nTailSamples = self.giveRandomTailLength_(verbosity)
                
        """define the y and t arrays"""

        # before putting pulses in, yValues in thisIsland are just zeros
        thisIsland = np.append(np.zeros_like(self.spline.splineY),
                               np.zeros(self.nTailSamples))
        
        # the time array is the original one, plus
        # however much else we want based on nTailSamples
        theseTimes = np.append(splineTimes,
                               np.array([splineTimes[splineTimes.size-1] + \
                                        splineSamplingRate*i \
                                        for i in range(1,self.nTailSamples+1)]))
        
        """if nPulses is zero, we can skip almost everything"""

        if(self.nPulses == 0):
            # get the sampling right
            thisIsland, theseTimes = self.digitize_(verbosity,
                                                    thisIsland,
                                                    theseTimes)

            # add in our noise
            if(self.noise):
                thisIsland = self.giveNoise_(verbosity, thisIsland)

            # and now we're done
            self.time = theseTimes
            self.yValues = thisIsland

        else:
            """If we do have pulses, they need energy values"""

            if (useUniformEscale):
                self.energyScaleFactors = \
                self.giveUniformEnergyScaleFactors_(verbosity, 
                                                    minEscale, 
                                                    maxEscale)
            else:
                self.energyScaleFactors = \
                self.giveEnergyScaleFactors_(verbosity, 
                                             energyPeak, 
                                             energyScale,
                                             self.energyCalibrationVal)

            """create the time offsets for each pulse """

            self.timeOffsets = self.giveTimeOffsets_(verbosity,
                                                     deltaTmin,
                                                     deltaTmax)

            """put the pulses together into thisIsland"""

            for pulseIndex, deltaT in enumerate(self.timeOffsets):
                deltaT /= 2
                sample_offset = int(np.floor(deltaT))
                splineI = np.interp(splineTimes+deltaT-sample_offset,
                                    splineTimes, splineShape) * \
                          self.energyScaleFactors[pulseIndex]
                
                thisIsland[sample_offset:sample_offset + \
                           len(splineI)] += splineI

            """get the sampling right"""

            thisIsland, theseTimes = self.digitize_(verbosity,
                                                    thisIsland,
                                                    theseTimes)

            """add gaussian noise to each sample"""

            if(self.noise):
                thisIsland = self.giveNoise_(verbosity, thisIsland)
            
                
            """compute this island's integral and energy"""

            self.integral = self.integrate_(theseTimes, thisIsland)
            self.totalEnergy = self.energize_(self.energyCalibrationVal)
            
            """
            normalize if needed
            chopping threshold needs to be scaled accordingly
            """

            if(self.normalize):
                # if we are normalizing, our chopping threshold needs to
                # match our normalization
                thisIsland /= self.integral
                chopThreshold /= self.integral
            else:
                # if we aren't normalizing, we should make it look like
                # it came from an ADC, with a pedestal
                thisIsland = self.scaleADC_(verbosity, thisIsland)
                chopThreshold += self.pedestal
            
            """simulate island chopping, cut to around the pulses"""

            if (chop):
                chopTime, chopEn = self.chop_(theseTimes,
                                              thisIsland,
                                              chopThreshold,
                                              nPreSamples,
                                              nPostSamples)

                self.choppedTime = chopTime
                self.choppedYValues = chopEn

            """now we're good to go"""

            self.time = theseTimes
            self.yValues = thisIsland
            
            """lastly we can print a few parameters if requested"""

            if(verbosity):
                print("Spline template crystal number: "+ \
                      "{0}".format(self.spline.xtalNum))
                
                print("Spline template peak: " + \
                      "{0:.2f}".format(self.spline.getPeak()))
                
                print("Crystal energy calibration constant: " + \
                      "{0:.2f}".format(self.energyCalibrationVal))
                
                print("Pulse integral: " + \
                      "{0:.2f}".format(self.integral))
                
                print("Island total energy: " + \
                      "{0:.2f} MeV".format(self.totalEnergy))

"""
-------------------- Behold SplIland --------------------

contains the spline and island classes for making artificial
pileup islands. This code is the best thing to have ever been
made by any human being ever and if you disagree you can go eat
a cookie because you are probably right.

Written by Lars Borchert with extensive help from Josh LaBounty

---------------------------------------------------------
"""


import numpy as np
from scipy import stats


class Spline:
    """
    Take a root spline and make a simple,
    pythonified object with two np arrays.
    

    ------------ Attributes ------------
    
    splineY      : ndarray, the y values in the spline. They are
                   normalized such that their integral in time
                   is equal to one
                   
    time         : ndarray, the x values in the spline, which 
                   are time
    
    samplingRate : float, the sampling rate of the undigitized
                   spline
    """
    
    def peakTime(self):
        """
        Gives the time when a peak occurs
        The idea is to have a time for the center of a pulse
        """
        index = np.argmax(self.splineY)
        return(self.times[index])
    
    
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

    
    def __init__(self, rSpline):
        """
        rSpline : a ROOT spline, representing a normalized pulse
                  template
        """
        
        spline = []
        times = []
        
        for i in np.linspace(rSpline.GetXmin(), rSpline.GetXmax(), rSpline.GetNpx()):
            spline.append(rSpline.Eval(i))
            times.append(i)

        self.splineY = np.array(spline)
        self.times = np.array(times)
        
        self.samplingRate = times[2] - times[1]



class Island:
    """
    ------------ Attributes ------------
    
    
    
    """

    def digitize_(self, verbosity, thisIsland, theseTimes):
        """
        convert the timing to match that of the ADC
        by default this should be 1.25ns sampling
        """
        
        if(self.samplingRateArtificial is None):
            return thisIsland, theseTimes
        
        
        assert type(self.samplingRateArtificial) is float, print(
            """Error: Sampling rate must be a float.
            Currently {0}""".format(type(self.samplingRateArtificial)))

        if(verbosity):
            print(
            """Sampling this spline with a deltaT of {0} ns
            """.format(self.samplingRateArtificial) )
        
        sampledTimes = [theseTimes[0]]
    
        while(sampledTimes[-1] < [theseTimes[-1]]):
            sampledTimes.append(sampledTimes[-1] + self.samplingRateArtificial)
        
        theseSamples = np.interp(sampledTimes, theseTimes, thisIsland)
        
        thisIsland = np.array(theseSamples)
        theseTimes = np.array(sampledTimes)
        
        return thisIsland, theseTimes

    
    def makeDF(self):
        """
        return a simple pandas dataframe
        with time and energy columns for 
        this island
        """
        d = {'time': self.time,
             'energy': self.energy}
        df = pd.DataFrame(data = d)
        return(df)
    
    
    def giveNPulses_(self, verbosity, minPulses, maxPulses):
        """
        figure out how many pulses there will be
        never run this outside of initialization
        """
        self.nPulses = np.random.randint(minPulses, maxPulses+1)
        if(verbosity): print(self.nPulses, "pulses in this island")
            
            
    def giveRandomTailLength_(self, verbosity):
        """
        figure out how long the tail will be
        never run this outside of initialization
        """
        self.nTailSamples = self.nTailSamples + np.random.randint(-50, 100)
        if(verbosity): print(self.nTailSamples, " tail samples in this island")
            
    def giveTimeOffsets_(self, verbosity, deltaTmin, deltaTmax):
        """
        for the love of kale make deltaTmax >> minTimeOffset
        """
        offsets = np.random.uniform(deltaTmin, deltaTmax, size=self.nPulses)
        
        check = 0
        while(check < self.nPulses ** 2):
            offsets = np.random.uniform(deltaTmin, deltaTmax, size=self.nPulses)
            for i in range(0, len(offsets)):
                for j in range(0, len(offsets)):
                    if(i != j):
                        separation = offsets[i] - offsets[j]
                        if(abs(separation) < self.minTimeOffset):
                            check = 0
                        else:
                            check += 1
                    else:
                        check += 1
                        
        
        self.timeOffsets = offsets
        if(verbosity): print("The time offsets are: {0}".format(self.timeOffsets))
            
        
    def giveUniformEnergyScaleFactors_(self, verbosity, minEscale, maxEscale):
        """
        """
        self.energyScaleFactors = np.random.uniform(minEscale, maxEscale, size=self.nPulses)
        if(verbosity):
            print("The energy scale factors are: {0}".format(self.energyScaleFactors))
            
            
    def giveEnergyScaleFactors_(self, verbosity, energyPeak, energyScale, calibConstant):
        """
        """
        factors = []
        normalIntegral = self.spline.getNormalIntegral()
        
        for pulse in range(0, self.nPulses):
            # generate an energy, in Mev
            # we ignore stuff under a threshhold
            thisEnergy = 0
            while(thisEnergy < 100):
                thisEnergy = abs(np.random.normal(loc=energyPeak, scale=energyScale))

            thisHeight = thisEnergy / (calibConstant * normalIntegral)
            
            factors.append(thisHeight)
            
        self.energyScaleFactors = factors
        
        if(verbosity):
            print("The energy scale factors are: {0}".format(self.energyScaleFactors))
        
        
        # 170 MeV = calibration constant * pulse integral, by definition
        # calibIntegral = 170 / calibConstant
        
        
        # Energy is directly proportional to integral, is directly proportional to height
        # calibHeight / calibIntegral = normalHeight / normalIntegral
        # normalHeight == 1, by definition
        # calibHeight = calibIntegral / normalIntegral
        
        # Energy is directly proportional to integral, is directly proportional to height
        # thisHeight / calibHeight = thisEenrgy / calibEnergy
        # calibEnergy == 170 by definition
        # thisHeight = (thisEnergy *calibHeight) / (170)
        
        # calibHeight = calibIntegral / normalIntegral, from earlier
        # thisHeight = (thisEnergy*calibIntegral/NormalIntegral) / (170)
        
        # calibIntegral = 170 / calibConstant, from earlier
        # thisHeight = (thisEnergy*(170 / calibConstant)/normalIntegral) / (170)
        # thisHeight = thisEnergy / (calibConstant * normalIntegral)
        
        # Energy is directly proportional to integral, is directly proportional to height
        # thisHeight / normalHeight = thisEnergy / normalEnergy
        
        
        
    def giveNoise_(self, thisIsland):
        """
        """
        if(self.noise):
            thisIsland += np.random.normal(0, self.noiseLevel, size=thisIsland.size)
            return thisIsland
        
        
    def chop_(self, times, island, critValue, preSamples, postSamples):
        """
        """
        startIndex = None
        endIndex = None
        peakIndex = None
        cuts = []
        for index, val in enumerate(island):
            if(val > critValue):
                peakIndex = index
                cuts.append(index + postSamples)
                
                if(startIndex is None):
                    if(index - preSamples > 0):
                        startIndex = index - preSamples
                    else:
                        startIndex = 0
            """
            else:
                if((peakIndex is not None) and 
                   (index == postSamples + peakIndex)):
                    endIndex = index
                    cuts.append(endIndex)
            """
                    
        return times[startIndex:cuts[-1]], island[startIndex:cuts[-1]]
    
    
    def scaleADC_(self, verbosity, island):
        """
        """
        self.pedestal = np.random.normal(loc=1750, scale=2.5)
        
        if(verbosity):
            print("The pedestal value of this island is {0}".format(self.pedestal))
            
        return (island - self.pedestal)
    
    
    def integrate_(self, time, island):
        """
        """
        return (np.sum(island) * (time[1]-time[0]))
    
    
    def energize_(self, calibVal):
        """
        """
        return (self.integral * calibVal)


    def __init__(self, referenceSpline,
                 minPulses = 0, maxPulses = 4,
                 useUniformEscale = False,
                 minEscale = 50, maxEscale=2400,
                 energyCalibrationVal = 0.3,
                 energyPeak = 0, energyScale = 600,
                 deltaTmin = 0, deltaTmax = 25,
                 minTimeOffset = 0,
                 gainSag = None,
                 verbosity = False,
                 normalize = True,
                 chop = True,
                 chopThreshhold = 50,
                 nPreSamples = 8, nPostSamples = 18,
                 noise = False, noiseLevel = 4,
                 nTailSamples = 150, randomizeTailLength = True,
                 samplingRateArtificial=1.25):
        
        """
        referenceSpline:
        minPulses:
        maxPulses:
        minEscale:
        maxEscale:
        deltaTmin:
        deltaTmax:
        minRelativeTime:
        gainSag: short term double pulse correction, well mapped function, 
                 fix for non-recovered pixels, if super close and pixels
                 haven't recovered, we might see a lower e value. Not super
                 important, unless you want to also extract energies, because
                 this has the main effect of just changing energy scale factors
        verbosity:
        normalize:
        noise:
        noiseLevel:
        nTailSamples:
        randomizeTailLengte:
        samplingRateArtificial:
        """
        self.spline = referenceSpline
        self.minTimeOffset = minTimeOffset
        self.normalize = normalize
        self.noise = noise
        self.noiseLevel = noiseLevel
        self.randomizeTailLength = randomizeTailLength
        self.samplingRateArtificial = samplingRateArtificial
        self.nTailSamples = nTailSamples
        self.pedestal = 0
        
        splineTimes = self.spline.times
        splineShape = self.spline.energy / self.spline.getPeak()
        splineSamplingRate = self.spline.samplingRate
        
        # figure out how many pulses there will be
        self.giveNPulses_(verbosity, minPulses, maxPulses)
            
        # randomize the length of the tail samples
        if(randomizeTailLength):
            self.giveRandomTailLength_(verbosity)
                
        # define the y and t arrays that you see in the plot before filling them
        thisIsland = np.append(np.zeros_like(self.spline.energy), np.zeros(self.nTailSamples))
        
        # the time array is the original one, plus
        # however much else we want based on nTailSamples
        theseTimes = np.append(splineTimes,
                               np.array([splineTimes[splineTimes.size-1] + \
                                        splineSamplingRate*i \
                                        for i in range(1, self.nTailSamples+1)]))
        
        # if nPulses is zero, we can skip pulsing
        if(self.nPulses == 0):
            thisIsland, theseTimes = self.digitize_(verbosity, thisIsland, theseTimes)
            thisIsland = self.giveNoise_(thisIsland)
            self.energyScaleFactors = []
            self.timeOffsets = []
            self.time = theseTimes
            self.energy = thisIsland

        else:
            
            # create the energy scaling values
            if (useUniformEscale):
                self.giveUniformEnergyScaleFactors_(verbosity, minEscale, maxEscale)
            else:
                self.giveEnergyScaleFactors_(verbosity, energyPeak, energyScale, energyCalibrationVal)

            # create the time offsets
            self.giveTimeOffsets_(verbosity, deltaTmin, deltaTmax)

            # put the pulses together into thisIsland
            for pulseIndex, deltaT in enumerate(self.timeOffsets):
                sample_offset = int(np.floor(deltaT))
                splineI = np.interp(splineTimes+deltaT-sample_offset,
                                    splineTimes, splineShape) * \
                          self.energyScaleFactors[pulseIndex]
                
                thisIsland[sample_offset:sample_offset + len(splineI)] += splineI

            # add digitization correction, interpolate to 1.25 ns by default
            thisIsland, theseTimes = self.digitize_(verbosity, thisIsland, theseTimes)

            # add gaussian noise to each sample
            thisIsland = self.giveNoise_(thisIsland)
            
                
            # compute this island's integral and energy
            self.integral = self.integrate_(theseTimes, thisIsland)
            self.totalEnergy = self.energize_(energyCalibrationVal)
            
            # normalize if needed:
            if(self.normalize):
                integral = (np.sum(thisIsland) * (theseTimes[1]-theseTimes[0]))
                thisIsland /= integral
                chopThreshhold /= integral
            else:
                thisIsland = self.scaleADC_(verbosity, thisIsland)
                chopThreshhold -= self.pedestal
            
            # simulate island chopping, cut to around the pulses
            if (chop):
                chopTime, chopEn = self.chop_(theseTimes, thisIsland, chopThreshhold, nPreSamples, nPostSamples)

                self.choppedTime = chopTime
                self.choppedEnergy = chopEn

            # now we are good to go
            self.time = theseTimes
            self.energy = thisIsland
            
            
            
            
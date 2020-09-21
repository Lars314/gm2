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


class Spline:
    """
    Take a root spline and make a simple,
    pythonified object with two np arrays.
    
    energy : the energy axis
    time : the time axis
    """
    
    def peakTime(self):
        """
        Gives the time when a peak occurs
        The idea is to have a time for the center of a pulse
        """
        index = np.argmax(self.spline)
        return(self.times[index])

    
    def __init__(self, rSpline):
        spline = []
        times = []
        
        for i in np.linspace(rSpline.GetXmin(), rSpline.GetXmax(), rSpline.GetNpx()):
            spline.append(rSpline.Eval(i))
            times.append(i)

        self.energy = np.array(spline)
        self.times = np.array(times)
        
        self.samplingRate = times[2] - times[1]



class Island:
    """
    ------------ primary attributes ------------
    nPulses: int, the number of pulses in this island
    time: np.array, the time axis values
    energy: np.array, the energy axis values
    pulseTimes: np.array, the time offset of each pulse
    pulseEnegies: np.array, the relative energies of each pulse
    
    ----------- secondary attributes -----------
    noise: boolean, artificial noise off->false, on->true
    
    """

    def digitize_(self, verbosity, thisIsland, theseTimes):
        """
        convert the timing to what we want for actually
        applicable physics
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
            
        
    def giveEnergyScaleFactors_(self, verbosity, minEscale, maxEscale):
        """
        """
        self.energyScaleFactors = np.random.uniform(minEscale, maxEscale, size=self.nPulses)
        if(verbosity):
            print("The energy scale factors are: {0}".format(self.energyScaleFactors))
            
        
    def giveNoise_(self, thisIsland):
        """
        """
        if(self.noise):
            thisIsland += np.random.normal(0, self.noiseLevel, size=thisIsland.size)
            return thisIsland


    def __init__(self, referenceSpline,
                 minPulses = 0, maxPulses = 4,
                 minEscale = 1, maxEscale=10,
                 deltaTmin = 0, deltaTmax = 25,
                 minTimeOffset = 0,
                 gainSag = None,
                 verbosity = False,
                 normalize = True,
                 noise = False, noiseLevel = 0.001,
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
        gainSag:
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
        
        splineTimes = self.spline.times
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
            self.pulseEnergies = []
            self.pulseTimes = []
            self.time = theseTimes
            self.energy = thisIsland

        else:
            
            # create the energy scaling values
            self.giveEnergyScaleFactors_(verbosity, minEscale, maxEscale)

            # create the time offsets
            self.giveTimeOffsets_(verbosity, deltaTmin, deltaTmax)

            # put the pulses together into thisIsland
            for pulseIndex, deltaT in enumerate(self.timeOffsets):
                sample_offset = int(np.floor(deltaT))
                splineI = np.interp(splineTimes+deltaT-sample_offset,
                                    splineTimes, self.spline.energy) * \
                          self.energyScaleFactors[pulseIndex]
                thisIsland[sample_offset:sample_offset + len(splineI)] += splineI

            # normalize if needed:
            if(self.normalize):
                thisIsland /= np.sum(thisIsland)

            # add digitization correction, interpolate to 1.25 ns by default
            thisIsland, theseTimes = self.digitize_(verbosity, thisIsland, theseTimes)

            # add gaussian noise to each sample
            thisIsland = self.giveNoise_(thisIsland)

            # now we are good to go
            self.time = theseTimes
            self.energy = thisIsland
            
            
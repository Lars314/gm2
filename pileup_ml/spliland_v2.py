"""
----------------------------- Behold  CaloSplIland ----------------------------

contains the calo, spline, and island classes for making artificial
pileup islands. This code is the best thing to have ever been
made by any human being ever.

Written by Lars Borchert with extensive help from Josh LaBounty

This dashed line is 79 characters long
-------------------------------------------------------------------------------
"""
import numpy as np
import fclParse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import ROOT as r


class Spline:
    """
    Take a root spline and make a simple,
    pythonified object with two np arrays.


    ------------ Attributes ------------

    splineY      : [ndarray] the y values in the spline. They are
                   normalized such that their integral in time
                   is equal to one

    time         : [ndarray] the x values in the spline, which are time

    samplingRate : [float] the sampling rate of the undigitized spline

    xtalNum      : [int] the crystal number of the template
    """
# -----------------------------------------------------------------------------
# ----------------------------- Public  Functions -----------------------------
# -----------------------------------------------------------------------------

    def peakTime(self):
        """
        Gives the time when a peak occurs
        The idea is to have a time for the center of a pulse
        """
        index = np.argmax(self.trace)
        return(self.time[index])

    def getPeak(self):
        """
        Gives the value of the peak of the spline
        """
        return np.max(self.trace)

    def getNormalIntegral(self):
        """
        Gives the pulse integral of the peak-normalized spline
        Essentially, a unit pulse integral. peak-normalized in
        this case means that the peak of the spline is equal
        to one, rather than the integral of the spline equal to one
        """
        scaledVals = self.trace / self.getPeak()
        return np.sum(scaledVals) * self.samplingRate


# -----------------------------------------------------------------------------
# ----------------------------- Private Functions -----------------------------
# -----------------------------------------------------------------------------

    def __init__(self, pySpline=None,
                 rSpline=None, xtalNum=None, caloNum=None):

        """
        pySpline : a row of a pandas dataframe. I think this is a dictonary. It
                   can be indexed like one anyway.

        rSpline : a ROOT spline, representing a normalized pulse
                  template

        xtalNum : [int] the crystal number this spline is based from

        caloNum : [int] the calo number this spline is based from
        """
        if(pySpline is not None):
            self.trace = pySpline['trace']
            self.times = np.linspace(0, len(pySpline['trace']), 1)
            self.samplingRate = 1
            self.xtalNum = pySpline['xtalNum']
            self.caloNum = pySpline['caloNum']

        elif(rSpline is not None):
            trace = []
            times = []

            for i in np.linspace(rSpline.GetXmin(),
                                 rSpline.GetXmax(),
                                 rSpline.GetNpx()):
                trace.append(rSpline.Eval(i))
                times.append(i)

            self.trace = np.array(trace)
            self.time = np.array(times)

            self.samplingRate = times[2] - times[1]

            self.xtalNum = xtalNum
            self.caloNum = caloNum


class Xtal:
    """

    ------------ Attributes ------------

    caloNum : [int] the number of the calorimeter where this crystal lives

    x       : [int] the x position of this crystal, in crystal units

    y       : [int] the y position of this crystal, in crystal units

    eCalVal : [float] the energy calibration value for this crystal

    trace   :

    times   :

    impacts :
    """
# -----------------------------------------------------------------------------
# ----------------------------- Public  Functions -----------------------------
# -----------------------------------------------------------------------------
    def energize(self, impact_time, xtal_energy):
        """
        Put energy into this crystal. The energy is calculated at the
        calorimeter level
        """
        self.impacts.append([impact_time, xtal_energy])

    def build_trace(self):
        """
        Finalizes the crystal, gets the trace set up. In version one this
        was done in the constructor, but we need to be able to add energy
        values to this crystal before building it.
        """

        """randomize the length of the tail samples"""

        if(self.randomizeTailLength):
            self.nTailSamples = self.giveRandomTailLength_()

        # before putting pulses in, trace values are just zero
        thisTrace = np.append(np.zeros_like(self.template.trace),
                              np.zeros(self.nTailSamples))

        """setup the template variables"""

        templateTimes = self.template.time
        templateSamplingRate = self.template.samplingRate

        """define the trace and time arrays"""

        # the time array is the original one, plus however much else we want
        # based on nTailSamples
        theseTimes = np.append(templateTimes,
                               np.array([templateTimes[templateTimes.size-1] +
                                         templateSamplingRate*i
                                         for i in range(1,
                                                        self.nTailSamples+1)]))

        """Now we can put the pulses together"""

        templateShape = self.template.trace / self.template.getPeak()
        emptyTrace = thisTrace.copy()

        pulses = []
        for event in self.impacts:
            deltaT = event[0] / 2
            energy = event[1]
            height = self.convertToADC_(energy)

            sample_offset = int(np.floor(deltaT))     # redundant?
            splineI = np.interp(templateTimes + deltaT - sample_offset,
                                templateTimes, templateShape) * height

            thisTrace[sample_offset:sample_offset + len(splineI)] += splineI

            thisPulse = emptyTrace.copy()
            thisPulse[sample_offset:sample_offset + len(splineI)] += splineI
            pulses.append(thisPulse)

        """get the sampling right"""

        thisTrace, theseTimes = self.digitize_(self.verbosity,
                                               thisTrace,
                                               theseTimes)

        emptyTimes = theseTimes.copy()
        for pulse in pulses:
            thisPulse, someTimes = self.digitize_(False, pulse, emptyTimes)

        """add gaussian noise to each sample"""

        if (self.noise):
            thisTrace = self.giveNoise_(thisTrace)

        """compute this islands integral"""

        self.integral = self.integrate_(theseTimes, thisTrace)

        """normalize if needed"""

        if (self.normalize):
            # if we are normalizing, our chopping threshold needs to match
            # our normalization
            thisTrace /= self.integral
            self.chopThreshold /= self.integral
            self.pulses /= self.integral
        else:
            # if we don't want to normalize, we should make it look like it
            # came from an ADC, with a pedestal
            thisTrace = self.scaleADC_(thisTrace)

        """simulate island chopping, cut to around the pulses"""

        if (self.doChop):
            self.choppedTime, self.choppedTrace = self.chop_(theseTimes,
                                                             thisTrace)

        """now we're done!"""
        self.time = theseTimes
        self.trace = thisTrace


# -----------------------------------------------------------------------------
# ----------------------------- Private Functions -----------------------------
# -----------------------------------------------------------------------------

    def digitize_(self, verbosity, thisTrace, theseTimes):
        """
        convert the time sampling to match that of the ADC
        by default this should be 1.25ns sampling
        """
        if (self.samplingRateArtificial is None):
            return thisTrace, theseTimes

        assert type(self.samplingRateArtificial) is float, print(
            """Error: Sampling rate must be a float.
            Currently {0}""".format(type(self.samplingRateArtificial))
        )

        if (verbosity):
            print("Sampling this trace with a deltaT of " +
                  "{0:.2f} ns".format(self.samplingRateArtificial))

        sampledTimes = [theseTimes[0]]

        while(sampledTimes[-1] < [theseTimes[-1]]):
            sampledTimes.append(sampledTimes[-1] + self.samplingRateArtificial)

        theseSamples = np.interp(sampledTimes, theseTimes, thisTrace)

        thisTrace = np.array(theseSamples)
        theseTimes = np.array(sampledTimes)

        return thisTrace, theseTimes

    def giveRandomTailLength_(self):
        """
        figure out how long the tail will be
        never run this outside of initialization
        """
        nTailSamples = self.nTailSamples + np.random.randint(-50, 100)
        if(self.verbosity):
            print(nTailSamples, " tail samples in this island")
        return nTailSamples

    def convertToADC_(self, energy):
        """
        convert some energy to ADC units for this crystal
        """
        normalIntegral = self.template.getNormalIntegral()
        thisHeight = energy / (self.energyCalibrationVal * normalIntegral *
                               self.template.getPeak())
        return thisHeight

    def giveNoise_(self, thisTrace):
        """
        adds gaussian noise into the island, in ADC units. Do not use this
        function unless you are operating in ADC units
        """
        thisTrace += np.random.normal(0, self.noiseLevel, size=thisTrace.size)

        if (self.verbosity):
            print("The noise level standard deviation in this trace: " +
                  "{0:.2f}".format(self.noiseLevel))

        return thisTrace

    def chop_(self, times, trace):
        """
        Performs simulated island chopping. We don't want a bunch of zero
        values if we have many many data, so we cut islands
        down to just the interesting stuff. We keep preSamples before the
        first value over our threshold value, and postSamples after the
        last
        """
        startIndex = None
        endIndex = None
        cuts = []

        # loop through all the values
        for index, val in enumerate(trace):
            if (val > self.chopThreshold):
                # for something over our threshold, make an endIndex
                endIndex = index + self.nPostSamples

                # make sure we don't get out of bounds errors
                if (endIndex >= len(trace)):
                    cuts.append(len(trace))
                else:
                    cuts.append(endIndex)

                # if this is the first threshold passer, get the startIndex
                if (startIndex is None):
                    if(index - self.nPreSamples > 0):
                        startIndex = index - self.nPreSamples
                    else:
                        startIndex = 0

        # we only care about the endcaps, we have only one start value, but a
        # list of endIndex vales. But we only want the very last one in that
        # list
        return times[startIndex:cuts[-1]], trace[startIndex:cuts[-1]]

    def scaleADC_(self, trace):
        """
        This just shifts all the trace by the pedestal
        value, to make the y axis scale a little more realistic. For
        right now, the pedestal value is just made randomly. As far as
        I know there is no fcl file with the pedestal values for each
        crystal, but if there are that would be better.
        """
        self.pedestal = np.random.normal(loc=1750, scale=2.5) * (-1)

        if(self.verbosity):
            print("The pedestal value of this crystal is " +
                  "{0:.2f}".format(self.pedestal))

        return (trace + self.pedestal)

    def integrate_(self, time, trace):
        """
        Returns the integral of the trace, the area under the curve
        """
        return (np.sum(trace) * (time[1]-time[0]))

    def __init__(self, caloNum, x, y, xtalNum, normalize, chopThreshold,
                 nPreSamples, nPostSamples, noise, randomizeTailLength,
                 nTailSamples, samplingRateArtificial, pedestal, doChop,
                 verbosity):
        """
        caloNum : [int] the number of the calorimeter where this crystal lives

        x       : [int] the x position of this crystal, in crystal units

        y       : [int] the y position of this crystal, in crystal units

        eCalVal : [float] the energy calibration value for this crystal
        """
        self.caloNum = caloNum
        self.x = x
        self.y = y
        self.xtalNum = xtalNum
        self.impacts = []
        self.normalize = normalize
        self.chopThreshold = chopThreshold
        self.nPreSamples = nPreSamples
        self.nPostSamples = nPostSamples
        self.noise = noise
        self.randomizeTailLength = randomizeTailLength
        self.nTailSamples = nTailSamples
        self.samplingRateArtificial = samplingRateArtificial
        self.pedestal = pedestal
        self.doChop = doChop
        self.verbosity = verbosity

        self.time = []
        self.trace = []
        self.choppedTime = []
        self.choppedTrace = []

        # separated pulse traces, for debugging
        self.pulses = []

        self.integral = None

        self.noiseLevel = \
            fclParse.fclReader("gm2pedestals_run3formatRE.fcl")[
                               'pedestalConstantsLaserRun3']['calo15'][
                               'xtal' + xtalNum]['noiseLevel']

        self.energyCalibrationVal = \
            fclParse.fclReader("mipEnergyCalibration_" +
                               "PostDisk_CoincidencNumber_2.fcl")[
                               'absolute_calibration_constants'][
                               'calo'+str(caloNum)]['xtal'+str(xtalNum)]

        f = r.TFile("./calotemplate15.root")
        self.template = Spline(rspline=f.Get("masterSpline_xtal"+self.xtalNum),
                               xtalNum=self.xtalNum, caloNum=self.caloNum)


class Calorimeter:
    """
    Where the magic happens. Essentially just a 6x9 matrix of crystal objects.
    Each of those has a spline or none at some time, which is what goes into an
    island object. We probably don't need Calo and Island to be separate, but
    it makes it a litte easier for me to think about, and makes expansion easy

    ------------ Attributes ------------

    caloNum  : [int] the number of this calorimeter, 1 to 24

    xtalGrid : [2D list] a list of crystal objects

    """
# -----------------------------------------------------------------------------
# ----------------------------- Public  Functions -----------------------------
# -----------------------------------------------------------------------------

    def impact(self, p):
        """
        A positron hits. Send time and energy information to the crystals,
        making sure to send the correct energy to each crystal. Energies are
        determined by overlapping an area onto the calorimeter face. In this
        area, there is an energy density with some sort of falloff. The total
        energy is the energy passed to this function. Integrate the energy
        density over the area overlapping this energy impact area with each
        crystal, and tell that crystal that it got a signal with that energy at
        this time
        """
        impact_radius = 2 * self.moliere_radius
        self.impacts.append({"x": p.x,
                             "y": p.y,
                             "time": p.time,
                             "energy": p.energy,
                             "ring_x": impact_radius
                             * np.cos(np.linspace(0, 2 * np.pi, 17))
                             + p.x,
                             "ring_y": impact_radius
                             * np.sin(np.linspace(0, 2 * np.pi, 17))
                             + p.y,
                             "small_ring_x": self.moliere_radius
                             * np.cos(np.linspace(0, 2 * np.pi, 17))
                             + p.x,
                             "small_ring_y": self.moliere_radius
                             * np.sin(np.linspace(0, 2 * np.pi, 17))
                             + p.y})

    def clear(self):
        """
        Deletes all the impacts, only for debugging!!!
        """
        self.impacts = []
        self.reset_xtalGrid_()

    def getNonZeroXtals(self):
        """
        return the crystals in this calo that are not empty, to be stored in
        the island
        """

    def draw(self, legend=False, show_moliere_radius=False,
             show_hit_xtals=False, show_xtal_energies=False):
        """
        Give a matplotlib representation of this calorimeter, and show where
        the particles are, as well as their energies. This function is for
        debugging purposes
        """
        fig, ax = plt.subplots(1, 1)
        fig.set_size_inches(9, 6)

        for x in range(1, 9):
            ax.axvline(x=x, ymin=0, ymax=6, color='xkcd:grey')

        for y in range(1, 6):
            ax.axhline(y=y, xmin=0, xmax=9, color='xkcd:grey')

        for impact in self.impacts:
            this_color = next(ax._get_lines.prop_cycler)['color']
            ax.plot(impact["x"], impact["y"],
                    marker='o', markersize=10, color=this_color,
                    label="{0:.2f} MeV {1:.2f} ns".format(impact["energy"],
                                                          impact['time']))

            # ax.plot(impact["ring_x"], impact["ring_y"], color=this_color)
            # ax.plot(impact["small_ring_x"], impact["small_ring_y"],
            #        color=this_color)

            moliere_circle = plt.Circle((impact["x"], impact["y"]),
                                        self.moliere_radius, color=this_color,
                                        linewidth=3, fill=False)
            moliere_circle_2 = plt.Circle((impact["x"], impact["y"]),
                                          2*self.moliere_radius,
                                          color=this_color,
                                          linewidth=3, fill=False)

            if (show_hit_xtals):
                for hit_crystal in impact["hit_crystals"]:
                    rect = patches.Rectangle((hit_crystal[0], hit_crystal[1]),
                                             1, 1, linewidth=3,
                                             edgecolor='none',
                                             facecolor='xkcd:light grey')
                    ax.add_patch(rect)

            if (show_xtal_energies):
                for xtal in impact["hit_crystals"]:
                    xtal_events = self.xtalGrid[xtal[0]][xtal[1]].impacts
                    energies_list = []
                    for event in xtal_events:
                        energies_list.append(event[1])
                    energy_string = ""
                    for val in energies_list:
                        energy_string += "{0} MeV \n".format(int(val))
                    ax.text(x=xtal[0]+0.1, y=xtal[1],
                            s=energy_string)

            if (show_moliere_radius):
                ax.add_artist(moliere_circle)
                ax.add_artist(moliere_circle_2)

        ax.set_xlim(0, 9)
        ax.set_ylim(0, 6)
        if (legend):
            ax.legend()
        # return fig

# -----------------------------------------------------------------------------
# ----------------------------- Private Functions -----------------------------
# -----------------------------------------------------------------------------
    def get_energy_density_(self, r):
        """
        The energy density of the impacted particle over the region of impact.
        This will change, 1/r is definitely not right but easy for testing
        """
        return 1 / (4 * np.pi * r**2)

    def get_hit_crystals_(self):
        """
        """
        # first, let's just flag where the moliere circle overlaps crystals
        # not all of these crystals will end up mattering, as it is likely
        # that the edge ones will get less than 50 MeV

        # loop through the impacts
        for impact in self.impacts:
            impact["hit_crystals"] = []
            # loop for the 2nd Moliere radius
            for i in range(0, len(impact["ring_x"])):
                # for this location on the Moliere ring, take note of the
                # crystal coordinate
                this_x = int(impact["ring_x"][i])
                this_y = int(impact["ring_y"][i])
                # check that this crystal actually exists on the calo face
                if ((this_x >= 9) or (this_y >= 6)):
                    continue
                # check if the crystal has already been counted, count it
                if ((this_x, this_y) not in impact["hit_crystals"]):
                    impact["hit_crystals"].append((this_x, this_y))
            # loop for the 1st Moliere radius
            for i in range(0, len(impact["small_ring_x"])):
                # for this location on the Moliere ring, take note of the
                # crystal coordinate
                this_x = int(impact["small_ring_x"][i])
                this_y = int(impact["small_ring_y"][i])
                # check that this crystal actually exists on the calo face
                if ((this_x >= 9) or (this_y >= 6)
                   or (this_x < 0) or (this_y < 0)):
                    continue
                # check if the crystal has already been counted, count it
                if ((this_x, this_y) not in impact["hit_crystals"]):
                    impact["hit_crystals"].append((this_x, this_y))

    def get_crystal_energy_(self, xtal_loc, impact_loc, impact_energy):
        """
        For some crystal, hit location, and hit energy, we calculate how much
        energy was deposited into this crystal.

        Currently this method is inaccurate. We are assuming that the energy
        is distributed evenly accross the face of the crystal, with the energy
        density of the center of the crystal. In reality, we would need to
        integrate the energy density of the impact over the surface of the
        crystal. Additional complications come from energy loss at the borders
        of crystals. These are ignored. We can do this for now because we are
        mostly interested in the number of impacts, not the energy.
        """
        imp_x = impact_loc[0]
        imp_y = impact_loc[1]
        xtal_x = xtal_loc[0] + 0.5
        xtal_y = xtal_loc[1] + 0.5
        r = np.sqrt(((xtal_x - imp_x)**2) + ((xtal_y - imp_y)**2))
        e_density = self.get_energy_density_(r)

        # crystal has area of 1
        return e_density * impact_energy

    def build_crystals_(self):
        """
        Using our lists of impacted crystals, create crystal objects with the
        correct energy depositions. For multiple impacts, crystals may appear
        twice, once in each impact list. But this is what we want. The energy
        of two impacts on one crystal will add in that crystal, so we
        definitely want to look at both situations.
        """
        # go through each impact
        for impact in self.impacts:
            # if we have already sent information from this impact to the
            # crystals, we should not do it again
            if impact['time'] not in self.impacts_registered_by_xtals:
                self.impacts_registered_by_xtals.append(impact['time'])
                # go through each crystal
                for xtal in impact['hit_crystals']:
                    impact_loc = (impact['x'], impact['y'])
                    xtal_energy = self.get_crystal_energy_(xtal, impact_loc,
                                                           impact['energy'])
                    self.xtalGrid[xtal[0]][xtal[1]].energize(impact['time'],
                                                             xtal_energy)

    def reset_xtalGrid_(self):
        """
        Creates a fresh grid of crystals. Use when initializing the calo, or
        when clearing it of impacts.
        """
        xtalGrid = []
        for column in range(0, 9):
            this_column = []
            for row in range(0, 6):
                this_num = (column * 6) + (row)
                this_xtal = Xtal(self.caloNum, x=column+1, y=row+1,
                                 xtalNum=this_num)
                this_column.append(this_xtal)
            xtalGrid.append(this_column)
        self.xtalGrid = xtalGrid

    def __init__(self, caloNum):
        """
        caloNum  : [int] the number of this calorimeter, 1 to 24
        """
        self.caloNum = caloNum
        self.moliere_radius = 1.8 / 2.54
        self.impacts = []
        self.impacts_registered_by_xtals = []

        # create our grid
        self.reset_xtalGrid_()

    def __repr__(self):
        summary = ""

        for columm in self.xtalGrid:
            for xtal in columm:
                this_report = "Crystal {0} at {1} has \
[time (ns), energy (MeV)]: {2}".format(xtal.xtalNum,
                                       (xtal.x, xtal.y),
                                       xtal.impacts)
                summary += this_report + "\n\n"

        return summary


class Island:
    """

    ------------ Attributes ------------

    """
# -----------------------------------------------------------------------------
# ----------------------------- Public  Functions -----------------------------
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# ----------------------------- Private Functions -----------------------------
# -----------------------------------------------------------------------------

    def giveNPulses_(self, verbosity, minPulses, maxPulses):
        """
        figure out how many pulses there will be
        never run this outside of initialization
        """
        nPulses = np.random.randint(minPulses, maxPulses+1)
        if(verbosity):
            print("The number of pulses in this island: " +
                  "{0}".format(nPulses))
        return nPulses

    def giveTimeOffsets_(self, verbosity, deltaTmin, deltaTmax):
        """
        Creates the time offset values for each pulse
        for the love of kale make deltaTmax >> minTimeOffset, otherwise
        this will take absolutely forever to run
        """

        check = 0
        # keep looking until every combination of offsets has passed the check
        while(check < self.nParticles ** 2):
            # make some random time offsets
            offsets = np.random.uniform(deltaTmin, deltaTmax,
                                        size=self.nParticles)
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
            print("The artificial minimum time offset is: " +
                  "{0:.2f} ns".format(self.minTimeOffset))
        return offsets

    def __init__(self,
                 caloNum=1,
                 minPulses=1,
                 maxPulses=4,
                 deltaTmin=0,
                 deltaTmax=25,
                 minTimeOffset=0,
                 gainSag=None,
                 verbosity=False,
                 normalize=True,
                 chop=True,
                 chopThreshold=50,
                 nPreSamples=8,
                 nPostSamples=18,
                 noise=True,
                 nTailSamples=150,
                 randomizeTailLength=False,
                 samplingRateArtificial=1.25,
                 sample_dataframe=None):
        """

        """

        # decide on particle number
        self.nParticles = self.giveNPulses_(verbosity, minPulses, maxPulses)

        # get particle times
        """the particle class can get its own times, but it does not know about
        any restrictions on the times by minTimeOffset. It can only get its
        own times for debugging purposes"""
        self.minTimeOffset = minTimeOffset
        self.timeOffsets = self.giveTimeOffsets_(verbosity,
                                                 deltaTmin, deltaTmax)

        # get particles
        self.particles = []
        for this_time in self.timeOffsets:
            this_particle = Particle(time=this_time, e_df=sample_dataframe)
            self.particles.append(this_particle)

        # make a calorimeter
        self.calo = Calorimeter(caloNum)

        # particles impact the calorimeter
        for particle in self.particles:
            self.calo.impact(x=particle.x, y=particle.y,
                             energy=particle.energy)

        # get the crystals from the calorimeter, save positions and traces


class Particle:
    """
    Represents a particle that hits the calorimeter, and the information
    relevant to that event. Position on the calo face and energy are considered
    independent here, even though that is not true in the real experiment. But
    it is a good enough assumption here for now

    ------------ Attributes ------------
    x      : [float] the x position of the particle on the calo face, in calo
             crystal units

    y      : [float] the y position of the particle on the calo face, in calo
             crystal units

    energy : [float] the energy of this particle

    time   : [float] the time this particle hits the calorimeter
    """
# -----------------------------------------------------------------------------
# ----------------------------- Public  Functions -----------------------------
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# ----------------------------- Private Functions -----------------------------
# -----------------------------------------------------------------------------

    def __init__(self,
                 x=None,
                 y=None,
                 energy=None,
                 time=None,
                 deltaTmin=0,
                 deltaTmax=25,
                 e_df=None):
        """
        x      : [float] the x position of the particle on the calo face,
                 in calo crystal units

        y      : [float] the y position of the particle on the calo face,
                 in calo crystal units

        energy : [float] the energy of this particle

        time   : [float] the time this particle hits the calorimeter
        """
        if (x is None):
            self.x = np.random.uniform(low=0.0, high=9.0)
        else:
            self.x = x

        if (y is None):
            self.y = np.random.uniform(low=0.0, high=6.0)
        else:
            self.y = y

        if (energy is None):
            if (e_df is None):
                print("No particle energy provided, " +
                      "defaulting to energy of 3.094 GeV")
                self.energy = 3094
            else:
                index = np.random.randint(low=0, high=len(e_df))
                self.energy = e_df[index]
        else:
            self.energy = energy

        if (time is None):
            self.time = np.random.uniform(low=deltaTmin, high=deltaTmax)
        else:
            self.time = time

    def __repr__(self):
        summary = "\n----- New Particle -----\n"

        summary += "(x, y) : ({0:.2f}, {1:.2f})\n".format(self.x, self.y)
        summary += "Energy : {0:.2f} MeV\n".format(self.energy)
        summary += "Time : {0:.2f} ns\n".format(self.time)

        return summary
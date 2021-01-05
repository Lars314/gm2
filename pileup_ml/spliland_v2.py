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
from scipy import integrate


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

    trace

    times

    energies :
    """
# -----------------------------------------------------------------------------
# ----------------------------- Public  Functions -----------------------------
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# ----------------------------- Private Functions -----------------------------
# -----------------------------------------------------------------------------

    def digitize():
        """

        """

    def giveRandomTailLength_():
        """

        """

    def giveNoise_():
        """

        """

    def chop_():
        """

        """

    def scaleADC():
        """

        """

    def integrate_():
        """

        """

    def energize_():
        """

        """

    def __init__(self, caloNum, x, y, xtalNum):
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

        self.energyCalibrationVal = \
            fclParse.fclReader("mipEnergyCalibration_" +
                               "PostDisk_CoincidencNumber_2.fcl")[
                               'absolute_calibration_constants'][
                               'calo'+str(caloNum)]['xtal'+str(xtalNum)]


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

    def getNonZeroXtals(self):
        """
        return the crystals in this calo that are not empty, to be stored in
        the island
        """

    def draw(self, legend=False, do_moliere_radius=False,
             show_hit_xtals=False):
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
                    label="{0:.2f} MeV".format(impact["energy"]))

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

            if (do_moliere_radius):
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
    def energy_density_(self, r):
        """
        The energy density of the impacted particle over the region of impact.
        This will change, 1/r is definitely not right but easy for testing
        """
        return 1 / r

    def get_hit_crystals_(self):
        """
        """
        # first, let's just flag where the moliere circle overlaps crystals
        # not all of these crystals will end up mattering, as it is likely
        # that the edge ones will get less than 50 MeV
        for impact in self.impacts:
            impact["hit_crystals"] = []
            for i in range(0, len(impact["ring_x"])):
                this_x = int(impact["ring_x"][i])
                this_y = int(impact["ring_y"][i])
                if ((this_x, this_y) not in impact["hit_crystals"]):
                    impact["hit_crystals"].append((this_x, this_y))
            for i in range(0, len(impact["small_ring_x"])):
                this_x = int(impact["small_ring_x"][i])
                this_y = int(impact["small_ring_y"][i])
                if ((this_x, this_y) not in impact["hit_crystals"]):
                    impact["hit_crystals"].append((this_x, this_y))
            """
            # now that we have our overlaps, we need to go into each one and
            # integrate with our energy density over the area of the overlap
            f = lambda y, x: 1 / np.sqrt((x**2) + (y**2))
            for crystal in impact["hit_crystals"]:
                const_A = np.sqrt(((crystal[0] - impact["x"])**2)
                                  + ((crystal[1] - impact["y"])**2))

                energy = integrate.dblquad(f, 0, x_c,
                                           lambda x: 0, lambda x: y_c)
            """

    def __init__(self, caloNum):
        """
        caloNum  : [int] the number of this calorimeter, 1 to 24
        """
        # create our grid
        xtalGrid = []
        for column in range(0, 9):
            this_column = []
            for row in range(0, 6):
                this_xtal = Xtal(caloNum, x=column+1, y=row+1,
                                 xtalNum=1)
                this_column.append(this_xtal)
            xtalGrid.append(this_column)

        self.xtalGrid = xtalGrid
        self.caloNum = caloNum
        self.moliere_radius = 1.8 / 2.54
        self.impacts = []


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
                print("No particle energy provided," +
                      "defaulting to energy of 1MeV")
                self.energy = 2400
            else:
                index = np.random.randint(low=0, high=len(e_df))
                self.energy = e_df[index]
        else:
            self.energy = energy

        if (time is None):
            self.time = np.random.uniform(low=deltaTmin, high=deltaTmax)
        else:
            self.time = time

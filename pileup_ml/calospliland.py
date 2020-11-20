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


class Calo:
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

# -----------------------------------------------------------------------------
# ----------------------------- Private Functions -----------------------------
# -----------------------------------------------------------------------------

    def __init__(self, caloNum):
        """
        caloNum  : [int] the number of this calorimeter, 1 to 24
        """
        xtalGrid = []
        for column in range(0, 9):
            this_column = []
            for row in range(0, 6):
                this_xtal = Xtal(caloNum, x=column+1, y=row+1)
                this_column.append(this_xtal)
            xtalGrid.append(this_column)

        self.xtalGrid = xtalGrid
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
    """
# -----------------------------------------------------------------------------
# ----------------------------- Public  Functions -----------------------------
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# ----------------------------- Private Functions -----------------------------
# -----------------------------------------------------------------------------

    def __init__(self, caloNum, x, y, eCalVal):
        """
        caloNum : [int] the number of the calorimeter where this crystal lives

        x       : [int] the x position of this crystal, in crystal units

        y       : [int] the y position of this crystal, in crystal units

        eCalVal : [float] the energy calibration value for this crystal
        """
        self.caloNum = caloNum
        self.x = x
        self.y = y


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

    def __init__(self):
        """

        """

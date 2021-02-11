# this file has no function with respect to being run. The sole purpose is to
# be a blank space for me to copy and paste stuff to for reference of originals
# while editing functions. 
   
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

    if (self.verbosity):
        print("This trace was chopped with threshold " +
                "{0:.2f}\n".format(self.chopThreshold) +
                "and with chop indices {0}".format(startIndex) +
                " and {0}".format(cuts[-1]))

    # we only care about the endcaps, we have only one start value, but a
    # list of endIndex vales. But we only want the very last one in that
    # list
    try:
        return times[startIndex:cuts[-1]], trace[startIndex:cuts[-1]]
    except(IndexError):
        # this means no end cut was found ever, so no value over the
        # threshold exists. We don't want to keep this crystal
        self.clear_()
        return [], []
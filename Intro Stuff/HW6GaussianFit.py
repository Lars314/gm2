import ROOT as r
import pandas as pd
import matplotlib.pyplot as plt

def gaus_fit(hist, nbins, step):
    """
    Performs a Gaussian fit on 
    """
    stats = []
    
    for index in range(0, hist.GetNbinsX()-step, step):
        proj = hist.ProjectionY("_py", index, index+step)
        the_fit = proj.Fit("gaus", "SN")
        stats.append({
            'TimeBin':   index*(hist.GetXaxis().GetXmax()/nbins),
            'GausConst': the_fit.Parameter(0),
            'GausMean':  the_fit.Parameter(1),
            'GausSD':    the_fit.Parameter(2),
            'HistMean':  proj.GetMean(),
            'HistSD':    proj.GetStdDev()
        })

        
def plot_fit_results(data):
    fig, ax = plt.subplots(1,2)
    fig.set_size_inches(10,5)

    ax[0].plot(data['TimeBin'], data['GausMean'], linestyle='none', marker='o',
               color='xkcd:blue', markersize=4, label='Gaus Fit');
    ax[0].plot(data['TimeBin'], data['HistMean'], linestyle='none', marker='o',
               color='xkcd:red', markersize=4, label='Histogram Stat');
    ax[0].set_ylabel('Mean');
    ax[0].set_xlabel('Time');
    ax[0].legend(loc=0,shadow=True);

    ax[1].plot(data['TimeBin'], data['GausSD'], linestyle='none', marker='o',
               color='xkcd:blue', markersize=4, label='Gaus Fit');
    ax[1].plot(data['TimeBin'], data['HistSD'], linestyle='none', marker='o',
               color='xkcd:red', markersize=4, label='Histogram Stat');
    ax[1].set_ylabel('Standard Deviation');
    ax[1].set_xlabel('Time');
    ax[1].legend(loc=0,shadow=True);

    fig.tight_layout()
    plt.savefig('test.png', bbox_inches='tight')
    
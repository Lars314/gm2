import numpy as np
import pandas as pd
import lmfit


def linear(t, c, m):
    return(c + m * t)


def single_exponential(t, c, a, tau, m):
    return(c - a * (np.exp(-t/tau)) + m*t)


def double_exponential(t, c, a, tauA, b, tauB):
    return(c - a * (np.exp(-t/tauA)) + b * (np.exp(t/tauB)))


class DataSet:
    """
    """
    def one_linear_fit_(self, df, fit_lim):
        """
        """
        model = lmfit.Model(linear)

        mean_params = {}

        fit_data = df[(df['TimeBin'] > fit_lim[0])
                      & (df['TimeBin'] < fit_lim[1])][:]

        # Mean
        mean_result = model.fit(fit_data['Mean'],
                                t=fit_data['TimeBin'],
                                c=0.5, a=8, tau=7, m=0.0001,
                                weights=1/fit_data['MeanError'])
        mean_cov = mean_result.covar

        for name, param in mean_result.params.items():
            mean_params[name] = param.value

        for index, name in enumerate(mean_result.params):
            mean_params[name + '_err'] = np.sqrt(mean_cov[index][index])

        mean_params['redchi'] = mean_result.redchi

        return mean_params

    def calo_df_(self, hist, step):
        """
        """
        data = pd.DataFrame(columns=['TimeBin', 'Mean', 'MeanError', 'SD',
                                     'SDError', 'Var', 'VarError', 'Npoints',
                                     'RMS', 'RMSError'])

        for index in range(0, hist.GetNbinsX()-step, step):
            a = dict()

            proj = hist.ProjectionY("_py", index, index+step)

            a['TimeBin'] = hist.GetXaxis().GetBinCenter(index)
            a['Mean'] = (proj.GetMean() - 3) * 25.2
            a['MeanError'] = proj.GetMeanError() * 25.2
            a['SD'] = proj.GetStdDev() * 25.2
            a['SDError'] = proj.GetStdDevError() * 25.2
            a['Var'] = (proj.GetStdDev() ** 2) * (25.2 ** 2)
            a['VarError'] = 2 * abs(proj.GetStdDev()) \
                * proj.GetStdDevError() * (25.2 ** 2)
            a['Npoints'] = proj.GetEntries()
            a['RMS'] = proj.GetRMS() * 25.2
            a['RMSError'] = proj.GetRMSError() * 25.2

            data.loc[index] = a

        return data

    def build_hists_(self):
        calo_hists = {}
        rand_hists = {}

        for i in range(1, 25):
            calo = self.file.Get("verticalPosition/clusters"
                                 + str(i)).Clone(self.name + "_calo_" + str(i))
            rand = self.file.Get("verticalPosition/randclusters"
                                 + str(i)).Clone(self.name
                                                 + "_rand_calo_" + str(i))

            calo.SetAxisRange(self.energy_range[0], self.energy_range[1], "y")
            rand.SetAxisRange(self.energy_range[0], self.energy_range[1], "y")
            calo_hists[i] = calo.Project3D("zx")
            rand_hists[i] = rand.Project3D("zx")

        return {"unrandomized": calo_hists, "randomized": rand_hists}

    def build_dfs_(self, step):
        unrandomized = {}
        randomized = {}

        for caloNum, hist in self.hists['unrandomized'].items():
            unrandomized[caloNum] = self.calo_df_(hist, step)

        for caloNum, hist in self.hists['randomized'].items():
            randomized[caloNum] = self.calo_df_(hist, step)

        return {"unrandomized": unrandomized, "randomized": randomized}

    def linear_fit(self, fit_lim):
        """
        """
        unrandomized = {}
        randomized = {}

        for caloNum, df in self.dfs['unrandomized'].items():
            unrandomized[caloNum] = self.one_linear_fit_(df, fit_lim)

        for caloNum, df in self.dfs['randomized'].items():
            randomized[caloNum] = self.one_linear_fit_(df, fit_lim)

        self.linear_parameters = {"unrandomized": unrandomized,
                                  "randomized": randomized}

    def __init__(self,
                 name=None,
                 long_name=None,
                 file=None,
                 energy_range=[0, 99999999999],
                 step=10,
                 fit_lim=[50, 200],
                 output_directory=None):
        self.name = name
        self.long_name = long_name
        self.file = file
        self.step = step
        self.energy_range = energy_range
        self.linear_parameters = None
        self.output_directory = output_directory

        self.hists = self.build_hists_()
        self.dfs = self.build_dfs_(step=self.step)

        # self.linear_fit(fit_lim)

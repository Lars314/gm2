import ROOT as r
import pandas as pd
import math
import numpy as np
import os
import lmfit
import fit_functions as funcs

    
    
class CalorimeterX:
    """
    
    """
    def __init__(self, number=None, step=None,
                 hists={'randomized':None, 'unrandomized':None}):
        self.number          = number
        self.step            = step
        self.hists           = hists
        self.data            = None
        self.fft_data        = None
        self.fit_lim_0       = []
        self.fit_lim_1       = []
        self.fit_lim_2       = []
        self.func0_params    = {}
        self.func1_params    = {}
        self.func2_params    = {}
        self.func0_residuals = {}
        self.func1_residuals = {}
        self.func2_residuals = {}
        self.func2_result = None
        
     
    def _hist_to_df(self, hist):
        """
        """
        df = pd.DataFrame(columns=['TimeBin', 'Mean', 'MeanError', 'SD', 
                                   'SDError', 'Var', 'VarError','Npoints', 
                                   'RMS', 'RMSError'])
        
        for index in range(0, hist.GetNbinsX()-self.step, self.step):
            a = dict()

            proj = hist.ProjectionY("_py", index, index+self.step)

            a['TimeBin']      = hist.GetXaxis().GetBinCenter(index)
            a['Mean']         = (proj.GetMean() - 3) * 25.2
            a['MeanError']    = proj.GetMeanError() * 25.2
            a['SD']           = proj.GetStdDev() * 25.2
            a['SDError']      = proj.GetStdDevError() * 25.2
            a['Var']          = (proj.GetStdDev() ** 2) * (25.2 ** 2)
            a['VarError']     = 2 * abs(proj.GetStdDev()) * proj.GetStdDevError() * (25.2 ** 2)
            a['Npoints']      = proj.GetEntries()
            a['RMS']          = proj.GetRMS() * 25.2
            a['RMSError']     = proj.GetRMSError() * 25.2

            df.loc[index] = a
                
        return df
        

    def build_dfs(self):
        """
        """
        nora_df = self._hist_to_df(self.hists['unrandomized']).reset_index()
        rand_df = self._hist_to_df(self.hists['randomized']).reset_index()
        
        self.data = {'unrandomized':nora_df, 'randomized':rand_df}
    
    
    def _do_fit_0(self, df, guess, hints):
        """
        """
        # Get things ready
        model = lmfit.Model(funcs.func0)
        rms_params = {}
        fit_lim = self.fit_lim_0
        
        model.set_param_hint('c', min=hints['c'][0], max=hints['c'][1])
        model.set_param_hint('m', min=hints['m'][0], max=hints['m'][1])
       
        fit_data = df[(df['TimeBin'] > fit_lim[0]) & (df['TimeBin'] < fit_lim[1]) ][:]
        
        # Actually do the fit
        rms_result = model.fit(fit_data['RMS'], t=fit_data['TimeBin'],
                                c=guess['c'], m=guess['m'],
                                weights=1/fit_data['RMSError'])
        rms_cov = rms_result.covar
        
        # Go through the results, and format them into a dictionary
        for name, param in rms_result.params.items():
            rms_params[name] = param.value
    
        for index, name in enumerate(rms_result.params):
            rms_params[name + '_err'] = np.sqrt(rms_cov[index][index])

        rms_params['redchi'] = rms_result.redchi
        
        # return the formatted dictionary
        return rms_params, rms_result
    
    
    def _do_fit_1(self, df, guess, hints):
        """
        """
        # Get things ready
        model = lmfit.Model(funcs.func1)
        rms_params = {}
        fit_lim = self.fit_lim_1
        
        model.set_param_hint('c', min=hints['c'][0], max=hints['c'][1])
        model.set_param_hint('a', min=hints['a'][0], max=hints['a'][1])
        model.set_param_hint('tau', min=hints['tau'][0], max=hints['tau'][1])
        model.set_param_hint('m', min=hints['m'][0], max=hints['m'][1])
       
        
        fit_data = df[(df['TimeBin'] > fit_lim[0]) & (df['TimeBin'] < fit_lim[1]) ][:]
        
        # Actually do the fit
        rms_result = model.fit(fit_data['RMS'], t=fit_data['TimeBin'],
                                c=guess['c'], a=guess['a'], tau=guess['tau'], m=guess['m'],
                                weights=1/fit_data['RMSError'])
        rms_cov = rms_result.covar
        
        # Go through the results, and format them into a dictionary
        for name, param in rms_result.params.items():
            rms_params[name] = param.value
    
        for index, name in enumerate(rms_result.params):
            rms_params[name + '_err'] = np.sqrt(rms_cov[index][index])

        rms_params['redchi'] = rms_result.redchi
        
        # return the formatted dictionary
        return rms_params, rms_result
    
    
    def _do_fit_2(self, df, guess, hints):
        """
        """
        # Get things ready
        model = lmfit.Model(funcs.func2)
        mean_params = {}
        fit_lim = self.fit_lim_2
        
        #model.set_param_hint('c', min=hints['c'][0], max=hints['c'][1])
        model.set_param_hint('a', min=hints['a'][0], max=hints['a'][1])
        model.set_param_hint('tauA', min=hints['tauA'][0], max=hints['tauA'][1])
        model.set_param_hint('b', min=hints['b'][0], max=hints['b'][1])
        model.set_param_hint('tauB', min=hints['tauB'][0], max=hints['tauB'][1])
        
        fit_data = df[(df['TimeBin'] > fit_lim[0]) & (df['TimeBin'] < fit_lim[1]) ][:]
        
        # Actually do the fit
        mean_result = model.fit(fit_data['Mean'], t=fit_data['TimeBin'],
                                a=guess['a'], tauA=guess['tauA'], b=guess['b'], tauB=guess['tauB'],
                                weights=1/fit_data['MeanError'])
        mean_cov = mean_result.covar
        
        # Go through the results, and format them into a dictionary
        for name, param in mean_result.params.items():
            mean_params[name] = param.value
    
        for index, name in enumerate(mean_result.params):
            mean_params[name + '_err'] = np.sqrt(np.abs(mean_cov[index][index]))

        mean_params['redchi'] = mean_result.redchi
        
        # Return the formatted dictionary
        return mean_params, mean_result
    
    
    def build_params_0(self, fit_lim, param_guesses, param_hints):
        """
        """
        self.fit_lim_0 = fit_lim
        nora_params = self._do_fit_0(df=self.data['unrandomized'], guess=param_guesses, hints=param_hints)
        rand_params = self._do_fit_0(df=self.data['randomized'], guess=param_guesses, hints=param_hints)
        
        self.func0_params = {'unrandomized':nora_params[0], 'randomized':rand_params[0]}
        self.func0_result = {'unrandomized':nora_params[1], 'randomized':rand_params[1]}
    
     
    def build_params_1(self, fit_lim, param_guesses, param_hints):
        """
        """
        self.fit_lim_1 = fit_lim
        nora_params = self._do_fit_1(df=self.data['unrandomized'], guess=param_guesses, hints=param_hints)
        rand_params = self._do_fit_1(df=self.data['randomized'], guess=param_guesses, hints=param_hints)
        
        self.func1_params = {'unrandomized':nora_params[0], 'randomized':rand_params[0]}
        self.func1_result = {'unrandomized':nora_params[1], 'randomized':rand_params[1]}
    
    
    def build_params_2(self, fit_lim, param_guesses, param_hints):
        """
        """
        self.fit_lim_2 = fit_lim
        nora_params = self._do_fit_2(df=self.data['unrandomized'], guess=param_guesses, hints=param_hints)
        rand_params = self._do_fit_2(df=self.data['randomized'], guess=param_guesses, hints=param_hints)
        
        self.func2_params = {'unrandomized':nora_params[0], 'randomized':rand_params[0]}
        self.func2_result = {'unrandomized':nora_params[1], 'randomized':rand_params[1]}
        

    def _find_residuals_1(self, df, params):
        """
        """
        residuals = pd.DataFrame(columns=['TimeBin', 'Residual'])
        low  = self.fit_lim_1[0]
        high = self.fit_lim_1[1]
        
        residuals['TimeBin'] = df[(df['TimeBin'] > low) & (df['TimeBin'] < high)]['TimeBin']
                
        expected = residuals['TimeBin'].apply(funcs.func1,
                                              args=(params['c'], params['a'], 
                                                    params['tau'], params['m']))
        
        residuals['Residual'] = df['Mean'] - expected
        
        return residuals
    
    
    def build_residuals_1(self):
        """
        """
        nora_residuals = self._find_residuals_1(df=self.data['unrandomized'], params=self.func1_params['unrandomized']).reset_index()
        rand_residuals = self._find_residuals_1(df=self.data['randomized'], params=self.func1_params['randomized']).reset_index()
        
        self.func1_residuals = {'unrandomized':nora_residuals,
                                'randomized':rand_residuals}
    
    
    def do_fft(self):
        """
        """



class DataSetX:
    """
    Basically a fancy dictionary of calorimeter objects
    """
    
    def __init__(self, name = None, long_name = None, file = None, energy_range = [0, 99999999999]):
        self.name         = name
        self.long_name    = long_name
        self.file         = file
        self.energy_range = energy_range
        self.calos        = {}    # a dictionary of Calorimeter objects
        
    
    def set_up_calos(self, step_length, verbose=False):
        """
        """
        # Loop through all calos around the ring
        for caloNum in range(1, 25):
                        
            # Grab the histograms from the TFile
            nora_hist = self.file.Get("verticalPosition/xtals" + \
                                     str(caloNum)).Clone(self.name + "_Xcalo_" + str(caloNum))
            rand_hist = self.file.Get("verticalPosition/randxtals" + \
                                      str(caloNum)).Clone(self.name + "_rand_Xcalo_" + str(caloNum))
            
            # Set the energy axis range
            nora_hist.SetAxisRange(self.energy_range[0], self.energy_range[1], "y")
            rand_hist.SetAxisRange(self.energy_range[0], self.energy_range[1], "y")
            
            # Project to y vs time
            nora_proj = nora_hist.Project3D("zx")
            rand_proj = rand_hist.Project3D("zx")
            
            # Pack up into a dictionary
            this_calo_hists = {'randomized':rand_proj, 'unrandomized':nora_proj}
            
            # Construct the calo object
            this_calo = CalorimeterX(number=caloNum, step=step_length, hists=this_calo_hists)
            this_calo.build_dfs()
            
            # Put this calo into the dataset calo dictionary
            self.calos[caloNum] = this_calo
            
            if (verbose == True):
                print("Finished " + self.name + " calorimeter " + str(caloNum))
                
                
    def set_up_one_calo(self, step_length, caloNum, verbose=False):
        """
        """

        # Grab the histograms from the TFile
        nora_hist = self.file.Get("verticalPosition/xtals" + \
                                 str(caloNum)).Clone(self.name + "_Xcalo_" + str(caloNum))
        rand_hist = self.file.Get("verticalPosition/randxtals" + \
                                  str(caloNum)).Clone(self.name + "_rand_Xcalo_" + str(caloNum))

        # Set the energy axis range
        nora_hist.SetAxisRange(self.energy_range[0], self.energy_range[1], "y")
        rand_hist.SetAxisRange(self.energy_range[0], self.energy_range[1], "y")

        # Project to y vs time
        nora_proj = nora_hist.Project3D("zx")
        rand_proj = rand_hist.Project3D("zx")

        # Pack up into a dictionary
        this_calo_hists = {'randomized':rand_proj, 'unrandomized':nora_proj}

        # Construct the calo object
        this_calo = CalorimeterX(number=caloNum, step=step_length, hists=this_calo_hists)
        this_calo.build_dfs()

        # Put this calo into the dataset calo dictionary
        self.calos[caloNum] = this_calo

        if (verbose == True):
            print("Finished " + self.name + " calorimeter " + str(caloNum))
    
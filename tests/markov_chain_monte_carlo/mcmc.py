# MCMC simulation script for CHIME Catalog 1
# Run 'python [options] mcmc.py' in the terminal. Default number of CPU cores in parallel is 125.
import numpy as np
import pandas as pd
from scipy.stats import anderson_ksamp
import re
import time
import emcee
import corner
from copy import deepcopy
from joblib import Parallel, delayed
from tqdm import tqdm
import argparse

from frbpoppy import Survey, CosmicPopulation, SurveyPopulation, unpickle, pprint, merge_pop, paths
from frbpoppy.tns import get_chimefrb_catalog1

parser = argparse.ArgumentParser(description='Get output file serial Number.')
parser.add_argument('-d', '--delay_time', dest='delay_time', action='store',
                    type=float, default=0,
                    help='Provide the delay time (default 0)')
parser.add_argument('-n', '--number', dest='output_number', action='store',
                    type=str, default='1',
                    help='Provide the output file serial number (default 1)')
parser.add_argument('--nwalkers', dest='nwalkers', action='store',
                    type=int, default=60,
                    help='Provide the MCMC nwalkers (default 60)')
parser.add_argument('--nsteps', dest='nsteps', action='store',
                    type=int, default=500,
                    help='Provide the MCMC nsteps (default 500)')
parser.add_argument('--pop_size', dest='pop_size', action='store',
                    type=float, default=1e6,
                    help='Provide the pop size (default 1e6)')
parser.add_argument('--ncpus', '--ncores', dest='ncpus', action='store',
                    type=int, default=125,
                    help='Provide the number of CPU cores (default 125)')
parser.add_argument('--samp_size', '--sample_size', dest='samp_size', action='store',
                    type=int, default=1000,
                    help='Provide the redshift distribution model (default powerlaw)')
parser.add_argument('--zmodel', '--model', dest='zmodel', action='store',
                    type=str, default='vol_co',
                    help='Provide the redshift distribution model (default powerlaw)')
pargs = parser.parse_args()

chimefrbcat = pd.read_csv(paths.home + '/frbpoppy/data/frbcat/chimefrbcat1.csv', delimiter=',')
for i in range(len(chimefrbcat['scat_time'])):
    chimefrbcat['scat_time'][i] = str(chimefrbcat['scat_time'][i])
    chimefrbcat['scat_time'][i] = re.sub('<', '', chimefrbcat['scat_time'][i])
    chimefrbcat['scat_time'][i] = float(chimefrbcat['scat_time'][i])

def get_snr(pop, position):
    return pop.frbs.snr[position]

def get_w_eff(pop, position):
    return pop.frbs.w_eff[position]

def get_t_dm(pop, position):
    return pop.frbs.t_dm[position]

def get_dm(pop, position):
    return pop.frbs.dm[position]

def iter_pop(pop_size, survey, j, *args):
    
    pop = CosmicPopulation.complex(pop_size, mute=True)
    if pargs.zmodel=='vol_co':
        alpha, li, log_w_mean, w_std, dm_igm_slope, dm_host_mean, dm_host_std = locals()['args'][0]
        pop.set_dist(model='vol_co', z_max=1.5, alpha=alpha)
    elif pargs.zmodel=='sfr':
        li, log_w_mean, w_std, dm_igm_slope, dm_host_mean, dm_host_std = locals()['args'][0]
        pop.set_dist(model='sfr', z_max=1.5)
    elif pargs.zmodel=='delayed_sfr':
        li, log_w_mean, w_std, dm_igm_slope, dm_host_mean, dm_host_std = locals()['args'][0]
        pop.set_dist(model='delayed_sfr', z_max=1.5, delay_time=pargs.delay_time)
    pop.set_si(model='constant', value=-1.5)
    pop.set_lum(model='powerlaw', low=10**41.0, high=10**46.0, power=li)
    pop.set_w(model='lognormal', mean=10**log_w_mean, std=w_std)
    pop.set_dm_igm(model='ioka', slope=dm_igm_slope)
    pop.set_dm_host(model='lognormal', mean=dm_host_mean, std=dm_host_std)
    pop.generate()
    
    surv_pop = SurveyPopulation(pop, survey, scat=True, mute=True)
    
    i = 1
    npercore = round(pargs.samp_size/pargs.ncpus)
    while (len(surv_pop.frbs.z) < npercore) and (i < 5*npercore) and (len(surv_pop.frbs.z) > 5 if i > 25 else True):
        i += 1
        #pprint(len(surv_pop.frbs.z))
        pop.set_lum(model='powerlaw', low=10**41.0, high=10**46, power=li)
        pop.set_w(model='lognormal', mean=10**log_w_mean, std=w_std)
        pop.gen_lum()
        pop.gen_w()
        t_pop = SurveyPopulation(pop, survey, scat=True, mute=True)
        surv_pop = merge_pop(surv_pop, t_pop)
    return surv_pop

#survey_name = self.survey_names
#survey = Survey(name=survey_name)
#survey.set_beam(model=survey_name)
global survey
survey = Survey(name='chime-frb')
survey.set_beam(model='chime-frb')
        
class MCMC:

    def __init__(self):

        self.zmodel = pargs.zmodel
        self.nwalkers = pargs.nwalkers
        self.nsteps = pargs.nsteps
        self.pop_size = pargs.pop_size
        self.samp_size = pargs.samp_size
        
        self.survey_names = 'chime-frb'

        self.chime_catlog_1 = get_chimefrb_catalog1(repeater=False)

    def lnlike(self, theta):
        '''
        survey_name = self.survey_names
        survey = Survey(name=survey_name)
        survey.set_beam(model=survey_name)
        '''
        t_pop = Parallel(n_jobs=pargs.ncpus)(delayed(iter_pop)(self.pop_size, survey, j, theta) for j in tqdm(range(pargs.ncpus)))
        
        surv_pop = merge_pop(*t_pop)
        pop = surv_pop

        N = self.samp_size
        if len(pop.frbs.z) >= N:
            pprint(N)
            #snrN = get_snr(pop, range(N))
            #wN = get_w_eff(pop, range(N))
            #t_dmN = get_t_dm(pop, range(N))
            #dmN = get_dm(pop, range(N))
            try:
                snr_gof = -anderson_ksamp([pop.frbs.snr, self.chime_catlog_1.snr])[0]
                w_gof = -anderson_ksamp([np.round(pop.frbs.w_eff / 0.983) * 0.983, np.round(self.chime_catlog_1.w_eff / 0.983) * 0.983])[0]
                dm_gof = -anderson_ksamp([pop.frbs.dm, self.chime_catlog_1.dm])[0]
            except ValueError:
                snr_gof = -200
                w_gof = -200
                dm_gof = -200
        
        else:
            snr_gof = -200
            w_gof = -200
            dm_gof = -200
            pprint(len(pop.frbs.z))

        lp = snr_gof + w_gof + dm_gof 
        
        pprint(f'lp is {lp, snr_gof, w_gof, dm_gof}')
        
        return lp

    def lnprior(self, theta):

        if self.zmodel=='vol_co':
            alpha, li, log_w_mean, w_std, dm_igm_slope, dm_host_mean, dm_host_std = theta
            if (-2.5 < alpha < -0.5 and -2.0 < li < 0 and -1.0 < log_w_mean < 0.4 and 0 < w_std < 4 and 600 < dm_igm_slope < 1300 and 200 < dm_host_mean < 800 and 200 < dm_host_std < 900):
                return 0.0
            return -np.inf
        else:
            li, log_w_mean, w_std, dm_igm_slope, dm_host_mean, dm_host_std = theta
            if (-2.0 < li < 0 and -1.0 < log_w_mean < 0.4 and 0 < w_std < 4 and 600 < dm_igm_slope < 1300 and 200 < dm_host_mean < 800 and 200 < dm_host_std < 900):
                return 0.0
            return -np.inf

    def lnprob(self, theta):
    
        lp = self.lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.lnlike(theta)

if __name__ == '__main__':

    mcmc = MCMC()
    zmodel = mcmc.zmodel
    if zmodel=='vol_co':
        alpha_true = -1.4
        li_true = -1.5
        log_w_mean_true = -0.4
        w_std_true = 1.7
        dm_igm_slope_true = 900
        dm_host_mean_true = 570
        dm_host_std_true = 570
    else:
        li_true = -1.5
        log_w_mean_true = -0.4
        w_std_true = 1.4
        dm_igm_slope_true = 850
        dm_host_mean_true = 500
        dm_host_std_true = 500
    
    if zmodel=='vol_co':
        params = ['alpha', 'li', 'log_w_mean', 'w_std', 'dm_igm_slope', 'dm_host_mean', 'dm_host_std']
        init = np.array([alpha_true, li_true, log_w_mean_true, w_std_true, dm_igm_slope_true, dm_host_mean_true, dm_host_std_true])
        ndim = 7
        dflag = ''
    elif zmodel=='sfr':
        params = ['li', 'log_w_mean', 'w_std', 'dm_igm_slope', 'dm_host_mean', 'dm_host_std']
        init = np.array([li_true, log_w_mean_true, w_std_true, dm_igm_slope_true, dm_host_mean_true, dm_host_std_true])
        ndim = 6
        dflag = ''
    elif zmodel=='delayed_sfr':
        params = ['li', 'log_w_mean', 'w_std', 'dm_igm_slope', 'dm_host_mean', 'dm_host_std']
        init = np.array([li_true, log_w_mean_true, w_std_true, dm_igm_slope_true, dm_host_mean_true, dm_host_std_true])
        ndim = 6
        dflag = str(pargs.delay_time) + 'Gyr_'
    nwalkers = mcmc.nwalkers
    #p0 = [init + ([1e-3]*(ndim-3) + [1]*3)*np.random.randn(ndim) for i in range(nwalkers)]
    rng = np.random.default_rng()
    p0 = [init + ([1e-3]*(ndim-3) + [1]*3)*rng.normal(size=ndim) for i in range(nwalkers)]
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, mcmc.lnprob)
    
    nsteps = mcmc.nsteps
    sampler.run_mcmc(p0, nsteps, progress=True)

    np.savetxt('mcmc_sampler_' + zmodel + '_' + dflag + str(nwalkers*nsteps) + '_' + pargs.output_number + '.txt', sampler.flatchain, fmt='%5.8f')
    fig = corner.corner(sampler.flatchain, labels=params, dpi_plot=400, color="blue", smooth=2.0, show_titles=False)
    fig.savefig('mcmc_' + zmodel + '_' + dflag + str(nwalkers*nsteps) + '_' + pargs.output_number + '.pdf')

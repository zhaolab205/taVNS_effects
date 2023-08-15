
"""
@author: Siyu Zhu; Qi Liu; Weihua Zhao zarazhao@uestc.edu.cn
"""
import sys, os, pdb 
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns 
import hddm 
from joblib import Parallel, delayed
# import Parallel, delayed from IPython 
# import embed as shell 
import matplotlib.patches as mpatches

# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# import hddm from joblib 
# import Parallel, delayed

# def get_choice(row):
    
#     if row.condition == 'present':
#         if row.response == 1:
#             return 1
#         else:
#             return 0
#     elif row.condition == 'absent':
#         if row.response == 0:
#             return 1
#         else:
#             return 0
def fit_subject(data, quantiles):
    subj_idx = np.unique(data['subj_idx'])
    
    m = hddm.HDDMStimCoding(data, stim_col='response', split_param='v', bias=True, p_outlier=0.05,
            depends_on={ 'v':'condition'})
            
                            
    
    m.optimize('gsquare', quantiles=quantiles, n_runs=3) 
    
    # m = hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True, p_outlier=0,
    #                         depends_on={'v':'condition', 'a':'condition', 't':'condition', 'z':'condition', 'dc':'condition', })
    # m.optimize('gsquare', quantiles=quantiles, n_runs=8)
    
    res = pd.concat((pd.DataFrame([m.values], index=[subj_idx]), pd.DataFrame([m.bic_info], index=[subj_idx])), axis=1)
    
    return res

  

def summary_plot(df_group, df_sim_group=None, quantiles=[0, 0.1, 0.3, 0.5, 0.7, 0.9,], xlim=None, condition=None):
    
    total_nr_subjects = len(np.unique(df_group['subj_idx']))
    
    nr_subjects = 5
    
    fig = plt.figure(figsize=(16,nr_subjects*2))
    plt_nr = 1
    for ii,s in enumerate(np.unique(df_group['subj_idx'])):
        
        if ii == total_nr_subjects:
            break
        
        if ii%5 == 0: 
            if ii > 0:
                plt.suptitle('Subjects %d to %d, condition=%s'%(ii-5,ii-1,condition),fontsize=20)
                plt.savefig('fits_%d-%d_%s.png'%(ii-5,ii-1,condition), dpi=600)
                
            
            plt.show()
        
            fig = plt.figure(figsize=(16,nr_subjects*2))
            plt_nr=1
    
        print(s)
        plt_nr
    
        
    
        # Fetch subject s's entry from data
        df = df_group.copy().loc[(df_group['subj_idx']==s),:]
    
        # Fetch subject s's simulated data
        df_sim = df_sim_group.copy().loc[(df_sim_group['subj_idx']==s),:]
    
        # Incorrect are assigned neg RT (obs data)
        df['rt_acc'] = df['rt'].copy()
        df.loc[df['correct']==0, 'rt_acc'] *= -1
    
        # Incorrect are assigned neg RT (sim data)
        df_sim['rt_acc'] = df_sim['rt'].copy()
        df_sim.loc[df_sim['correct']==0, 'rt_acc'] = df_sim.loc[df_sim['correct']==0, 'rt_acc'] * -1
        #df_sim['stimulus'] = np.array((np.array(df_sim['response']==1) & np.array(df_sim['correct']==1)) + (np.array(df_sim['response']==0) & np.array(df_sim['correct']==0)), dtype=int)      

        # Set up RT bins based on observed data
        max_rt = np.percentile(df_sim.loc[~np.isnan(df_sim['rt']), 'rt'], 99)
        bins = np.linspace(-max_rt,max_rt,21)
    
        # Histogram of correct and incorrect observed RT
        ax = fig.add_subplot(nr_subjects,4,plt_nr)
        N, bins, patches = ax.hist(df.loc[:, 'rt_acc'], bins=bins,density=True, color='green', alpha=0.5)
    
        # Set patches of negative bins to red
        for bin_size, bin, patch in zip(N, bins, patches):
            if bin < 0:
                plt.setp(patch, 'facecolor', 'r')
                
            
        # Histogram of correct and incorrect from simulated RT
        if df_sim is not None:
            ax.hist(df_sim.loc[:, 'rt_acc'], bins=bins, density=True,histtype='step', color='k', alpha=1, label='Sim')
        
        red_patch = mpatches.Patch(color='red', label='Obs incorr', alpha=0.5)
        green_patch = mpatches.Patch(color='green', label='Obs corr', alpha=0.5)
        white_patch = mpatches.Patch(facecolor='none',edgecolor='black', label='Simulated')
        plt.legend(handles=[red_patch, green_patch, white_patch])
            
        
        ax.set_xlabel('RT(s)')
        ax.set_ylabel('Density')
        plt_nr += 1
        
        # condition accuracy plots:
        ax = fig.add_subplot(nr_subjects,4,plt_nr)
    
        df.loc[:,'rt_bin'] = pd.qcut(df['rt'], quantiles, labels=False, duplicates='drop')
        d = df.groupby(['rt_bin']).mean().reset_index()
        ax.errorbar(d.loc[:, "rt_bin"], d.loc[:, "correct"], fmt='-o', color='orange', markersize=10)
        ax.set_xticks(range(len(quantiles)))
        ax.set_xticklabels(map(str,quantiles))
    
        if df_sim is not None:
            df_sim.loc[:,'rt_bin'] = pd.qcut(df_sim['rt'], quantiles, labels=False, duplicates='drop')
            d = df_sim.groupby(['rt_bin']).mean().reset_index()
            ax.errorbar(d.loc[:, "rt_bin"], d.loc[:, "correct"], fmt='x', color='k', markersize=6)
            ax.set_xticks(range(len(quantiles)))
            ax.set_xticklabels(map(str,quantiles))
            
        
        
        plt.legend(['Obs','Sim'])
    
        if xlim:
            ax.set_xlim(xlim)
        ax.set_ylim(0, 1.2)  
        ax.set_xlabel('RT quantiles')
        ax.set_ylabel('Fraction correct')
        plt_nr += 1
    
        # Fetch subject s's entry from data 
        df = df_group.copy().loc[(df_group['subj_idx']==s),:]
    
        # Fetch subject s's simulated data
        df_sim = df_sim_group.copy().loc[(df_sim_group['subj_idx']==s),:]
    
        df['rt_resp'] = df['rt'].copy()
        df.loc[df['condition']==0, 'rt_resp'] *= -1 # Nogo RT set to neg (obs data)
    
        df_sim['rt_resp'] = df_sim['rt'].copy()
        df_sim.loc[df_sim['condition']==0, 'rt_resp'] *= -1 # Nogo RT set to neg (sim data)
    
        # rt distributions : 
        ax = fig.add_subplot(nr_subjects,4,plt_nr)
        N, bins, patches = ax.hist(df.loc[:, 'rt_resp'], bins=bins, density=True, color='cyan', alpha=0.5)
    
        for bin_size, bin, patch in zip(N, bins, patches):
            if bin < 0:
                plt.setp(patch, 'facecolor', 'm')
           
            
        if df_sim is not None:
            ax.hist(df_sim.loc[:,'rt_resp'], bins=bins, density=True, 
                 histtype='step', color='k', alpha=1, label=None)
       
        
        ax.set_xlabel('RT(s)')
        ax.set_ylabel('Density')
        plt_nr += 1
    
        mag_patch = mpatches.Patch(color='magenta', label='Obs no-go resp', alpha=0.5)
        cyan_patch = mpatches.Patch(color='cyan', label='Obs go resp', alpha=0.5)
        white_patch = mpatches.Patch(facecolor='none',edgecolor='black', label='Simulated')
        plt.legend(handles=[mag_patch, cyan_patch, white_patch])
    
        # condition accuracy plots:
        ax = fig.add_subplot(nr_subjects,4,plt_nr)
    
        df.loc[:,'rt_resp'] = pd.qcut(df['rt'], quantiles, labels=False, duplicates='drop')
        d = df.groupby(['rt_resp']).mean().reset_index()
        ax.errorbar(d.loc[:, "rt_resp"], d.loc[:, "response"], fmt='-o', color='orange', markersize=10)
        ax.set_xticks(range(len(quantiles)))
        ax.set_xticklabels(map(str,quantiles))
    
        if df_sim is not None:
            df_sim.loc[:,'rt_resp'] = pd.qcut(df_sim['rt'], quantiles, labels=False, duplicates='drop')
            d = df_sim.groupby(['rt_resp']).mean().reset_index()
            ax.errorbar(d.loc[:, "rt_resp"], d.loc[:, "response"], fmt='x', color='k', markersize=6)
            ax.set_xticks(range(len(quantiles)))
            ax.set_xticklabels(map(str,quantiles))
            
        
        
        plt.legend(['Obs','Sim'])
    
        if xlim:
             ax.set_xlim(xlim)
        ax.set_ylim(0, 1.2)
        ax.set_xlabel('RT quantiles')
        ax.set_ylabel('Fraction go')
        plt_nr += 1
       
    
 
    sns.despine(offset=3, trim=True)

    plt.suptitle('Subjects %d to %d, condition=%s'%(ii-ii%5,ii,condition),fontsize=20)
    plt.savefig('fits_%d-%d_%s.png'%(ii-ii%5,ii,condition), dpi=600)

    plt.show()

    return fig   
    

# Fit G-square:   
dataset = hddm.load_csv('hddm_taVNS_gsquare.csv')
n_subjects = len(np.unique(dataset.subj_idx))
# dataset.insert(loc=5, column='correct', value=(dataset.response == dataset.condition).astype(int))

quantiles = [.1, .2, .3, .4, .5, .6, .7, .8, .9]

params_fitted = pd.concat(Parallel(n_jobs=n_subjects)(delayed(fit_subject)(data[1], quantiles)
                                     for data in dataset.groupby('subj_idx')))
params_fitted.drop(['bic', 'likelihood', 'penalty'], axis=1, inplace=True)

# n_subjects = len(np.unique(df_emp.subj_idx))
# quantiles = [.1, .3, .5, .7, .9]
# params_fitted = pd.concat(Parallel(n_jobs=n_subjects)(delayed(fit_subject)(data[1], quantiles)
#                                                       for data in df_emp.groupby('subj_idx')))
# if not os.path.exists('params_fitted_5.pkl'):
#     params_fitted = pd.concat(Parallel(n_jobs=n_subjects)(delayed(fit_subject)(data[1], quantiles)
#                                     for data in dataset.groupby('subj_idx')))
#     pd.DataFrame.to_pickle(params_fitted, 'params_fitted_5.pkl')

# else:
#     params_fitted = pd.read_pickle('params_fitted_5.pkl')
    
# params_fitted.to_csv('params_fitted_5.csv')


def simulate_data(subj, params_fitted, nr_trials=1000,rate_0_to_1=1.0):
    
    # Fetch Treatment-go parameters from params_fitted
    go_index = ['a','v(0)','t','z']
    params_go = {ii[0]:params_fitted.loc[subj,ii].values for ii in go_index}
    
    # Generate Treatment-go simulated data
    df_Treatment_go,_ = hddm.generate.gen_rand_data(params=params_go, size=nr_trials, subjs=1, subj_noise=0)
    df_Treatment_go['Treatment'] = 'taVNS'
    df_Treatment_go['condition'] = 1
    
    # Fetch safe-no go parameters from params_fitted
    nogo_index = ['a','v(1)','t','z']
    params_nogo = {ii[0]:params_fitted.loc[subj,ii].values for ii in nogo_index}
    
    # params_nogo['sv'] = params_fitted.loc[subj].sv
    # params_go['sv'] = params_fitted.loc[subj].sv
    # params_nogo['sz'] = params_fitted.loc[subj].sz
    # params_go['sz'] = params_fitted.loc[subj].sz
    # params_nogo['st'] = params_fitted.loc[subj].st
    # params_go['st'] = params_fitted.loc[subj].st
    
    # Generate safe no-go simulated data
    df_Treatment_nogo,_ = hddm.generate.gen_rand_data(params=params_nogo, size=nr_trials*rate_0_to_1, subjs=1, subj_noise=0)
    df_Treatment_nogo['Treatment'] = 'taVNS'
    df_Treatment_nogo['condition'] = 0
    
    # # Fetch threat-go parameters from params_fitted
    # params_TH1 = {ii[0]:params_fitted.loc[subj,ii].values for ii in params_fitted.keys() if 'TH.1' in ii or not '.0' in ii and not 'SA' in ii}
    
    # # Generate threat-go simulated data
    # df_threat_go,_ = hddm.generate.gen_rand_data(params=params_TH1, size=nr_trials, subjs=1, subj_noise=0)
    # df_threat_go['block'] = 'TH'
    # df_threat_go['condition'] = 1
    
    # # Fetch threat no-go parameters from params_fitted
    # params_TH0 = {ii[0]:params_fitted.loc[subj,ii].values for ii in params_fitted.keys() if 'TH.0' in ii or not '.1' in ii and not 'SA' in ii}
    # params_TH0['v']=-params_TH0['v']
    
    # params_TH0['sv'] = params_fitted.loc[subj].sv
    # params_TH1['sv'] = params_fitted.loc[subj].sv
    # params_TH0['sz'] = params_fitted.loc[subj].sz
    # params_TH1['sz'] = params_fitted.loc[subj].sz
    # params_TH0['st'] = params_fitted.loc[subj].st
    # params_TH1['st'] = params_fitted.loc[subj].st
    
    # # Generate threat no-go simulated data
    # df_threat_nogo,_ = hddm.generate.gen_rand_data(params=params_TH0, size=nr_trials*rate_0_to_1, subjs=1, subj_noise=0)
    # df_threat_nogo['block'] = 'TH'
    # df_threat_nogo['condition'] = 0
    
    # Dataframe with simulated data
    df_sim = pd.concat((df_Treatment_go, df_Treatment_nogo))
    
    # Add column with correct (0 or 1)
    df_sim.insert(loc=5, column='correct', value=(df_sim.response == df_sim.condition).astype(int))
    
    
    return df_sim

# simulate data based on fitted params:
dfs = []
sub_id = [0]*n_subjects
for i in range(n_subjects):
    sub_id[i] = dataset.loc[i*192,'subj_idx']
for i in sub_id:
    
    rate_0_to_1 = len(dataset.loc[(dataset.subj_idx==i)
    &(dataset.condition==0)])/len(dataset.loc[(dataset.subj_idx==i)])
    
    df = simulate_data(i,params_fitted,nr_trials=20000,rate_0_to_1=rate_0_to_1)                           
    df['subj_idx'] = i
    dfs.append(df)
    
    
df_sim = pd.concat(dfs)

# # plot data with model fit on top:
c_id = [0]*2
# for c in np.unique(dataset['condition']):
# print('BLOCK {}'.format(c))
    
summary_plot(df_group=dataset,
          df_sim_group=df_sim, quantiles=quantiles,
          condition='taVNS' )
    
plt.show() 
pdb.set_trace()    
      

    

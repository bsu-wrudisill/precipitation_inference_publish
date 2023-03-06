import numpy as np
import sys
import seaborn as sb
import scipy.stats as stats
import matplotlib.pyplot as plt
import pickle
import argparse
import pymc3 as pm
import pathlib
# import some cool fortran modules that are precompiled and fast
from model_paths import fortran_path, base_path, data_path
sys.path.append(fortran_path)
sys.path.append(base_path)
from inference import *
from ForcingModel import ForcingModel
from model_wrapper_v2 import *
import argparse
import copy

def kge(mod, obs):
    # mean ratio
    b = np.mean(mod) / np.mean(obs)
    # std
    a = np.std(mod) / np.std(obs)
    # corr coeff
    r = np.corrcoef(mod, obs)[0, 1]  # corrcoef returns the correlation matrix...
    # the diagonals are 1, the off-diags are the 'r'
    # value that we want
    kgeval = 1 - np.sqrt((r - 1.)**2 + (a - 1.)**2 + (b - 1)**2)
    return kgeval

# -------------------------------- #
parser = argparse.ArgumentParser()
parser.add_argument("trace", type=str)
parser.add_argument("prior", type=str)
parser.add_argument("hydro_file", type=str)
parser.add_argument("plot_name", type=str)
parser.add_argument("--param_file", type=str, help="parameter (yml) file", default="./prior_parameters/param.yml")
parser.add_argument("--silent", type=bool, help="yes/no to log output", default=True)
parser.add_argument("--output_loc", type=str, help="output location directory", default="./")
args = parser.parse_args()


# get important stuff

trace = pathlib.Path(args.trace)
prior = pathlib.Path(args.prior)
print(trace.name)
start_date = trace.name.split("_")[1]
end_date = trace.name.split("_")[2].split(".")[0]
print(start_date, end_date)

# open up the trace and priors
with open(prior, 'rb') as buff:
    prior_dict = pickle.load(buff)

with open(trace, 'rb') as buff:
    trace = pickle.load(buff)
#        pm.traceplot(trace)
tp = trace.points()

# configure model
nlayers = 100
model = model_wrapper()
model.read_parameter_files(param_file=hydro_file)
model.create_forcings(start_date, end_date, nlayers)

# open the parameter dictionaries...
with open(args.param_file, "r") as p:
    dct = yaml.load(p, yaml.FullLoader)


# log the initial condition values...
model.snow_model_args.update(rvs=1, opg_method=2)
model.hydro_model_args.update(dict(sm1i=.1, sm2i=.1, hydroop=0))

# RUN MODEL
kgebest = -100.
lkbest = -1000

# now plot the posterior
kge_list = []

# for i,tr in enumerate(tp):
#     update_snow  = {}
#     for key in dct['snow'].keys():
#         update_snow[key] = tr.get(key)

#     update_hydro  = {}
#     for key in dct['hydro'].keys():
#         update_hydro[key] = tr.get(key)

#     # update the model..
#     model.snow_model_args.update(update_snow)
#     model.hydro_model_args.update(update_hydro)

#     # run it.
#     qchan = model.runmain([])['qchan']
#     kge_i = kge(qchan, model.q_obs.values)
#     kge_list.append(kge_i)
#     #plt.plot(qchan, color='blue', alpha=.05)

#     if kge_i > kgebest:
#         kgebest = kge_i
#         trbest = tr


#kge_list = np.array(kge_list)
# plot the best value


def avg_from_trace(trace, dct):
    # return copy of the model object
    # update_snow  = {}
    # for key in dct['snow'].keys():
    #     update_snow[key] = trace[key].mean()

    update_hydro  = {}
    for key in dct['snow'].keys():
        update_hydro[key] = trace[key].mean()

    # get the sigma...
    # update the model and return it
    # model.snow_model_args.update(update_snow)
    #model.hydro_model_args.update(update_hydro)
    #return copy.deepcopy(model)
    return update_hydro, trace["sigma"].mean(), trace["alpha"].mean()

def std_from_trace(trace, dct):
    update_hydro  = {}
    for key in dct['snow'].keys():
        update_hydro[key] = trace[key].std()

    # update the model and return it
    # model.snow_model_args.update(update_snow)
    #return copy.deepcopy(model)
    #model.hydro_model_args.update(update_hydro)
    return update_hydro

def get_whole_trace(trace, dct, model, n_samples, t):
    # return copy of the model object
    tpl = list(trace.points())
    l = len(tpl)
    #Vqcarr = np.empty((l,t))
    qcarr = np.empty((n_samples,t))

    random_selection = np.random.randint(0,l-1,n_samples)
    for i,r in enumerate(random_selection):

        # get the random ppint...
        tp = tpl[r]

        # run the model...
        # update_snow  = {}
        # for key in dct['snow'].keys():
        #     update_snow[key] = tp[key]
        update_hydro  = {}
        for key in dct['snow'].keys():
            update_hydro[key] = tp[key]

        # update the model and return it
        # model.snow_model_args.update(update_snow)
        #model.hydro_model_args.update(update_hydro)
        #qc=model.runmain([])['qchan']
        qc = model.runmain(list(update_hydro.values()))["qchan"]
        #ax.plot(qc, alpha=.1, color='blue')
        qcarr[i,:] = qc
    return qcarr

        #return copy.deepcopy(model)

def plot_whole_prior(prior_dct, dct, model, ax):
    len_of_dict = len(prior_dict[[k for k in prior_dict][0]])
    for i in range(len_of_dict):
        # update_snow = {}
        # for key in dct['snow'].keys():
        #     update_snow[key] = prior_dict.get(key)[i]

        update_hydro  = {}
        for key in dct['snow'].keys():
            update_hydro[key] = prior_dict.get(key)[i]

        # model.snow_model_args.update(update_snow)
        model.hydro_model_args.update(update_hydro)
        qc = model.runmain([])['qchan']

        plt.plot(qc, color='green')

# model_best = run_from_trace(trbest, dct, model)
# qchan_best = model_best.runmain([])['qchan']

### ############## PLOTTING BELOW HERE ################

# GET THE MEAN AND STD OF THE TRAC PARAMETERS
inferred_params_mean_dct, sigma_mean, alpha_mean  = avg_from_trace(trace, dct)
inferred_params_std_dct =  std_from_trace(trace, dct)
model.substitute_arguments(list(inferred_params_mean_dct.keys()), [])

inferred_params_mean = list(inferred_params_mean_dct.values())
inferred_params_std  = list(inferred_params_std_dct.values())

qchan_avg =   model.runmain(inferred_params_mean)['qchan']
qchan_1std  = model.runmain([(inferred_params_mean[i]+x) for i,x in enumerate(inferred_params_std)])['qchan']
qchan_m1std = model.runmain([(inferred_params_mean[i]-x) for i,x in enumerate(inferred_params_std)])['qchan']


### GET THE ENTIRE TRACE STREAMFLOW ######
n = 2000
qc = get_whole_trace(trace, dct, model, n, len(qchan_avg))
fig,ax = plt.subplots()
fig.set_size_inches(12,4)
x = model.daily_swe.index

#LOOP through and plot n of the trace lines...
best_kge = 0. 
best_kge_i = 0
for i in range(n):
    ax.plot(x, qc[i,:], color='grey', alpha=.2, linewidth=.1)
    kgeval = kge(qc[i,:], model.q_obs.values)
    if kgeval > best_kge:
            best_kge = kgeval
            best_kge_i = i

# Plot the average
ax.plot(x, qchan_avg, label='Posterior', linewidth=.8)
ax.plot(x, qc[best_kge_i, :], label='Best-KGE', linewidth=.8, color='green')
conf = alpha_mean + qchan_avg*sigma_mean
ax.fill_between(x, qchan_avg + conf, qchan_avg-conf, color='blue', alpha=.2)

# Plot +/- 1 std
# ax.plot(x, qchan_1std,   linewidth=.4, color='black')
# ax.plot(x, qchan_m1std,  linewidth=.4, color='black')

# # NOW PLOT
ax.plot(x, model.q_obs.values, color='red', label='data', linewidth=.6)
ax.set_ylim(0,14)
ax.set_ylabel("mm")
plt.legend()
print(kge(model.q_obs.values, qchan_avg))
print(best_kge)
plt.savefig(args.plot_name, dpi=400)
#plt.show()



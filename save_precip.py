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


def kge(mod, obs):
    # mean ratio
    b = np.mean(mod) / np.mean(obs)
    # std
    a = np.std(mod) / np.std(mod)
    # corr coeff
    r = np.corrcoef(mod, obs)[0, 1]  # corrcoef returns the correlation matrix...
    # the diagonals are 1, the off-diags are the 'r'
    # value that we want
    kgeval = 1 - np.sqrt((r - 1.)**2 + (a - 1.)**2 + (b - 1)**2)
    return kgeval

# End helper functions


parser = argparse.ArgumentParser()
parser.add_argument("trace", type=str)
parser.add_argument("prior", type=str)
parser.add_argument("--param_file", type=str, help="parameter (yml) file", default="./prior_parameters/param.yml")
parser.add_argument("--silent", type=bool, help="yes/no to log output", default=True)
parser.add_argument("--output_loc", type=str, help="output location directory", default="./")
args = parser.parse_args()


# get important stuff
trace = pathlib.Path(args.trace)
prior = pathlib.Path(args.prior)
start_date = trace.name.split("_")[1]
end_date = trace.name.split("_")[2].split(".")[0]
print(start_date, end_date)

# open up the trace and priors
with open(prior, 'rb') as buff:
    prior_dict = pickle.load(buff)

with open(trace, 'rb') as buff:
    trace = pickle.load(buff)
#        pm.traceplot(trace)

# get all of the trace points...

# configure model
nlayers = 100
model = model_wrapper()
model.read_parameter_files()
model.create_forcings(start_date, end_date, nlayers)
model.substitute_arguments([],[])

# open the parameter dictionaries...
with open(args.param_file, "r") as p:
    dct = yaml.load(p, yaml.FullLoader)

# update initial conditions
model.hydro_model_args.update(dct['initial_conditions'])
elev_areas = model.elevation_bins_areas.reshape(nlayers,1)


# DO the prior precipitation
def p_from_dict(trace_dict, dct, model):

    precip_list = []
    len_of_dict = len(trace_dict[[k for k in trace_dict][0]])
    for i in range(len_of_dict):
        update_snow = {}
        for key in dct['snow'].keys():
            update_snow[key] = trace_dict.get(key)[i]

        update_hydro  = {}
        for key in dct['hydro'].keys():
            update_hydro[key] = trace_dict.get(key)[i]

        # now run the model..
        model.snow_model_args.update(update_snow)
        model.hydro_model_args.update(update_hydro)
        precip_list.append(model.runmain([])['qin'].sum())
    return precip_list


inferred_precip_list = p_from_dict(trace, dct, model)
prior_precip_list = p_from_dict(prior_dict, dct, model)

# write out the precipitation
inf_p_name="inferred_precip_%s-%s.pkl"%(start_date, end_date)
p_p_name="prior_precip_%s-%s.pkl"%(start_date, end_date)

# write the files..
with open(inf_p_name, "wb") as handle:
        pickle.dump(inferred_precip_list, handle)

with open(p_p_name, "wb") as handle:
        pickle.dump(prior_precip_list, handle)




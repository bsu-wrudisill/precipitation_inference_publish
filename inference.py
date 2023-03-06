"""Summary
"""
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas as pd
import pickle
import sys
import argparse
import pymc3 as pm
import theano
import theano.tensor as tt
#from model_wrapper import *
from model_wrapper_v2 import model_wrapper
import time
import yaml
import pathlib


def sls_like(func, parameters, data):
    model = func(parameters[:-2])['qchan']
    sigma0 = parameters[-1]
    alpha = parameters[-2]
    sigma = alpha + model*sigma0
    return -1*np.sum(0.5 * np.log(2 * np.pi * sigma ** 2) + ((data - model) / sigma) ** 2)
    #return kgeval




# weighted least squares like
def wls_like(func, parameters, data):
    model = func(parameters[:-2])['qchan']
    alpha = parameters[-1]
    beta = parameters[-2]
    warm_up = 0
    model = model[warm_up:]
    data = data[warm_up:]

    sigma = data*alpha + beta
    # return (-.5 / sigma**2)*(np.sum(qpred - data)**2)
    # return -np.sum(0.5 * np.log(2 * np.pi * sigma ** 2) + ((data - model) / sigma) ** 2)
    return -np.sum(0.5 * np.log(2 * np.pi * sigma ** 2) + ((data - model) / sigma) ** 2)


class likelihood_evaluate(tt.Op):
    itypes = [tt.dvector] # expects a vector of parameter values when called
    otypes = [tt.dscalar] # outputs a single scalar value (the log likelihood)

    def __init__(self, func, loglike, obs_data):
        # add inputs as class attributes
        self.func = func
        self.likelihood = loglike
        self.obs_data = obs_data

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        parameters, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(self.func, parameters, self.obs_data)
        outputs[0][0] = np.array(logl)
        #print(outputs)


def main(func,
         dict_of_parameters,
         obs_q,
         tracename="trace.pkl",
         priorname="prior.pkl"):

    # create the logliklihood function
    logl = likelihood_evaluate(func, sls_like, obs_q)

    #logl = likelihood_evaluate(func, wls_like, obs_q)
    with pm.Model() as model:
        # This list collects the model paramters that are passed in by "dict of parameters"
        # this is just a convenience
        list_of_parameters = []
        for key in dict_of_parameters.keys():
            p0 = dict_of_parameters.get(key)
            if p0[0] == "normal":
                p1 = pm.Normal(key, mu=p0[1], sigma=p0[2])

            if p0[0] == "uniform":
                p1 = pm.Uniform(key, lower=p0[1], upper=p0[2])#, transform=None)

            # append it to the list ..
            list_of_parameters.append(p1)

        # Append the "sigma" parameter to the end of the list. this is the sigma
        # in the error model
        list_of_parameters.append(pm.Normal("alpha2", 0, .15))
        list_of_parameters.append(pm.HalfNormal("sigma", .5))
        #list_of_parameters.append(pm.Uniform("sigma", lower=.1, upper=1.0))
        parameters = tt.as_tensor_variable(list_of_parameters)
        pmpot = pm.Potential('loglike', logl(parameters))

        prior = pm.sample_prior_predictive(samples=1000)

        step1 = pm.DEMetropolis()
        trace = pm.sample(draws=15000, step = step1, chains=20, tune=5000, discard_tuned_samples=True, cores=28)
        #trace = pm.sample(draws=5000, step = step1, chains=8, tune=1000, discard_tuned_samples=True, cores=10)

        with open(tracename, 'wb') as buff:
            pickle.dump(trace, buff)

        with open(priorname, 'wb') as buff:
             pickle.dump(prior, buff)

        pm.traceplot(trace)
        plt.savefig("most_recent_trace_long")
        print(pm.summary(trace))




if __name__ == '__main__':
    __spec__ = None


    print("pymc3 version: ", pm.__version__)

    # Set up the argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("start_date", type=str, help="WaterYear (YYYY)")
    parser.add_argument("end_date", type=str, help="WaterYear (YYYY)")
    parser.add_argument("--param_file", type=str, help="parameter (yml) file", default="./prior_parameters/param.yml")
    parser.add_argument("--silent", type=bool, help="yes/no to log output", default=True)
    parser.add_argument("--output_loc", type=str, help="output location directory", default="./")
    args = parser.parse_args()


    # set up the logger
    if not args.silent:
        suffix = datetime.datetime.now().strftime("%m-%d_%H%M%S")
        logfile= '{}{}log_{}.log'.format(args.wateryear, month_double_pad, suffix)
        file_handler = logging.FileHandler(filename=logfile)
        stdout_handler = logging.StreamHandler(sys.stdout)
        logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)15s %(levelname)-8s %(message)s',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    handlers=[file_handler, stdout_handler]
                    )
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)

    else:
        print("quiet")

    # print or log some stuff at the top of the script
    open_msg = "BEGIN PRECIPITATION INFERENCE\n Start Date: %s \n End Date: %s"%(args.start_date, args.end_date)
    logger.info(open_msg) if not args.silent else print(open_msg)

    # setup the trace and prior output names
    tracename = "%strace_%s_%s.pkl"%(args.output_loc, args.start_date, args.end_date)
    priorname = "%sprior_%s_%s.pkl"%(args.output_loc, args.start_date, args.end_date)
    savemsg = "saving outputs to...\n %s \n %s"%(tracename, priorname)
    logger.info(savemsg) if not args.silent else print(savemsg)

    # set up the forward model
    nlayers = 100
    model = model_wrapper()
    model.read_parameter_files()

    # make the spinup period...
    model.create_forcings(args.start_date, args.end_date, nlayers)

    # open the parameter dictionaries...
    with open(args.param_file, "r") as p:
        dct = yaml.load(p, yaml.FullLoader)

    # log the initial condition values...
    model.snow_model_args.update(rvs=1, opg_method=2)
    model.hydro_model_args.update(dict(sm1i=.1, sm2i=.1, hydroop=0, nr=3., sm1max=400.))

    # combine the args that are read in from the bayesian parameter dict
    if dct['snow'] is not None:
        joined_dct = dct['snow'] | dct['hydro']

    else:
        joined_dct = dct["hydro"]


    # parse a message to print out
    options_msg = "PRIOR PARAMETERS \n"
    options_msg = options_msg + "File: %s \n"%args.param_file
   # options_msg = options_msg + "\n".join(["%s %s"%(k,v) for i in ['snow', 'hydro'] for k,v in dct[i].items()])
    options_msg = options_msg + "\n".join(["%s %s"%(k,v) for k,v in joined_dct.items()])

    # optionally log or print the parameters....
    logger.info(options_msg) if not args.silent else print(options_msg)

    # parse the output and pass into model
    #model.substitute_arguments(dct['snow'].keys(), dct['hydro'].keys())
    model.substitute_arguments([], dct['hydro'].keys())

    # join the dictionaries

    # model
    model.Verbose = True

    # join the model args into one dict
    main(model.runmain, joined_dct, model.q_obs.values,
                   tracename=tracename,
                   priorname=priorname)








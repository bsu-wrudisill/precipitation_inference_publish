import numpy as np
import sys
import pickle
import pandas as pd
import matplotlib.pyplot as plt

# supress numpy/base python conflict about boolean comparisions of strings an np vectors
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


from model_paths import fortran_path, base_path, data_path
sys.path.append(fortran_path)
sys.path.append(base_path)
from driver import driver as dr
from snowmodule17 import snowmodule17 as sm
from ForcingModel import ForcingModel
import yaml



class model_wrapper:
    # schofield snotel station...

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.hydro_model_args  = {"hydroop" : None,
                                  "ntimes" : None,
                                  "pet" : None,
                                  "tair" : None,
                                  "qin" : None,
                                  "sm1i" : None,
                                  "sm2i" : None,
                                  "sm1max" : None,
                                  "sm2max" : None,
                                  "fracten": None,
                                  "ku" : None,
                                  "c" : None,
                                  "ks" : None,
                                  "lowercasen" : None,
                                  "beta" : None,
                                  "nr" : 3.2,
                                  "kr" : .5}

        self.snow_model_args   = {"ntimes" : None,
                                  "jdayvec" : None,
                                  "precipvec" : None,
                                  "tairvec" : None,
                                  "nlayers" : None,
                                  "dz" : None,
                                  "dt" : 24,
                                  "rvs" : 1,
                                  "opg_method": 1,
                                  "pm" : 0,
                                  "opg" : None,
                                  "bias" : None,
                                  "uadj" : None,
                                  "mbase" : None,
                                  "mfmax" : None,
                                  "mfmin" : None,
                                  "tipm" : None,
                                  "nmf" : None,
                                  "plwhc" : None,
                                  "pxtemp" : None,
                                  "pxtemp1" : None,
                                  "pxtemp2" : None,
                                  "alpha" : None,
                                  "ecutoff" : None,
                                  "t_lapse" : None,
                                  "t_bias" : None}

        self.snow_sub_args_len = 0
        self.hydro_sub_args_len = 0
        self.other_sub_args_len = 0
        self.arg_order = []

    def create_forcings(self, start_date, end_date, nlayers):
        frc = ForcingModel()
        frc.read()
        frc.PrepareDEM()
        snotel_elevation = 3096.768

        # read in forcings for the desired time period
        daily_temp, daily_precip, q_obs, day_len_hrs, daily_swe = frc(start_date, end_date)
        self.q_obs = q_obs
        self.daily_swe = daily_swe
        self.daily_precip = daily_precip
        self.daily_temp = daily_temp

        ### Do the elevation computing
        elevation_bins = np.array_split(frc.hypsometry, nlayers)
        elevation_means = np.array([x.mean() for x in elevation_bins])
        elevation_bins_areas = np.array([len(x)/len(frc.hypsometry) for x in elevation_bins])

        self.elevation_bins = elevation_bins
        self.elevation_bins_areas = elevation_bins_areas
        self.elevation_means = elevation_means
        dz = elevation_means - snotel_elevation

        # process dateintrs
        dates = pd.date_range(start_date, end_date, freq='D')
        jdays = dates.dayofyear.values * 1.0
        ntimes=len(dates)

        # create PET forcing
        PET = np.zeros_like(jdays)
        i=0
        for l,t in zip(day_len_hrs, daily_temp):
           PET[i] = 1.0*dr.pet(l,t)
           i+=1

        # Now update the model arguments
        self.forcing_dictionary = {"ntimes":ntimes,
                                   "tair": daily_temp.values,
                                   "dz": dz,
                                   "tairvec" : daily_temp.values,
                                   "precipvec": daily_precip.values,
                                   "ntimes" : ntimes,
                                   "jdayvec"  : jdays,
                                   "nlayers" : nlayers,
                                   "pet" : PET}

        # UPDATE THE MODEL ARGUMENTS
        for k,v in self.forcing_dictionary.items():
            # check if the key is in the dict we want..
            # we don't want to add new things to the
            # respective paramter dictionaries
            if k in self.hydro_model_args.keys():
                self.hydro_model_args.update({k:v})

            if k in self.snow_model_args.keys():
                self.snow_model_args.update({k:v})

    def read_parameter_files(self, param_file="basic_parameters.yml"):
        if self.verbose:
            print("Reading saved parameter files...")

        # read in some parameters
        param_base = base_path + "/" + param_file    #"/parameter_files/snow17params.pkl"

        with open(param_base, "r") as p:
            dct = yaml.load(p, yaml.FullLoader)

        # update self with these
        self.snow_model_args.update(dct['snow'])
        self.hydro_model_args.update(dct['hydro'])

        # let the user know what happened
        if self.verbose:
            [print(k, " : ", v) for k,v in dct['snow'].items()]
            [print(k, " : ", v) for k,v in dct['hydro'].items()]

#    def update_run_args_general(self, update_args):


    def substitute_arguments(self,
                             snow_sub_args = [],
                             hydro_sub_args = [],
                             other_sub_args = []):
        # Use this function to set up the run parameters
        # dictionary with 'xs'. This will tell the next
        # function that these values are to be substituted
        # with vectorized values (x[0], x[1],.... etc).
        arg_order = []
        for key in snow_sub_args:
            if key in self.snow_model_args.keys():
                self.snow_model_args.update({key:'x'})
                arg_order.append(key)

        for key in hydro_sub_args:
            if key in self.hydro_model_args.keys():
                self.hydro_model_args.update({key:'x'})
                arg_order.append(key)

        self.snow_sub_args_len  = len(snow_sub_args) if snow_sub_args != None else None
        self.hydro_sub_args_len = len(hydro_sub_args) if snow_sub_args != None else None

        # now print out the order of the args...
        self.arg_order = arg_order
        print("Argument order: ", self.arg_order)

########## NOT WORKING ###########
    def run_snow_only(self, x):
        snow_model_args_copy = self.snow_model_args.copy()
        # pmult = snow_model_args_copy.get("pm")
        # snow_model_args_copy.update({"precipvec": snow_model_args_copy.get("precipvec")*pmult})

        if len(x) != len(self.arg_order):
            print(len(x), len(self.arg_order))
            print("args not correct len")
            print(x)
            print(print(self.arg_order))
            return

        snow_model_args_copy = self.snow_model_args.copy()
        hydro_model_args_copy = self.hydro_model_args.copy()

        for arg,xi in zip(self.arg_order, x):
            if arg in hydro_model_args_copy.keys():
                print(arg, '--> ', xi) if self.verbose else None
                hydro_model_args_copy.update({arg:xi})

            elif arg in snow_model_args_copy.keys():
                print(arg, '--> ', xi) if self.verbose else None
                snow_model_args_copy.update({arg:xi})
            else:
                if self.verbose:
                    print("error, arg not found ")
                    print(arg)
                return


        _snow = sm.snow17driver(**snow_model_args_copy)

        snow_out  = {"outflow":_snow[0],
                     "swe":_snow[1],
                     "rain_total":_snow[2],
                     "precip_total":_snow[3],
                     "tai_layers":_snow[4]}
        return snow_out

########## NOT WORKING ###########



    def runmain(self,x):
        # This can only be a function of X and self. other arguments
        # are not permitted which is why this looks weird
        # Run the snow and the hydrologic model

        snow_model_args_copy = self.snow_model_args.copy()
        hydro_model_args_copy = self.hydro_model_args.copy()

        # loop through the
        if len(x) != len(self.arg_order):
            print(len(x), len(self.arg_order))
            print("args not correct len")
            print(x)
            print(print(self.arg_order))
            return

        for arg,xi in zip(self.arg_order, x):
            if arg in hydro_model_args_copy.keys():
                print(arg, '--> ', xi) if self.verbose else None
                hydro_model_args_copy.update({arg:xi})

            elif arg in snow_model_args_copy.keys():
                print(arg, '--> ', xi) if self.verbose else None
                snow_model_args_copy.update({arg:xi})
            else:
                if self.verbose:
                    print("error, arg not found ")
                    print(arg)
                return

        # print(snow_model_args_copy.keys())

        # for k,v in self.snow_model_args.items():
        #     if v == 'x':
        #         print(k, " --> ", x[i])
        #         snow_model_args_copy.update({k:x[i]})
        #         i+=1
        # for k,v in self.hydro_model_args.items():
        #     if v == 'x':
        #         print(k, " --> ", x[i])
        #         hydro_model_args_copy.update({k:x[i]})
        #         i+=1



        # adjust the precip by the "pmult" value
        # get pmult first

        # run the snowmodel
        _snow = sm.snow17driver(**snow_model_args_copy)
        snow_out  = {"outflow":_snow[0],
                     "swe":_snow[1],
                     "rain_total":_snow[2],
                     "precip_total":_snow[3]}

        # adjust/fix the snow outflow
        nlayers = self.snow_model_args.get('nlayers')
        qin = (self.elevation_bins_areas.reshape(nlayers,1)  * snow_out['outflow']).sum(axis=0)

        # update the "qin"
        hydro_model_args_copy.update({"qin":qin})

        # run the hydro model
        hydro_out = dr.model_driver(**hydro_model_args_copy)
        hydro_out = {"qchan" : hydro_out[0],
                     "evec"  : hydro_out[1],
                     "sm"    : hydro_out[2],
                     "qin"    : qin}

        # return some stuff...
        return hydro_out





if __name__ == "__main__":

    start_date = "1988-10-01"
    end_date   = "1989-09-30"
    nlayers = 100
    model = model_wrapper()
    model.read_parameter_files()
    model.create_forcings(start_date, end_date, nlayers)


    basic_parameters = "basic_parameters.yml"
    # open the parameter dictionaries...
    with open(basic_parameters, "r") as p:
        dct = yaml.load(p, yaml.FullLoader)


    # set the initial conditions in the model
    model.snow_model_args.update(rvs=1, opg_method=2)
    model.hydro_model_args.update(dict(sm1i=.1, sm2i=.1, hydroop=0))


    # parse a message to print out
    options_msg = "PRIOR PARAMETERS \n"
    options_msg = options_msg + "\n".join(["%s %s"%(k,v) for i in ['snow', 'hydro'] for k,v in dct[i].items()])

    # #
    # test_conditions = {
    #                   "snow":
    #                   {"pm" :0.},
    #                    "hydro":
    #                    {
    #                    "sm1i": .1,
    #                    "sm2i":.1,
    #                    "sm1max":100.,
    #                    "sm2max":400.,
    #                    "c":1.,
    #                    "ku":10.,
    #                    "ks":20,
    #                    "lowercasen":1.,
    #                    "beta":1.}
    #                    }

    model.hydro_model_args.update(dct["hydro"])
    model.snow_model_args.update( dct["snow"])
    model.runmain([])







    # parse the output and pass into model
    #model.substitute_arguments(dct['snow'].keys(), dct['hydro'].keys())


    # join the model args into one dict
    #output = model.runmain([])


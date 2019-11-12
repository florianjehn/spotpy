# -*- coding: utf-8 -*-
"""
Created on Thu Oct 02 13:13:19 2019

@author: Florian Ulrich Jehn
"""

import spotpy
from spotpy.examples.hymod_python.hymod import hymod
import pandas as pd
import os
import sys

# add the whole package to the path
file_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.sep.join(file_dir.split(os.sep)[:-1]))

    
def go_through_all_catchments(dataframes):
    """
    Cycles through all catchments (but uses all years) seperately and saves the KGE to
    a new dataframes (and then a csv)
    """
    all_catchments = pd.DataFrame(index = list(dataframes.keys()), 
                                      columns = ["ns_calib", "ns_valid", "bias_calib", "bias_valid"])
    for catch in dataframes.keys():
        mymodel= spot_setup_hymod(dataframes[catch].E, dataframes[catch].P, dataframes[catch].Q)
        sampler = spotpy.algorithms.lhs(mymodel,dbname="bla", dbformat="ram")
        sampler.sample(100)
        all_catchments.loc[catch, "ns_calib"] = max(sampler.getdata()["like1"])
        all_catchments.loc[catch, "ns_valid"] = max(sampler.getdata()["like2"])
        all_catchments.loc[catch, "bias_calib"] = max(sampler.getdata()["like3"])
        all_catchments.loc[catch, "bias_valid"] = max(sampler.getdata()["like4"])        
        
    all_catchments.to_csv("kge_all_catchments.csv", sep=";")     


class spot_setup_hymod(object):
    cmax  = spotpy.parameter.Uniform(low=1.0 , high=500,  optguess=250)
    bexp  = spotpy.parameter.Uniform(low=0.1 , high=2.0,  optguess=1)
    alpha = spotpy.parameter.Uniform(low=0.1 , high=0.99, optguess=0.5)
    Ks    = spotpy.parameter.Uniform(low=0.001 , high=0.10, optguess=0.05)
    Kq    = spotpy.parameter.Uniform(low=0.1 , high=0.99, optguess=0.5)
        
    def __init__(self, PET, Precip, trueObs):
        self.PET = PET
        self.Precip = Precip
        self.trueObs = trueObs
        
    def simulation(self,x):
        return hymod(self.Precip, self.PET, x[0], x[1], x[2], x[3], x[4])
        
    def evaluation(self):
        return self.trueObs
    
    def objectivefunction(self,simulation,evaluation, params=None):
        sim = pd.DataFrame(simulation, index=pd.date_range(start="01.01.1991", 
                                                           periods = len(simulation), 
                                                           freq="d"))
        evalu = pd.DataFrame(evaluation, index=pd.date_range(start="01.01.1991", 
                                                           periods = len(evaluation), 
                                                           freq="d"))
        # Remove empty days
        evalu.dropna(inplace=True)
        # Use only the sim days where we have evaluation data
        sim = sim.loc[evalu.index]
        # Split in calibration/validation
        evalu_calib = evalu.loc["2000-01-01": "2009-12-31",:]
        evalu_valid = evalu.loc["2009-12-31":,:]
        sim_calib = sim.loc["2000-01-01": "2009-12-31",:]
        sim_valid = sim.loc["2009-12-31":,:]
        ns_calib = spotpy.objectivefunctions.nashsutcliffe(evalu_calib, sim_calib)
        ns_valid = spotpy.objectivefunctions.nashsutcliffe(evalu_valid, sim_valid)
        bias_calib = spotpy.objectivefunctions.bias(evalu_calib, sim_calib)
        bias_valid = spotpy.objectivefunctions.bias(evalu_valid, sim_valid)
        
        return [ns_calib, ns_valid, bias_calib, bias_valid]
    

if __name__ == "__main__":
   import preprocessing.cleaned_data.create_cleaned_data_table as ccdt
   dataframes = ccdt.get_table_dict(calc_water_year=False)
   dataframes = {k: dataframes[k] for k in list(dataframes.keys())[:2]}
   go_through_all_catchments(dataframes)

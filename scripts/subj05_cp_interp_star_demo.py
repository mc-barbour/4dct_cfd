#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 11:38:24 2023

@author: mbarb1
"""


import pyvista as pv
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt


from scipy.interpolate import CubicSpline
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots

from source.star_control_points import *


pio.renderers.default = "browser"
colors = px.colors.sequential.Jet


#%% Load the cp locations as calulcated from slicer. Prior to this the control points are defined using remeshing, manually tweaked, and loaded into slicer with the dynamic registration computed

work_dir = '/Users/mbarb1/OneDrive - UW/RobinSequence/Data/Area and Volume Data/4DCT_analysis/Subj05/dynamic_motion_testing/'
cmap =  plt.cm.get_cmap("jet", 10)
cp_files = sorted(glob.glob(work_dir + "transformed_controlPoints_18_24_2mmSpacing_fullAirway/*.fcsv"))

dfs_save = []
p = pv.Plotter()
polys = []
for count,file in enumerate(cp_files):
    
    print(file)
    df_mod = pd.read_csv(file, skiprows=3, 
                     names=['id','X','Y','Z','ow','ox','oy','oz','vis','sel','lock','label','desc','associatedNodeID'])

    df_star = df_mod[["X", "Y", "Z"]]
    dfs_save.append(df_star)

    poly = pv.PolyData(np.array([df_mod.X, df_mod.Y, df_mod.Z]).T)
    polys.append(poly)

    p.add_mesh(poly, color=cmap(count)[0:3])
p.show()



#%% overwrite the last timestep with the first

df_periodic = dfs_save.copy()
df_periodic[-1] = df_periodic[0]


#%% Stack all of the control points into arrays

n_images = len(df_periodic)
n_points = len(df_periodic[0])

dt = 0.1
time = np.linspace(0, (n_images-1)*dt, n_images)

x_all = np.zeros((n_images, n_points))
y_all = np.zeros((n_images, n_points))
z_all = np.zeros((n_images, n_points))


for count,df in enumerate(df_periodic):
    x = df.X
    y = df.Y
    z = df.Z
    
    x_all[count,:] = x
    y_all[count,:] = y
    z_all[count,:] = z


#%% Run the interpolation
       
dt_cfd = 1e-4 
period_length = 0.7
new_time = np.arange(0, period_length+dt_cfd, dt_cfd)    
x_new, y_new, z_new = interpolate_controlPoints_time(df_periodic, 0.1, new_time, show=True)



#%% create new df structure for repeating and saving - Total Displcement


df_periodic_star = periodic_star_table_totalDisp_fromArrays(x_new, y_new, z_new, dt=0.01, n_cycles=3)  
df_periodic_star = df_periodic_star*1e-3
df_periodic_star.to_csv(work_dir +"StarControlPoinstFull_2mmSpacing_Periodic_5Cycles_Total_interpolated_01s.csv", index=False)

#%% Create the datafrom for star - incremental displacement

df_periodic_star = periodic_star_table_IncDisp_fromArraysV2(x_new, y_new, z_new, dt=dt_cfd, n_cycles=1, start_time=1.4)  
df_periodic_star = df_periodic_star*1e-3
df_periodic_star.to_csv(work_dir +"StarControlPoinstFull_2mmSpacing_Periodic_5Cycles_Inc_interpolated_0001s_3rdCycle.csv", index=False)



#%% Get the inlet information
surf = pv.read('/Users/mbarb1/OneDrive - UW/RobinSequence/Data/CFD/MovingMesh/Subj05/Airway_18_FE_simready.inlet.stl')
surf.plot()
surf.center()










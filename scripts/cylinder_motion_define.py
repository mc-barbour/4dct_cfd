#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 11:15:26 2024

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


#%% Read in control point motion for the cylinder test motion file


df = pd.read_csv("/Users/mbarb1/OneDrive - UW/RobinSequence/Data/CFD/MovingMesh/cylinder_control_points_time.csv")

p0 = pv.PolyData(df.iloc[:,0:3].values)
p1 = pv.PolyData(df.iloc[:,6:9].values)
p2 = pv.PolyData(df.iloc[:,9:12].values)
p3 = pv.PolyData(df.iloc[:,12:15].values)

cmap =  plt.cm.get_cmap("jet", 4)

cp_polys = [p0,p1,p2,p3]

p = pv.Plotter()
df_save = []



for count,poly in enumerate(cp_polys):
    
    p.add_mesh(poly, color=cmap(count)[0:3], label=str(count), render_points_as_spheres=True)
    df_star = pd.DataFrame(data = poly.points, columns=["X","Y","Z"])
    df_save.append(df_star)
    
p.background='white' 
p.add_legend()
p.show()


#%% structure the data and interpolate
df_save.append(df_save[0])

new_time = np.arange(0, 0.4+0.01, 0.01)   
x_new, y_new, z_new = interpolate_controlPoints_time(df_save, 0.1, new_time, show=True) # first and last points are the same



#%% convert to incrimental displacement 

df_periodic_star = periodic_star_table_IncDisp_fromArraysV2(x_new, y_new, z_new, dt=0.01, n_cycles=1, start_time = 0.0)  

df_periodic_star.to_csv("/Users/mbarb1/OneDrive - UW/RobinSequence/Data/CFD/MovingMesh/cylinder_control_points_0.01_incDisp_1stCycle.csv", index=False)


#%% convert to incrimental displacement 

df_periodic_star = periodic_star_table_IncDisp_fromArraysV2(x_new, y_new, z_new, dt=0.01, n_cycles=1, start_time = 0.4)  

df_periodic_star.to_csv("/Users/mbarb1/OneDrive - UW/RobinSequence/Data/CFD/MovingMesh/cylinder_control_points_0.01_incDisp_2ndCycle.csv", index=False)

#%% convert to incrimental displacement 

df_periodic_star = periodic_star_table_IncDisp_fromArraysV2(x_new, y_new, z_new, dt=0.01, n_cycles=1, start_time = 0.8)  

df_periodic_star.to_csv("/Users/mbarb1/OneDrive - UW/RobinSequence/Data/CFD/MovingMesh/cylinder_control_points_0.01_incDisp_3rdCycle.csv", index=False)






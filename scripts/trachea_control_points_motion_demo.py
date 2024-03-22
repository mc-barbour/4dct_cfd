#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:07:42 2024

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

#%% Load in the original and remeshed surfaces

work_dir = '/Volumes/files/RadResearch/Projects/TracheoMalacia/DynamicCTMalaciaExample_1/'

surface_file_dir = work_dir + "CFD_surfaces/"

surf = pv.read(surface_file_dir + "Airway_0_CFD_cap.stl")
surfCoarse = pv.read(work_dir+ "Airway_0_controlPoionts_subsample500.ply")

surf.scale(1e3, inplace=True)
surfCoarse.scale(1e3, inplace=True)

p = pv.Plotter()
p.add_mesh(surf)
p.add_mesh(surfCoarse.points, render_points_as_spheres=True, color='red')
p.show()

#%% Manually select points to exclude

coarse_pc = pv.PolyData(surfCoarse.points)
remove_points, remove_ids = pointPickerRemove(coarse_pc, surf)


cp_df = pd.DataFrame(data = coarse_pc.points, columns=["l", "p", "s"])
cp_df.drop(remove_ids, inplace=True)


new_cp = pv.PolyData(np.array([cp_df.l.values, cp_df.p.values, cp_df.s.values]).T)

p = pv.Plotter()
p.add_mesh(surf)
p.add_mesh(new_cp, render_points_as_spheres=True, point_size=15, color='blue')
p.show()

#%% Add control points

new_points, pt_ids = pointPickerAdd(new_cp, surf)

cp_df = pd.DataFrame(data = new_cp.points, columns=["l", "p", "s"])

df_add = pd.DataFrame(data = new_points, columns=["l", "p", "s"])

new_cp_df = cp_df.append(df_add, ignore_index=True)

new_cp = pv.PolyData(np.array([new_cp_df.l.values, new_cp_df.p.values, new_cp_df.s.values]).T)



p = pv.Plotter()
p.add_mesh(surf)
p.add_mesh(new_cp, render_points_as_spheres=True, point_size=15, color='blue')
p.show()



#%% save trimmed control points as cp file to load slicer

df = pd.DataFrame(data = new_cp.points, columns=["l", "p", "s"])
df['selected'] = np.ones(len(df))
df['visible'] = np.ones(len(df))
df['locked'] = np.ones(len(df))

df.index.name='label'
save_name = work_dir + 'control_points/Airway_0_controlPoionts_subsample500.fcsv'
print(save_name)
df.to_csv(save_name)

#%% Loop and save the cp locations !!! This gets dropped into the slicer python interpretor window

save_dir = '/Volumes/files/RadResearch/Projects/TracheoMalacia/DynamicCTMalaciaExample_1/control_points_start0_subsample1000/'
prefix = "Airway_controlPoints_0_subsample1000_"
transformSequenceID = 'vtkMRMLSequenceNode3' # name of sequence node with registration

controlPointNode = getNode("Airway_0_controlPoionts_subsample1000")

shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
itemIDToClone = shNode.GetItemByDataNode(controlPointNode)


transformSeq = slicer.mrmlScene.GetNodeByID(transformSequenceID)


for node_idx in range(17):
    # get the transform
    
    transformNode = transformSeq.GetNthDataNode(node_idx)
    transform = transformNode.GetTransformToParent()
    
    #copy the original node
    clonedItemID = slicer.modules.subjecthierarchy.logic().CloneSubjectHierarchyItem(shNode, itemIDToClone)
    clonedControlPointNode = shNode.GetItemDataNode(clonedItemID)
    
    #apply transform
    clonedControlPointNode.ApplyTransform(transform)
    save_id = node_idx
    # Write to file
    outputFileName = save_dir + f"{prefix}_{save_id:03}.fcsv"
    print(outputFileName)
    slicer.util.saveNode(clonedControlPointNode, outputFileName)





#%% Load the cp locations as calulcated from slicer. Prior to this the control points are defined using remeshing, manually tweaked, and loaded into slicer with the dynamic registration computed

cmap =  plt.cm.get_cmap("jet", 10)
cp_files = sorted(glob.glob(work_dir + 'control_points_start0_subsample500/*.fcsv'))

dfs_save = []
p = pv.Plotter()
polys = []
for count,file in enumerate(cp_files):
    
    print(file)
    df_mod = pd.read_csv(file, skiprows=3, usecols=[0,1,2,3], names=['id','X','Y','Z'])

    df_star = df_mod[["X", "Y", "Z"]]
    dfs_save.append(df_star)

    poly = pv.PolyData(np.array([df_mod.X, df_mod.Y, df_mod.Z]).T)
    polys.append(poly)

    p.add_mesh(poly, color=cmap(count)[0:3])
p.show()



#%% create plot of each set of control points

dfs_skip = []
counts = [0,2,6,8,10,12,14]
p = pv.Plotter()
for count in counts:

    p.add_mesh(polys[count], color=cmap(count)[0:3], label=str(count))
    
    df_mod = pd.read_csv(cp_files[count], skiprows=3, usecols=[0,1,2,3], names=['id','X','Y','Z'])

    df_star = df_mod[["X", "Y", "Z"]]
    dfs_skip.append(df_star)

p.add_legend()
p.show()


df_periodic = dfs_skip.copy()
df_periodic.append(dfs_skip[0])


#%% plot the first and last files
p = pv.Plotter()
p.add_mesh(polys[0], color='red', label = "0")
p.add_mesh(polys[15], color='blue', label = "15")
p.add_legend()
p.show()


#%% overwrite the last timestep with the first


df_periodic = dfs_save.copy()
df_periodic[-1] = df_periodic[0]


#%% Stack all of the control points into arrays

n_images = len(df_periodic)
n_points = len(df_periodic[0])

dt = 0.2
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

#%% save just the control points - no interpolation
df_periodic_star = periodic_star_table_IncDisp_fromArraysV2(x_all, y_all, z_all, dt=0.2, n_cycles=3, start_time=0.0)  
df_periodic_star = df_periodic_star*1e-3
df_periodic_star.to_csv(work_dir +"StarControlPoinstFull_0start_subsample500_everyother_Periodic_3Cycles_Inc_noInterp.csv", index=False)


#%% Run the interpolation

dt_image= 0.2
dt_cfd = 0.01
period_length = 1.6
new_time = np.arange(0, period_length+dt_cfd, dt_cfd)    
x_new, y_new, z_new = interpolate_controlPoints_time(df_periodic, dt_image, new_time, show=True)
de


#%% create new df structure for repeating and saving - Total Displcement

df_periodic_star = periodic_star_table_totalDisp_fromArrays(x_new, y_new, z_new, dt=0.01, n_cycles=3)  
df_periodic_star = df_periodic_star*1e-3
df_periodic_star.to_csv(work_dir +"StarControlPoinstFull_0start_subsample1000_Periodic_3Cycles_Total_interpolated_001s.csv", index=False)

#%% Create the datafrom for star - incremental displacement

df_periodic_star = periodic_star_table_IncDisp_fromArraysV2(x_new, y_new, z_new, dt=dt_cfd, n_cycles=3, start_time=0.0)  
df_periodic_star = df_periodic_star*1e-3
df_periodic_star.to_csv(work_dir +"StarControlPoinstFull_0start_subsample500_everyother_Periodic_3Cycles_Inc_interpolated_01.csv", index=False)



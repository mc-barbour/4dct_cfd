#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:49:05 2025

@author: mbarb1
"""

from functools import partial
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pycpd import DeformableRegistration
import numpy as np
import pandas as pd
import pyvista as pv
import time
import glob
from natsort import natsorted
import pymeshlab

"""
This requires pycpd. I had a hard time getting pycpd with my larger 4dct python environment. I simply downloaded the source files and exectute the pycpd functions 
(including this script) from that source directory.

"""

def compute_displacement_vec(cp1, cp0):
    """
    computed the displacement between two sets of points
    input: dataframes
    outpt: vector (numpy), and magnitude
    """

    vec = cp1 - cp0
    mag = np.linalg.norm(vec, axis = 1)
    return vec, mag


def remesh_surfaces(target_surface_files, initial_surface_file, wrk_dir, remesh_dir, targetfacenum=15000):
    """
    Remesh surfaces with a specific surface density.
    
    target_surface_files: list of target filenames (list/array)
    inital_surface_file: single file name of initial/starting surface (string)
    wrk_dir: string of working directory
    remesh_dir: string of new directory to save remeshed surface files
    targetfacenum: target number of faces after decimation (int)
    
    """

    for file in target_surface_files:
        ms = pymeshlab.MeshSet()
        ms.load_new_mesh(file)
        base_name = file.split("/")[-1]
        ms.meshing_decimation_quadric_edge_collapse(targetfacenum=targetfacenum, preservetopology=True)
        ms.save_current_mesh(remesh_dir + "remeshed_"+base_name)
    
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(initial_surface_file)
    base_name = initial_surface_file.split("/")[-1]
    ms.meshing_decimation_quadric_edge_collapse(targetfacenum=targetfacenum, preservetopology=True)
    ms.save_current_mesh(wrk_dir + "remeshed_" + base_name)

def visualize_surfaces(target_surface_files, initial_surface_file=None):
    """
    Visualizes one or more 3D surface meshes using PyVista in a grid layout.

    Parameters:
    target_surface_files (list of str): A list of file paths to 3D surface files 
                                        (e.g., .ply, .vtp, .stl) to be visualized.
    initial_surface_file (str, optional): An optional initial surface file to be 
                                          displayed first in the visualization grid.

    Output:Opens an interactive PyVista window displaying all the surfaces.
    """
    
    
    if initial_surface_file:
        target_surface_files.insert(0,initial_surface_file)
        
    ncols = int(np.ceil(np.sqrt(len(target_surface_files))))
    nrows = int(np.ceil(len(target_surface_files) / ncols))

    p = pv.Plotter(shape=(nrows,ncols))
    
    for count, file in enumerate(target_surface_files):
        surf = pv.read(file)

        row = int(count / nrows)
        col = count - row * nrows
        print("Number of Points: ", surf.n_points)

        p.subplot(row, col)
        p.add_mesh(surf)
        p.add_text(file.split("/")[-1])                        

    p.link_views()
    p.show()
    
    
def define_clipping_planes(surf, z0, z1):
    """
    Visualize Clipping Plane locations

    Parameters
    ----------
    surf : polydata surface
        
    z0 : z-position of plane
        
    z1 : z-position of second plane

    """
    center= surf.center
    center[2] = z0
    plane0 = pv.Plane(center = center, i_size=20, j_size=20)
    
    center[2] = z1
    plane1 = pv.Plane(center = center, i_size=20, j_size=20)
    
    p = pv.Plotter()
    p.add_mesh(surf)
    p.add_mesh(plane0, color='green', label='Z0' )
    p.add_mesh(plane1, color='red', label='Z1')
    p.add_legend()
    p.show()
    
    
def clip_surface(surf, zTop=None, zBottom=None, show=True):
    
    if zTop:
        clipped = surf.clip(normal="z", origin=(0,0,zTop),invert=True)
    
    else:
        clipped=surf
    
    if zBottom:
        clipped_out = clipped.clip(normal="z", origin=(0,0,zBottom),invert=False)
    
    else:
        clipped_out = clipped
        
    if show:
        p = pv.Plotter()
        p.add_mesh(clipped_out)
        p.show()
    
    return clipped_out


def deform_points_staticRef(reference_surf, reference_cps, target_surface_files, deformed_surfaces_dir, deformed_points_dir, zTop=None, zBottom=None):
    """
    Performs deformable registration of multiple target surface meshes to a static reference surface,
    and saves the deformed results.

    Parameters:
        reference_surf (pyvista.PolyData): The fixed reference surface .
        reference_cps (pyvista.PolyData): Control points on the reference surface.
        target_surface_files (list of str): List of file paths to target surfaces taht are registered from the refence surface.
        deformed_surfaces_dir (str): Directory path to save the deformed (registered) target surfaces.
        deformed_points_dir (str): Directory path to save the deformed control points.
        zTop (float, optional): Upper Z-bound for clipping the target surfaces. If None, no upper clip is applied.
        zBottom (float, optional): Lower Z-bound for clipping the target surfaces. If None, no lower clip is applied.


    Output:
        - Saves deformed surface meshes and control points as `.ply` files in the designated output directories.
        - Prints timing and basic information for each registration process.

    """
        
    for count, surf_file_name in enumerate(target_surface_files):
        
        target_surf = pv.read(surf_file_name)
        target_num = surf_file_name.split(".")[0].split("_")[-1]
        base_name = surf_file_name.split("/")[-1].split(".")[0]
        print(surf_file_name, target_num)
            
        
        clipped_target_surf = clip_surface(target_surf, zTop=zTop, zBottom=zBottom, show=False)
        
        # compute the registration
        print(len(clipped_target_surf.points), len(reference_surf.points))
    
        
        start_time = time.time()
        reg = DeformableRegistration(**{'X': np.array(clipped_target_surf.points), 'Y': np.array(reference_surf.points)}) # X:Target, Y:Source
        reg.register()
        
        print("--- %s seconds ---" % (time.time() - start_time))
        
        # perform deformation
        moved_surface_points = reg.transform_point_cloud(Y=reference_surf.points)
        moved_control_points = reg.transform_point_cloud(Y=reference_cps.points)
        
        #save and update
        save_control_points = pv.PolyData(moved_control_points)
        save_control_points.save(deformed_points_dir + "deformed_" + base_name + ".ply")
        save_surface_points = pv.PolyData(moved_surface_points)
        save_surface_points.save(deformed_surfaces_dir + "deformed_" + base_name + ".ply")
        
def visualize_deformed_points_and_vectors(deformed_cp_files, ref_cp_file, target_surface_files, ref_surface_file):
    """
    visualize morphed control points wth displacment vectors
    """
    
    deformed_points_list = []
    
    surf0 = pv.read(ref_cp_file)
    surf1 = pv.read(deformed_cp_files[0])
    disp, mag = compute_displacement_vec(surf1.points, surf0.points) #target, #source
    surf0['disp'] = disp
    surf0.set_active_vectors('disp')
    deformed_points_list.append(surf0)
    
    for count, file in enumerate(deformed_cp_files[0:-1]):
        print(file)
        surf0 = pv.read(file)
        surf1 = pv.read(deformed_point_files[count+1])
        disp, mag = compute_displacement_vec(surf1.points, surf0.points) #target, #source
        surf0['disp'] = disp
        surf0.set_active_vectors('disp')
        deformed_points_list.append(surf0)

     
    target_surface_files.insert(0, ref_surface_file)  
    ncols = int(np.ceil(np.sqrt(len(target_surface_files))))
    nrows = int(np.ceil(len(target_surface_files) / ncols))
    
    p = pv.Plotter(shape=(nrows,ncols))
    p = pv.Plotter(shape=(nrows,ncols))

    for count, file in enumerate(target_surface_files[0:-1]):
        surf = pv.read(file)
        deformed_points = deformed_points_list[count]
        row = int(count / nrows)
        col = count - row * nrows    
        
        
        p.subplot(row, col)
        p.add_mesh(surf)
        p.add_mesh(deformed_points, render_points_as_spheres=True, color='red')     
        p.add_mesh(deformed_points.arrows)                 
        p.add_text(file.split("/")[-1].split(".")[0])
        
    p.link_views()
    p.show()
        
def gather_deformed_cps(deformed_points_files, ref_point_file):

    ref_cp = pv.read(ref_point_file)
    
    points_list = []
    points_list.append(starting_points.points)
    for file in deformed_point_files:
        cps = pv.read(file)
        points_list.append(cps.points)
    points_list.append(starting_points.points)
    
    return points_list

def interpolate_cp_motion(points_list, dt_image=0.1, dt_cfd=1e-3 ):
    """
    Interpolate the motion between points_list at specified dt_cfd timespacing
    """

    morph_points_df = []
    
    for points in points_list:
    
        df = pd.DataFrame(points, columns=["X", "Y", "Z"])
        morph_points_df.append(df)
    
    period_length = (len(morph_points_df) - 1) * dt_image
    
    new_time = np.arange(0, period_length+dt_cfd, dt_cfd)    
    x_new, y_new, z_new = interpolate_controlPoints_time(morph_points_df, dt_image, new_time, show=True)
    
    return x_new, y_new, z_new

def save_star_motionFiles(x_new, y_new, z_new, period_length, save_dir, save_prefix, dt_cfd=1e-3):
    """
    Save control point locations for Star CCM morphing. One control point displaemnt per saved file
    """
    
    
    start_time = 0
    end_time = period_length * 2


    iter = np.arange(start_time, end_time, dt_cfd)

    for i in iter:
        df_periodic_star = reference_periodic_star_table_IncDisp_fromArrays_single_dt(x_new, y_new, z_new, dt=dt_cfd, start_time=i, cycle_len=period_length)
        df_periodic_star = df_periodic_star*1e-3
        df_periodic_star.to_csv(save_dir + save_prefix + "_" + f"{(i+dt_cfd):.5f}.csv", index=False)


    
#%% Let's first remesh - decimate the surfaces

wrk_dir = "/Users/mbarb1/OneDrive - UW/TracheaMalacia/CFD/Patients/TM9_debug/"
new_dir = wrk_dir + "remeshed_target_surfaces/"
target_surface_files = natsorted(glob.glob(wrk_dir + 'target_surface_v2/*.stl'))
init_surface_file = wrk_dir + "TM9_13.stl"

remesh_surfaces(target_surface_files, init_surface_file, wrk_dir, new_dir)

#%%
remeshed_init_surface_file = wrk_dir + "remeshed_TM9_13.stl"
remeshed_target_surface_files = natsorted(glob.glob(wrk_dir + 'remeshed_target_surfaces/*.stl'))
visualize_surfaces(remeshed_target_surface_files, remeshed_init_surface_file)

#%% define clipping planes
surf = pv.read(remeshed_init_surface_file)
print(surf.bounds) # looks at the last two elements of bounda to find max and min locations of z

define_clipping_planes(surf, z0=-85, z1=-110) # adjust these values until you get good planes for clipping.

#%% Load and clip inital surfaces
init_surface = pv.read(wrk_dir + "remeshed_TM9_13.stl")
init_surface_clipped = clip_surface(init_surface,zBottom=-112)


control_points = pv.read(wrk_dir + "TM9_13_downsample_coarse.ply")
control_points_clipped = clip_surface(control_points, zTop=-87, zBottom=-110)
control_points_clipped.save(wrk_dir + "inital_cps.ply")

#%% Deform cps


deformed_points_dir = wrk_dir +"deformed_points_staticStart/"
deformed_surfaces_dir = wrk_dir + "deformed_surfaces_staticStart/"
target_surface_files = natsorted(glob.glob(wrk_dir + 'remeshed_target_surfaces/*.stl'))

deform_points_staticRef(init_surface_clipped, control_points_clipped, target_surface_files, deformed_surfaces_dir, deformed_points_dir, zBottom=-112)


#%% visualize deformed points 
target_surface_files = natsorted(glob.glob(wrk_dir + 'remeshed_target_surfaces/*.stl'))
deformed_cp_files = natsorted(glob.glob(deformed_points_dir + "*.ply"))
init_cp_file = wrk_dir + "inital_cps.ply"


visualize_deformed_points_and_vectors(deformed_cp_files, init_cp_file, target_surface_files, init_surface_file)



#%% Gather deformed_points and interpolate cp motion - need to change to 4dct directory for remainder of script
from source.star_control_points import *

dt_image = 0.1
dt_cfd = 0.01

deformed_points_list = gather_deformed_cps(deformed_cp_files, init_cp_file)
x_new, y_new, z_new = interpolate_cp_motion(deformed_points_list, dt_cfd=dt_cfd, dt_image=dt_image)



#%% save motion vectors in Starccm format
save_dir = '/Users/mbarb1/OneDrive - UW/TracheaMalacia/CFD/Patients/TM9_debug/morpher_files_dt01/'
save_prefix = "TM9_motion_singleDT_01"
period_length = (len(deformed_points_list)-1)*dt_image
save_star_motionFiles(x_new, y_new, z_new, period_length, save_dir, save_prefix, dt_cfd=dt_cfd)



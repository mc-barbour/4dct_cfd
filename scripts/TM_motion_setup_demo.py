#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:49:05 2025

@author: mbarb1
"""

from functools import partial
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# from pycpd import DeformableRegistration
import numpy as np
import pandas as pd
import pyvista as pv
import time
import glob
from natsort import natsorted

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


#%% Let's start by viewing the surfaces. Airway_8.stl will be the inital surface

wrk_dir = "/Users/mbarb1/dev/pycpd/airway_test/"


target_surface_files = natsorted(glob.glob(wrk_dir + 'target_surfaces_v5/*.stl'))
nrows = 4
ncols = 4

p = pv.Plotter(shape=(nrows,ncols))



for count, file in enumerate(target_surface_files):
    surf = pv.read(file)
    
    row = int(count / nrows)
    col = count - row * nrows    

    p.subplot(row, col)
    p.add_mesh(surf)
    p.add_text(file.split("/")[-1])                        

p.link_views()
p.show()


#%% Might be helpful to clip the surfaces


moving_points = pv.read('airway_test/Airway_8_downsample.ply')
target1 = pv.read(target_surface_files[6])


clipped_pc = moving_points.clip(normal="z", origin=(0,0,-10),invert=True)
clipped_moving_points = clipped_pc.clip(normal="z", origin=(0,0,-50), invert=False)


clipped_surface = target1.clip(normal="z", origin=(0,0,-10),invert=True)
clipped_target_surf = clipped_surface.clip(normal="z", origin=(0,0,-50), invert=False)

p = pv.Plotter()
p.add_mesh(clipped_moving_points, render_points_as_spheres=True, color='red')
p.add_mesh(clipped_target_surf, color='grey')
p.add_axes()
p.show()



#%% Perform the surface registration


moving_points = pv.read('airway_test/Airway_8_v5_recon_downsample.ply')
moving_points_coarse = pv.read('airway_test/Airway_8_v5_recon_downsample_coarse.ply') # meshlab resampling 

output_dir = wrk_dir + "TM3_surfaces_v5_deformed_points"

target_surface_files = natsorted(glob.glob('airway_test/target_surfaces_v5/*.stl'))
target_surface_order = [6,7,8,9,10,11,12,13,0,1,2,3,4,5]
# target_surface_order = [7, 9, 11, 2, 4]

starting_surface = pv.read('airway_test/Airway_8_v5_recon.stl')
starting_surface_tmp = starting_surface.clip(normal="z", origin=(0,0,-10),invert=True)
starting_surface_clipped = starting_surface_tmp.clip(normal="z", origin=(0,0,-50), invert=False)


clipped_pc = moving_points.clip(normal="z", origin=(0,0,-10),invert=True)
clipped_moving_points = clipped_pc.clip(normal="z", origin=(0,0,-50), invert=False)

clipped_pc_coarse = moving_points_coarse.clip(normal="z", origin=(0,0,-10),invert=True)
clipped_moving_points_coarse = clipped_pc_coarse.clip(normal="z", origin=(0,0,-50), invert=False)


moving_points = clipped_moving_points.points
moving_points_coarse = clipped_moving_points_coarse.points
moving_points_fine = starting_surface_clipped.points

moving_points_list = []
for count, surf_num in enumerate(target_surface_order[0:3]):
    target_surf = pv.read(target_surface_files[surf_num])
    target_num = target_surface_files[surf_num].split(".")[0].split("_")[-1]
    print(target_surface_files[surf_num], target_num)
        
    clipped_surface = target_surf.clip(normal="z", origin=(0,0,-10),invert=True)
    clipped_target_surf = clipped_surface.clip(normal="z", origin=(0,0,-50), invert=False)
        
    # compute the registration    
    start_time = time.time()
    reg = DeformableRegistration(**{'X': np.array(clipped_target_surf.points), 'Y': np.array(moving_points_fine)}) # X:Target, Y:Source
    reg.register()
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
    # perform deformation
    moved_points = reg.transform_point_cloud(Y=moving_points)
    moved_points_coarse = reg.transform_point_cloud(Y=moving_points_coarse)
    moved_points_fine = reg.transform_point_cloud(Y=moving_points_fine)

            
    disp, mag = compute_displacement_vec(moved_points_fine, moving_points_fine)
    new_points = pv.PolyData(moving_points_fine)
    new_points['disp'] = disp
    new_points.set_active_vectors('disp')
    moving_points_list.append(new_points)

    # update moving points variable
    moving_points = moved_points
    moving_points_coarse = moved_points_coarse
    moving_points_fine = moved_points_fine
    
    # new_poly_fine = pv.PolyData(moved_points_fine)
    # new_poly_fine.save("airway_test/deformed_points_v5_fine_2skip_fullCPD/Airway_8_fine_reg_" + target_num + ".ply")
    
    new_poly_coarse = pv.PolyData(moved_points_coarse)
    new_poly_coarse.save(output_dir + "Airway_8_coarse_reg_" + target_num + ".ply")
    
    # new_poly = pv.PolyData(moved_points)
    # new_poly.save("airway_test/deformed_points_v5_2skip/Airway_8_reg_" + target_num + ".ply")
    
    
    # p = pv.Plotter()
    # p.add_mesh(target1, opacity = 0.5)
    # p.add_mesh(moved_points, color='blue', render_points_as_spheres=True)
    # p.add_mesh(moving_points, color='red', render_points_as_spheres=True)
    # # p.add_mesh(new_points.arrows)
    # p.show()




#%%Airway_8 is the starting point

starting_surface = pv.read('airway_test/Airway_8_v5_recon.stl')
starting_surface = pv.read('airway_test/Airway_8_v5_recon_downsample_coarse.ply')
starting_surface_tmp = starting_surface.clip(normal="z", origin=(0,0,-10),invert=True)
starting_surface_clipped = starting_surface_tmp.clip(normal="z", origin=(0,0,-50), invert=False)

starting_surface_clipped.plot()


#%% Now read in the deformed points 

target_surface_order = [6,7,8,9,10,11,12,13,0,1,2,3,4,5] 
# target_surface_order = [7, 9, 11, 2, 6]
# target_surface_order = [2,3,4,0,1] # fine_skip2
# deformed_point_files = natsorted(glob.glob("airway_test/deformed_points_v5_fine_2skip_fullCPD/*.ply"))
deformed_point_files = natsorted(glob.glob("airway_test/deformed_points_v5_coarse/*.ply"))
source_points = starting_surface_clipped.points


points_list = []
points_list.append(source_points)

for count, surf_num in enumerate(target_surface_order):
    deformed_points = pv.read(deformed_point_files[surf_num])
    print(deformed_point_files[surf_num], len(deformed_points.points))
    points_list.append(deformed_points.points)
    
points_list.append(source_points)

#%% view the deformed points

nrows=4
ncols=4

p = pv.Plotter(shape=(nrows,ncols))


for count, surf_num in enumerate(target_surface_order):
    surf_name = target_surface_files[surf_num].split("/")[-1]
    target_surf = pv.read(target_surface_files[surf_num])
    deformed_points = pv.read(deformed_point_files[surf_num])
    row = int(count / nrows)
    col = count - row * nrows    

    p.subplot(row, col)
    p.add_text(surf_name)
    p.add_mesh(target_surf, opacity = 0.5, color='cyan')
    p.add_mesh(deformed_points, render_points_as_spheres=True)                      

p.link_views()
p.show()



#%% convert from polydata to dataframe
morph_points_df = []

for points in points_list:

    df = pd.DataFrame(points, columns=["X", "Y", "Z"])
    morph_points_df.append(df)
    

#%% convert from dataframe to x,y,z arrays

n_images = len(morph_points_df)
n_points = len(df)

x_all = np.zeros((n_images, n_points))
y_all = np.zeros((n_images, n_points))
z_all = np.zeros((n_images, n_points))


for count,df in enumerate(morph_points_df):
    x = df.X
    y = df.Y
    z = df.Z
    
    x_all[count,:] = x
    y_all[count,:] = y
    z_all[count,:] = z


#%% interpolate motion between surfaces at desired cfd timestep - need to change to 4dct directory for remainder of script
from source.star_control_points import *



dt_image= 0.1
dt_cfd = 0.01
period_length = 1.5 

new_time = np.arange(0, period_length+dt_cfd, dt_cfd)    
x_new, y_new, z_new = interpolate_controlPoints_time(morph_points_df, dt_image, new_time, show=True)



#%% now convert the interpolated control points to excel files for starccm to read


end_time = 2.1
start_time = 0.0
period_length = 1.5
dt_files = 0.1 # frequency to split csv motion files
div_per_cycle = period_length / dt_files
iter = np.arange(start_time, end_time, dt_files)

for i in iter:
    df_periodic_star = star_motion_table_split_incDisp(x_new, y_new, z_new, dt=dt_cfd, start_time=i, div_per_cycle=div_per_cycle, cycle_len=period_length)
    df_periodic_star = df_periodic_star*1e-3  # convert from mm to m
    df_periodic_star.to_csv(wrk_dir + "starccm_motion_files/cp_files_01/" + f"test_new_cp_function_{i:.3f}.csv", index=False) 



















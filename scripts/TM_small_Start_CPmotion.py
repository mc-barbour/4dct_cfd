#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 10:02:13 2024

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


#%%

def compute_displacement_vec(cp1, cp0):
    """
    computed the displacement between two sets of points
    input: dataframes
    outpt: vector (numpy), and magnitude
    """

    vec = cp1 - cp0
    mag = np.linalg.norm(vec, axis = 1)
    return vec, mag


#%% Let's look at all of the target surfaces


target_surface_files = natsorted(glob.glob('airway_test/target_surfaces_v5/*.stl'))
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

#%% 

surf9 = pv.read(target_surface_files[6])
surf10 = pv.read(target_surface_files[7])
surf11 = pv.read(target_surface_files[8])
p = pv.Plotter()
p.add_mesh(surf9, opacity = 0.5, color='cyan')
p.add_mesh(surf10, opacity = 0.5)

p.add_mesh(surf11, opacity = 0.5)


#%% try smoothing it looks like the taubin smooth function doesnt exist in this version of pyvista

surf8 = pv.read(target_surface_files[5])
surf8_smooth = surf8.smooth_taubin(band_pass=0.1)

p = pv.Plotter(shape=(1,2))
p.subplot(0,0)
p.add_mesh(surf8)

p.subplot(0,1)
p.add_mesh(surf8_smooth)
p.link_views()
p.show()

#%% let's clip the source and target surfaces

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



start_time = time.time()
reg = DeformableRegistration(**{'X': np.array(clipped_target_surf.points), 'Y': np.array(clipped_moving_points.points)}) # X:Target, Y:Source
reg.register()

print("--- %s seconds ---" % (time.time() - start_time))


moved_points = reg.transform_point_cloud(Y=clipped_moving_points.points)
    
disp, mag = compute_displacement_vec(moved_points, clipped_moving_points.points)
new_points = pv.PolyData(clipped_moving_points)
new_points['disp'] = disp
new_points.set_active_vectors('disp')


p = pv.Plotter()
p.add_mesh(target1, opacity = 0.5)
p.add_mesh(moved_points, color='blue', render_points_as_spheres=True)
p.add_mesh(clipped_moving_points, color='red', render_points_as_spheres=True)
p.add_mesh(new_points.arrows)
p.show()
#%% let's try using the smallest diameter mesh as the startig point


moving_points = pv.read('airway_test/Airway_8_v5_recon_downsample.ply')
moving_points_coarse = pv.read('airway_test/Airway_8_v5_recon_downsample_coarse.ply')
target_surface_files = natsorted(glob.glob('airway_test/target_surfaces_v5/*.stl'))
target_surface_order = [6,7,8,9,10,11,12,13,0,1,2,3,4,5]
target_surface_order = [7, 9, 11, 2, 4]

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
for count, surf_num in enumerate(target_surface_order):
    target_surf = pv.read(target_surface_files[surf_num])
    target_num = target_surface_files[surf_num].split(".")[0].split("_")[-1]
    print(target_surface_files[surf_num], target_num)
        
    clipped_surface = target_surf.clip(normal="z", origin=(0,0,-10),invert=True)
    clipped_target_surf = clipped_surface.clip(normal="z", origin=(0,0,-50), invert=False)
        
        
    start_time = time.time()
    reg = DeformableRegistration(**{'X': np.array(clipped_target_surf.points), 'Y': np.array(moving_points_fine)}) # X:Target, Y:Source
    reg.register()
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
    
    moved_points = reg.transform_point_cloud(Y=moving_points)
    moved_points_coarse = reg.transform_point_cloud(Y=moving_points_coarse)
    moved_points_fine = reg.transform_point_cloud(Y=moving_points_fine)

            
    disp, mag = compute_displacement_vec(moved_points_fine, moving_points_fine)
    new_points = pv.PolyData(moving_points_fine)
    new_points['disp'] = disp
    new_points.set_active_vectors('disp')
    moving_points_list.append(new_points)

    moving_points = moved_points
    moving_points_coarse = moved_points_coarse
    moving_points_fine = moved_points_fine
    
    new_poly_fine = pv.PolyData(moved_points_fine)
    new_poly_fine.save("airway_test/deformed_points_v5_fine_2skip_fullCPD/Airway_8_fine_reg_" + target_num + ".ply")
    
    # new_poly_coarse = pv.PolyData(moved_points_coarse)
    # new_poly_coarse.save("airway_test/deformed_points_v5_coarse_2skip/Airway_8_coarse_reg_" + target_num + ".ply")
    
    # new_poly = pv.PolyData(moved_points)
    # new_poly.save("airway_test/deformed_points_v5_2skip/Airway_8_reg_" + target_num + ".ply")
    
    
    # p = pv.Plotter()
    # p.add_mesh(target1, opacity = 0.5)
    # p.add_mesh(moved_points, color='blue', render_points_as_spheres=True)
    # p.add_mesh(moving_points, color='red', render_points_as_spheres=True)
    # # p.add_mesh(new_points.arrows)
    # p.show()


#%%
nrows = 4
ncols = 4

p = pv.Plotter(shape=(nrows,ncols))



for count, surf in enumerate(moving_points_list):

    row = int(count / nrows)
    col = count - row * nrows    

    p.subplot(row, col)
    p.add_mesh(surf, opacity = 0.5)
    p.add_mesh(surf.arrows)                      

p.link_views()
p.show()

#%% let's plot the motion fields

nrows = 4
ncols = 4

p = pv.Plotter(shape=(nrows,ncols))



for count, surf_num in enumerate(target_surface_order):
    target_surf = pv.read(target_surface_files[surf_num])
    row = int(count / nrows)
    col = count - row * nrows    

    p.subplot(row, col)
    p.add_mesh(target_surf, opacity = 0.5)
    p.add_mesh(moving_points_list[count].arrows)                      

p.link_views()
p.show()

#%%
nrows = 1
ncols = 4

p = pv.Plotter(shape=(nrows,ncols))


previous_surf = pv.read('airway_test/Airway_8_taubin_smooth.stl')
for count, surf_num in enumerate(target_surface_order[0:4]):
    
    target_surf = pv.read(target_surface_files[surf_num])
    row = int(count / nrows)
    col = count - row * nrows    

    p.subplot(0, count)
    p.add_mesh(target_surf, opacity = 0.25)
    p.add_mesh(moving_points_list[count].arrows)
    p.add_mesh(previous_surf, opacity = 0.25, color='cyan')

    previous_surf = target_surf                    

p.link_views()
p.show()




#%% now let's create the array of points for interpolation

target_surface_order = [6,7,8,9,10,11,12,13,0,1,2,3,4,5]
target_surface_order = [7, 9, 11, 2, 6]
target_surface_order = [2,3,4,0,1]
deformed_point_files = natsorted(glob.glob("airway_test/deformed_points_v5_fine_2skip_fullCPD/*.ply"))
source_points = starting_surface_clipped.points


points_list = []
points_list.append(source_points)

for count, surf_num in enumerate(target_surface_order):
    deformed_points = pv.read(deformed_point_files[surf_num])
    print(deformed_point_files[surf_num], len(deformed_points.points))
    points_list.append(deformed_points.points)
    
points_list.append(source_points)




#%%
nrows = 3
ncols = 3

p = pv.Plotter(shape=(nrows,ncols))


for count, surf in enumerate(points_list):
    
    row = int(count / nrows)
    col = count - row * nrows    

    p.subplot(row, col)
    p.add_mesh(surf, render_points_as_spheres=True)
p.link_views()
p.show()

p = pv.Plotter()


for count, surf in enumerate(points_list):
    
    row = int(count / nrows)
    col = count - row * nrows    


    p.add_mesh(surf, render_points_as_spheres=True)

p.show()



#%% conver to dataframe
morph_points_df = []

for points in points_list:

    df = pd.DataFrame(points, columns=["X", "Y", "Z"])
    morph_points_df.append(df)
#%%  
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


#%%create a file with no interpolation
work_dir = '/Users/mbarb1/OneDrive - UW/TracheaMalacia/CFD/Motion_test/Small_start/'
df_periodic_star = periodic_star_table_IncDisp_fromArraysV2(x_all, y_all, z_all, dt=0.1, n_cycles=1, start_time=0)  
df_periodic_star = df_periodic_star
df_periodic_star.to_csv(work_dir +"StarControlPoints_8start_CPD_Inc_noInterp.csv", index=False)



#%% run the interpolation


dt_image= 0.2
dt_cfd = 0.0001
period_length = 1.2

new_time = np.arange(0, period_length+dt_cfd, dt_cfd)    
x_new, y_new, z_new = interpolate_controlPoints_time(morph_points_df, dt_image, new_time, show=True)


#%%
work_dir = '/Users/mbarb1/OneDrive - UW/TracheaMalacia/CFD/Motion_test/Small_start/'
df_periodic_star = periodic_star_table_IncDisp_fromArraysV2(x_new, y_new, z_new, dt=dt_cfd, n_cycles=2, start_time=0)  
df_periodic_star = df_periodic_star
df_periodic_star.to_csv(work_dir +"StarControlPoints_8start_CPD_v5modSeg_Inc_fine_2skip_fullCPD.csv", index=False)
#%%
work_dir = '/Users/mbarb1/OneDrive - UW/TracheaMalacia/CFD/Motion_test/Small_start/'
df_periodic_star = periodic_star_table_totalDisp_fromArraysV2(x_new, y_new, z_new, dt=dt_cfd, n_cycles=1, start_time=0)  
df_periodic_star = df_periodic_star
df_periodic_star.to_csv(work_dir +"StarControlPoints_8start_CPD_v5modSeg_tot_Coarse_2skip.csv", index=False)

#%% save the control points in small chunks, less than a single cycle
end_time = 2.4
start_time = 0.0
dt_ct = 0.1
iter = np.arange(start_time, end_time, dt_ct)

for i in iter:
    df_periodic_star = periodic_star_table_IncDisp_fromArrays_split_newPointSet(x_new, y_new, z_new, dt=dt_cfd, start_time=i, div_per_cycle=12, cycle_len=1.2)
    df_periodic_star = df_periodic_star*1e-3
    df_periodic_star.to_csv(work_dir + "Split_CPs_fullTime_realScale/" + f"StarControlPoints_8start_CPD_v5modSeg_Inc_fine_2skip_fullCPD_{i:.3f}.csv", index=False) 





#%% load the points where the solution has reached and compare to the actua target surfaces


surf9_cfd = pv.read(work_dir + "Airway_8_t01.stl")
surf10_cfd = pv.read(work_dir + "Airway_8_t02.stl")
surf11_cfd = pv.read(work_dir + "Airway_8_t03.stl")

surf9_target = pv.read(target_surface_files[6])
surf10_target = pv.read(target_surface_files[7])
surf11_target = pv.read(target_surface_files[8])

p = pv.Plotter(shape=(3,1))
p.subplot(0,0)
p.add_mesh(surf9_cfd, opacity = 0.5, color='cyan')

p.add_mesh(surf9_target, opacity = 0.5, color='white')


p.subplot(1,0)
p.add_mesh(surf10_cfd, opacity = 0.5, color='cyan')

p.add_mesh(surf10_target, opacity = 0.5, color='white')

p.subplot(2,0)
p.add_mesh(surf11_cfd, opacity = 0.5, color='cyan')

p.add_mesh(surf11_target, opacity = 0.5, color='white')


p.show()

#%% pick points you want to track

points, cp_ids = pointPickerRemove(points_list[0], starting_surface)


#%% visualize the tracks
p = pv.Plotter()
# p.add_mesh(surf10_cfd, opacity = 0.5, color='cyan')
# p.add_mesh(surf10_target, opacity = 0.5)
p.add_mesh(starting_surface, opacity = 0.5)





for count,cp_id in enumerate(cp_ids):


    poly_fine = pv.PolyData(np.array([x_new[:,cp_id], y_new[:,cp_id], z_new[:, cp_id]]).T)
    poly_image_data = pv.PolyData(np.array([x_all[:,cp_id], y_all[:,cp_id], z_all[:, cp_id]]).T)

    p.add_mesh(poly_fine, color='red', render_points_as_spheres=True, point_size=15)
    p.add_mesh(poly_image_data, color='blue', render_points_as_spheres=True, point_size=20)
    p.add_mesh(poly_image_data.points[0], color='green', render_points_as_spheres=True, point_size=20)
    p.background_color='white'
    # p.add_mesh(surf7, opacity = 0.5)
    # p.add_mesh(surf6, opacity = 0.5)
p.show()


#%% Let's compare the morpher displacement field in star against the prescribed motion field

df_CFD = pd.read_csv(work_dir + "fine_CP_morph_disp_05.csv")

interp_id = 44
scale = 5

cp_cfd = pv.PolyData(np.array([df_CFD['X (m)'], df_CFD['Y (m)'], df_CFD['Z (m)']]).T)
cp_cfd['vec'] = np.array([df_CFD['Morpher Displacement[i] (m)'], df_CFD['Morpher Displacement[j] (m)'], 
                          df_CFD['Morpher Displacement[k] (m)']]).T * scale
cp_cfd.set_active_vectors('vec')

dx = x_new[interp_id] - x_new[interp_id - 1]
dy = y_new[interp_id] - y_new[interp_id - 1]
dz = z_new[interp_id] - z_new[interp_id - 1]

vec = np.array([dx, dy, dz])

cp_input = pv.PolyData(np.array([x_new[interp_id], y_new[interp_id], z_new[interp_id]]).T)
cp_input['vec'] = vec.T * scale
cp_input.set_active_vectors('vec')




p = pv.Plotter(shape = (1,3))
p.subplot(0,0)
p.add_mesh(cp_input)
p.add_mesh(cp_input.arrows)

p.subplot(0,1)
p.add_mesh(cp_cfd)
p.add_mesh(cp_cfd.arrows)

p.set_background('white')

cp_input['delta'] = cp_input['vec'] - cp_cfd['vec']
cp_input.set_active_vectors('delta')
p.subplot(0,2)
p.add_mesh(cp_input)
p.add_mesh(cp_input.arrows)
p.link_views()




















p.show()



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 11:28:45 2023

@author: mbarb1
"""

import pyvista as pv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


from scipy.interpolate import CubicSpline



def periodic_star_table_totalDisp_fromArrays(X,Y,Z, dt=0.1, n_cycles=5):
    """
    Create periodic displacement table for starccm. Displacement is defined as total displacement: X_n - X_0
    
    Input is X,Y,Z arrays X[time, positions]
    dt: delta time between each defined cp
    n_cycles: number of cycles to repeat for periodic simulation

    """
    
    n_time_points = len(X[:,0])
    df_save = pd.DataFrame(data=np.array([X[0,:], Y[0,:], Z[0,:]]).T, columns=["X", "Y", "Z"])
    
    column_x = "X[t={:1.3f}s]".format(0.0)
    column_y = "Y[t={:1.3f}s]".format(0.0)
    column_z = "Z[t={:1.3f}s]".format(0.0)

    df_save[column_x] = 0.0
    df_save[column_y] = 0.0
    df_save[column_z] = 0.0
    
    time_count = 1
    for n_cycle in range(n_cycles):
        for count in range(1,n_time_points):
            time =  time_count* dt
            print(time)
            
            column_x = "X[t={:1.3f}s]".format(time)
            column_y = "Y[t={:1.3f}s]".format(time)
            column_z = "Z[t={:1.3f}s]".format(time)
        
            df_save[column_x] = X[count,:] - X[0,:]
            df_save[column_y] = Y[count,:] - Y[0,:]
            df_save[column_z] = Z[count,:] - Z[0,:]

            time_count = time_count+1
        
    return df_save


def periodic_star_table_totalDisp_fromArraysV2(X,Y,Z, dt=0.1, n_cycles=5, start_time=0.0):
    """
    Create periodic displacement table for starccm. Displacement is defined as total displacement: X_n - X_0
    
    Input is X,Y,Z arrays X[time, positions]
    dt: delta time between each defined cp
    n_cycles: number of cycles to repeat for periodic simulation

    """
    
    n_time_points = len(X[:,0])
    
    n_control_points = len(X[0,:])
    
    all_data_list = []

    all_data_list.append(pd.DataFrame(data=X[0,:], columns=['X']))
    all_data_list.append(pd.DataFrame(data=Y[0,:], columns=['Y']))
    all_data_list.append(pd.DataFrame(data=Z[0,:], columns=['Z']))

    column_x = ["X[t={:1.5f}s]".format(start_time)]
    column_y = ["Y[t={:1.5f}s]".format(start_time)]
    column_z = ["Z[t={:1.5f}s]".format(start_time)]

    
    all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_x))
    all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_y))
    all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_z))

    
    time_count = 1
    for n_cycle in range(n_cycles):
        for count in range(1,n_time_points):
            time =  time_count * dt

            
            column_x = ["X[t={:1.5f}s]".format(time + start_time)]
            column_y = ["Y[t={:1.5f}s]".format(time + start_time)]
            column_z = ["Z[t={:1.5f}s]".format(time + start_time)]
            
            dx = X[count,:] - X[0,:]
            dy = Y[count,:] - Y[0,:]
            dz = Z[count,:] - Z[0,:]
            
            all_data_list.append(pd.DataFrame(data=dx, columns=column_x))
            all_data_list.append(pd.DataFrame(data=dy, columns=column_y))
            all_data_list.append(pd.DataFrame(data=dz, columns=column_z))
            
            time_count = time_count+1
        
    return pd.concat(all_data_list, axis=1)



def periodic_star_table_IncDisp_fromArrays(X,Y,Z, dt=0.1, n_cycles=5, start_time=0.0):
    """
    Create periodic displacement table for starccm. Displacement is defined as incremnetal displacement: X_n - X_n-1
    
    Input is X,Y,Z arrays X[time, positions]
    dt: delta time between each defined cp
    n_cycles: number of cycles to repeat for periodic simulation
    
    * should re-write this such that the datafram isn't growing at each itertions
    """
    
    
    n_time_points = len(X[:,0])
    df_save = pd.DataFrame(data=np.array([X[0,:], Y[0,:], Z[0,:]]).T, columns=["X", "Y", "Z"])
    
    column_x = "X[t={:1.5f}s]".format(start_time)
    column_y = "Y[t={:1.5f}s]".format(start_time)
    column_z = "Z[t={:1.5f}s]".format(start_time)

    df_save[column_x] = 0.0
    df_save[column_y] = 0.0
    df_save[column_z] = 0.0
    
    all_data_list = []
    
    time_count = 1
    for n_cycle in range(n_cycles):
        for count in range(1,n_time_points):
            time =  time_count * dt
           
            
            column_x = "X[t={:1.5f}s]".format(time + start_time)
            column_y = "Y[t={:1.5f}s]".format(time + start_time)
            column_z = "Z[t={:1.5f}s]".format(time + start_time)
            
            dx = X[count,:] - X[count-1,:]
            dy = Y[count,:] - Y[count-1,:]
            dz = Z[count,:] - Z[count-1,:]
            
            df_save[column_x] = dx
            df_save[column_y] = dy
            df_save[column_z] = dz

            time_count = time_count+1
        
    return df_save

def periodic_star_table_IncDisp_fromArraysV2(X,Y,Z, dt=0.1, n_cycles=5, start_time=0.0):
    """
    Create periodic displacement table for starccm. Displacement is defined as incremnetal displacement: X_n - X_n-1
    
    Input is X,Y,Z arrays X[time, positions]
    dt: delta time between each defined cp
    n_cycles: number of cycles to repeat for periodic simulation
    
    """
    
    
    n_time_points = len(X[:,0])
    n_control_points = len(X[0,:])
    
    all_data_list = []

    all_data_list.append(pd.DataFrame(data=X[0,:], columns=['X']))
    all_data_list.append(pd.DataFrame(data=Y[0,:], columns=['Y']))
    all_data_list.append(pd.DataFrame(data=Z[0,:], columns=['Z']))

    column_x = ["X[t={:1.5f}s]".format(start_time)]
    column_y = ["Y[t={:1.5f}s]".format(start_time)]
    column_z = ["Z[t={:1.5f}s]".format(start_time)]

    
    all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_x))
    all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_y))
    all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_z))
    
    
    time_count = 1
    for n_cycle in range(n_cycles):
        for count in range(1,n_time_points):
            time =  time_count* dt
           
    
            column_x = ["X[t={:1.5f}s]".format(time + start_time)]
            column_y = ["Y[t={:1.5f}s]".format(time + start_time)]
            column_z = ["Z[t={:1.5f}s]".format(time + start_time)]
            
            dx = X[count,:] - X[count-1,:]
            dy = Y[count,:] - Y[count-1,:]
            dz = Z[count,:] - Z[count-1,:]
            
            all_data_list.append(pd.DataFrame(data=dx, columns=column_x))
            all_data_list.append(pd.DataFrame(data=dy, columns=column_y))
            all_data_list.append(pd.DataFrame(data=dz, columns=column_z))
            

            time_count = time_count+1
        
    return pd.concat(all_data_list, axis=1)



def interpolate_controlPoints_time(df_CT, dt, new_time_array, show=True, n_demo_images=5):
    """
    Perform temporal interpolation between a coarse set of control points.
    
    interploation is a cubic spline. first and lost positions should be the same

    """
    
    n_images = len(df_CT)
    n_points = len(df_CT[0])
    
    dt = 0.1
    time_CT = np.linspace(0, (n_images-1)*dt, n_images)
    
    x_all = np.zeros((n_images, n_points))
    y_all = np.zeros((n_images, n_points))
    z_all = np.zeros((n_images, n_points))
    
    
    for count,df in enumerate(df_CT):
        x = df.X
        y = df.Y
        z = df.Z
        
        x_all[count,:] = x
        y_all[count,:] = y
        z_all[count,:] = z

    
    #Interpolate Function
    cs_x = CubicSpline(time_CT, x_all, axis=0, bc_type='periodic')
    cs_y = CubicSpline(time_CT, y_all, axis=0, bc_type='periodic')
    cs_z = CubicSpline(time_CT, z_all, axis=0, bc_type='periodic')

    x_new = cs_x(new_time_array)
    y_new = cs_y(new_time_array)
    z_new = cs_z(new_time_array)


    if show:

        p = pv.Plotter(shape=(1,n_demo_images))
        cp_ids = np.random.randint(0, n_points, n_demo_images)
        
        for count,cp_id in enumerate(cp_ids):
        
            poly_image_data = pv.PolyData(np.array([x_all[:,cp_id], y_all[:,cp_id], z_all[:, cp_id]]).T)
            poly_fine = pv.PolyData(np.array([cs_x(new_time_array)[:,cp_id], cs_y(new_time_array)[:,cp_id], cs_z(new_time_array)[:,cp_id]]).T)
            
            p.subplot(0,count)
            p.add_mesh(poly_fine, color='red', render_points_as_spheres=True, point_size=15)
            p.add_mesh(poly_image_data, color='blue', render_points_as_spheres=True, point_size=20)
            p.add_mesh(poly_image_data.points[0], color='green', render_points_as_spheres=True, point_size=20)
            p.background_color='white'
        
        p.show()
        
    return x_new, y_new, z_new




def pointPickerRemove(init_points, surf):
    
    """
    Manually remove control points using a picker
    """
        
    
    picked_points = []
    picked_ids = []

    p = pv.Plotter(notebook=0)
    p.add_mesh(surf, pickable=False,opacity = 1, color='white')
    p.add_mesh(init_points, render_points_as_spheres=True, color='blue', point_size=15, pickable=True)
    p.add_text("\n Select Points to Remove")
    def point_callback(mesh, picked_id):

        point = mesh.points[picked_id]
        picked_points.append(point)
        picked_ids.append(picked_id)
        p.add_point_labels([point,], ['You picked me!'], point_color='red', point_size=15)

    p.enable_point_picking(point_callback, use_mesh=True)
    p.show()
    
    return picked_points, picked_ids


def pointPickerAdd(init_points, surf):
    
    """
    Manually add control points using a picker
    """
        
    
    picked_points = []
    picked_ids = []

    p = pv.Plotter(notebook=0)
    p.add_mesh(surf, pickable=True, opacity=1, color='white')
    p.add_mesh(init_points, render_points_as_spheres=True, color='blue', point_size=15, pickable=False)
    p.add_text("/n Select Points to Add")


    def point_callback(mesh, picked_id):

        point = mesh.points[picked_id]
        picked_points.append(point)
        picked_ids.append(picked_id)
        p.add_point_labels([point,], ['Adding Point!'], point_color='red', point_size=15)

    p.enable_point_picking(point_callback, use_mesh=True)
    p.show()
    
    return picked_points, picked_ids

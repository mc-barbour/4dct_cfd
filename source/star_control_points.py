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


def reference_periodic_star_table_IncDisp_fromArrays_single_dt(X,Y,Z, dt=0.1, start_time=0.0, cycle_len=0.7):
    """
    Create periodic displacement table for starccm. Displacement is defined as incremnetal displacement: X_n - X_n-1
    
    Input is X,Y,Z arrays X[time, positions]
    dt: delta time between each defined cp - this is also equal to dt_cfd
    cycle_len: period length
    

    
    """
    
    n_control_points = len(X[0,:])
    
    all_data_list = []
    
    
    if start_time >= cycle_len:
        start_time_local = start_time % cycle_len
        temp_index  = int(round(start_time_local*(1/dt), 5))
        print(temp_index)
    else:
        temp_index  = int(round(start_time*(1/dt), 5))
        print(temp_index)

    all_data_list.append(pd.DataFrame(data=X[temp_index,:], columns=['X']))
    all_data_list.append(pd.DataFrame(data=Y[temp_index,:], columns=['Y']))
    all_data_list.append(pd.DataFrame(data=Z[temp_index,:], columns=['Z']))
    
    if start_time == 0.0:

        column_x = ["X[t={:1.5f}s]".format(start_time)]
        column_y = ["Y[t={:1.5f}s]".format(start_time)]
        column_z = ["Z[t={:1.5f}s]".format(start_time)]
    

        all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_x))
        all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_y))
        all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_z))
    
    

    if start_time >= cycle_len:
        start_time_local = start_time % cycle_len
        count = int((start_time_local/cycle_len)*len(X[:,0]))
    else:
        count = int((start_time/cycle_len)*len(X[:,0]))


    column_x = ["X[t={:1.5f}s]".format(dt + start_time)]
    column_y = ["Y[t={:1.5f}s]".format(dt + start_time)]
    column_z = ["Z[t={:1.5f}s]".format(dt + start_time)]
    
    dx = X[count+1,:] - X[count,:]
    dy = Y[count+1,:] - Y[count,:]
    dz = Z[count+1,:] - Z[count,:]
    
    all_data_list.append(pd.DataFrame(data=dx, columns=column_x))
    all_data_list.append(pd.DataFrame(data=dy, columns=column_y))
    all_data_list.append(pd.DataFrame(data=dz, columns=column_z))
    
    
    return pd.concat(all_data_list, axis=1)




def star_motion_table_split_incDisp(X, Y, Z, start_time=0.0, dt=0.1, dim=3, n_cycles=3, div_per_cycle=10, cycle_len=1.0):
    """
    Modified function for exporting StarCCM motion table - 1/16/2025
    This version of the function exports a sparse version of the CP where the CP locations are updated as new rows for every time-step. StarCCM does not keep 
    row correspondence with CP ids

    Input is X,Y,Z arrays X[time, positions]
    dt: delta time between each defined cp
    n_cycles: number of cycles to repeat for periodic simulation
    cycle_len: period length (s)
    Start_time: start time of this table. recomend that tables are split to only contain 0.01s of data. Otherwise files get too large.


    """
    
    n_time_points = int(len(X[:,0]) * (1/div_per_cycle)) + 1
    n_control_points = len(X[0,:])
    print(n_time_points)
    
    
    if start_time >= cycle_len:
        start_time_local = start_time % cycle_len
        start_time_index = int((start_time_local/cycle_len)*len(X[:,0]))
        print(start_time_index)

    else:
        start_time_index = int((start_time/cycle_len)*len(X[:,0]))
        print(start_time_index)
    
    
    position_list = []
    position_list.append(pd.DataFrame(data=X[start_time_index:start_time_index+n_time_points-1,:].flatten(), columns=['X']))
    position_list.append(pd.DataFrame(data=Y[start_time_index:start_time_index+n_time_points-1,:].flatten(), columns=['Y']))
    position_list.append(pd.DataFrame(data=Z[start_time_index:start_time_index+n_time_points-1,:].flatten(), columns=['Z']))
    
    df_positions = pd.concat(position_list, axis=1)
    print(df_positions)
    
    disp_table = np.zeros(((n_time_points-1) * n_control_points, n_time_points * dim))
    column_list = []
    
    column_x = "X[t={:1.5f}s]".format(start_time)
    column_y = "Y[t={:1.5f}s]".format(start_time)
    column_z = "Z[t={:1.5f}s]".format(start_time)
    
    column_list.append(column_x)
    column_list.append(column_y)
    column_list.append(column_z)
    
    
    for count in range(1, n_time_points):
        
        if start_time >= cycle_len:
            start_time_local = start_time % cycle_len
            time_count = count + int((start_time_local/cycle_len)*len(X[:,0]))
        else:
            time_count = count + int((start_time/cycle_len)*len(X[:,0]))
        time =  count * dt
        
        
        time = count * dt
        cp_range_start = count  
        
        
        column_x = "X[t={:1.5f}s]".format(time + start_time)
        column_y = "Y[t={:1.5f}s]".format(time + start_time)
        column_z = "Z[t={:1.5f}s]".format(time + start_time)
        
        column_list.append(column_x)
        column_list.append(column_y)
        column_list.append(column_z)
        
        # calculate displacement
        dx = X[time_count,:] - X[time_count-1,:]
        dy = Y[time_count,:] - Y[time_count-1,:]
        dz = Z[time_count,:] - Z[time_count-1,:]
        
        # populate displacement table
        disp_table[(count-1)*n_control_points:(count-1)*n_control_points + n_control_points, count * dim] = dx
        disp_table[(count-1)*n_control_points:(count-1)*n_control_points + n_control_points, count * dim + 1] = dy
        disp_table[(count-1)*n_control_points:(count-1)*n_control_points + n_control_points, count * dim + 2] = dz
    
        
        
    df_motion = pd.DataFrame(disp_table, columns = column_list)
    
    df_full = pd.concat([df_positions, df_motion], axis=1)
    
    
    return df_full




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

def periodic_star_table_IncDisp_fromArrays_split_newPointSet(X,Y,Z, dt=0.1, start_time=0.0, div_per_cycle = 7, cycle_len=0.7):
    """
    Create periodic displacement table for starccm. Displacement is defined as incremnetal displacement: X_n - X_n-1
    
    Input is X,Y,Z arrays X[time, positions]
    dt: delta time between each defined cp
    n_cycles: number of cycles to repeat for periodic simulation
    
    """
    
    # if start_time % cycle_len == 0:
    #     n_time_points = int(len(X[:,0]) * (1/div)) + 1
    # else:
    #     n_time_points = int(len(X[:,0]) * (1/div))
    
    n_time_points = int(len(X[:,0]) * (1/div_per_cycle)) + 1

    n_control_points = len(X[0,:])
    
    if start_time >= cycle_len:
        start_time_local = start_time % cycle_len
        start_index = int((start_time_local/cycle_len)*len(X[:,0]))

    else:
        start_index = int((start_time/cycle_len)*len(X[:,0]))
    
    
    all_data_list = []

    all_data_list.append(pd.DataFrame(data=X[start_index,:], columns=['X']))
    all_data_list.append(pd.DataFrame(data=Y[start_index,:], columns=['Y']))
    all_data_list.append(pd.DataFrame(data=Z[start_index,:], columns=['Z']))

    column_x = ["X[t={:1.5f}s]".format(start_time)]
    column_y = ["Y[t={:1.5f}s]".format(start_time)]
    column_z = ["Z[t={:1.5f}s]".format(start_time)]

    if start_time % cycle_len == 0:
        all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_x))
        all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_y))
        all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_z))
    
    
    time_count = 1

    for count in range(1, n_time_points):
        
        if start_time >= cycle_len:
            start_time_local = start_time % cycle_len
            count = count + int((start_time_local/cycle_len)*len(X[:,0]))
        else:
            count = count + int((start_time/cycle_len)*len(X[:,0]))
        time =  time_count * dt
       

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




def periodic_star_table_IncDisp_fromArrays_split(X,Y,Z, dt=0.1, n_cycles=5, start_time=0.0, div=7, cycle_len=0.7):
    """
    Create periodic displacement table for starccm. Displacement is defined as incremnetal displacement: X_n - X_n-1
    
    Input is X,Y,Z arrays X[time, positions]
    dt: delta time between each defined cp
    n_cycles: number of cycles to repeat for periodic simulation
    
    """
    
    # if start_time % cycle_len == 0:
    #     n_time_points = int(len(X[:,0]) * (1/div)) + 1
    # else:
    #     n_time_points = int(len(X[:,0]) * (1/div))
    
    n_time_points = int(len(X[:,0]) * (1/div)) + 1

    n_control_points = len(X[0,:])
    
    
    all_data_list = []

    all_data_list.append(pd.DataFrame(data=X[0,:], columns=['X']))
    all_data_list.append(pd.DataFrame(data=Y[0,:], columns=['Y']))
    all_data_list.append(pd.DataFrame(data=Z[0,:], columns=['Z']))

    column_x = ["X[t={:1.5f}s]".format(start_time)]
    column_y = ["Y[t={:1.5f}s]".format(start_time)]
    column_z = ["Z[t={:1.5f}s]".format(start_time)]

    if start_time % cycle_len == 0:
        all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_x))
        all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_y))
        all_data_list.append(pd.DataFrame(data=np.zeros(n_control_points), columns=column_z))
    
    
    time_count = 1
    for n_cycle in range(n_cycles):
        for count in range(1,n_time_points):
            if start_time >= cycle_len:
                start_time_local = start_time % cycle_len
                count = count + int((start_time_local/cycle_len)*len(X[:,0]))
            else:
                count = count + int((start_time/cycle_len)*len(X[:,0]))
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

def periodic_star_table_TargetDisp_fromArraysV2(X,Y,Z, dt=0.1, n_cycles=5, start_time=0.0):
    """
    Create periodic displacement table for starccm. Displacement is defined as target position: X_n - X_n-1
    
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

    
    all_data_list.append(pd.DataFrame(data=X[0,:], columns=column_x))
    all_data_list.append(pd.DataFrame(data=Y[0,:], columns=column_y))
    all_data_list.append(pd.DataFrame(data=Z[0,:], columns=column_z))
    
    
    time_count = 1
    for n_cycle in range(n_cycles):
        for count in range(1,n_time_points):
            time =  time_count * dt
           
    
            column_x = ["X[t={:1.5f}s]".format(time + start_time)]
            column_y = ["Y[t={:1.5f}s]".format(time + start_time)]
            column_z = ["Z[t={:1.5f}s]".format(time + start_time)]
            
            
            all_data_list.append(pd.DataFrame(data=X[time_count,:], columns=column_x))
            all_data_list.append(pd.DataFrame(data=Y[time_count,:], columns=column_y))
            all_data_list.append(pd.DataFrame(data=Z[time_count,:], columns=column_z))
            

            time_count = time_count+1
        
    return pd.concat(all_data_list, axis=1)


def interpolate_controlPoints_time(df_CT, dt, new_time_array, show=True, n_demo_images=5):
    """
    Perform temporal interpolation between a coarse set of control points.
    
    interploation is a cubic spline. first and lost positions should be the same

    """
    
    n_images = len(df_CT)
    n_points = len(df_CT[0])
    

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

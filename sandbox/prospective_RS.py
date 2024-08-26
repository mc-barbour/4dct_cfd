#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:40:09 2023

@author: mbarbour

"""

import pyvista as pv
import numpy as np
import pandas as pd
import vtk
import glob
# import vedo as vd

import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots

from paths import *
from postprocess.process_utils import *
from postprocess.loader import *

pv.set_plot_theme('document')
pio.renderers.default = "browser"

colors = [px.colors.qualitative.Dark24[23],  px.colors.qualitative.Dark24[22],
          px.colors.qualitative.Dark2[7], px.colors.qualitative.Vivid[0]]



#%% working directory

data_dir = '/Users/mbarbour/OneDrive - UW/RobinSequence/Data/CFD/Subj01/Inspiration/'

vol_, wall_, inlet_ = load_casedata_to_poly(data_dir + 'Subj01_peakInspiration/')
cl_file = glob.glob(data_dir + "*CL.vtp")
print(cl_file)
cl = pv.read(cl_file[0])


df = centerline_metric_comp(cl, vol_)


#%% process all the inhalation data
data_dir = '/Users/mbarbour/OneDrive - UW/RobinSequence/Data/CFD/'

conditions = ['Inspiration', 'Exhalation']

patient_names = ['Subj01', 'Subj02', 'Subj03']

scale = [False, True, True]
regions = []
walls = []
inlets = []
centerlines = []

for count, patient in enumerate(patient_names):
    case_dir = glob.glob(data_dir + patient + '/Inspiration/CaseFiles/')
    vol_, wall_, inlet_ = load_casedata_to_poly(case_dir[0])
    cl_file = glob.glob(data_dir + patient + '/Inspiration/*CL.vtp')
    print(cl_file)
    cl = pv.read(cl_file[0])
    
    if scale[count]:
        cl.scale(1e-3)

    
    regions.append(vol_)
    walls.append(wall_)
    inlets.append(inlet_)
    centerlines.append(cl)
    
    
 #%% Plot pressure maps
p = pv.Plotter(shape=(1,3))

for count, (wall_,cl) in enumerate(zip(walls, centerlines)):
    
    p.subplot(0,count)
    
    p.add_mesh(wall_, scalars='Pressure', clim=[-20, 0],cmap='inferno', opacity=0.7)
    p.add_mesh(cl)
    p.add_mesh(pv.PolyData(cl.points[0,:]), color='red', point_size=20, label='start')
    p.add_mesh(pv.PolyData(cl.points[-1,:]), color='yellow', point_size=20, label='end')
    p.add_text(patient_names[count])
p.add_legend()
p.show()


#%% flip centerlines - want centerline to be in the direction of flow

for count,cl in enumerate(centerlines):
    centerlines[count] = flip_centerlines(cl)
    
    
    
#%% main function exectution
dfs_ex = []
for count, (region,cl) in enumerate(zip(regions, centerlines)):
    print("Computing flow metrics for condition: {}".format(patient_names[count]))
    dfs_ex.append(centerline_metric_comp(cl, region))

#%%

for (patient,df) in zip(patient_names,dfs):
    save_dir = data_dir + patient
    print(save_dir)
    df.to_csv(save_dir + 'CenterlineMetrics_peakInhalation.csv')
    
#%% save exhalation metrics
for (patient,df) in zip(patient_names,dfs_ex):
    save_dir = data_dir + patient
    print(save_dir)
    df.to_csv(save_dir + patient + 'CenterlineMetrics_peakExpiration.csv')
    

#%% load inhalation metrics


dfs_in = []
for patient in patient_names:
    patient_dir = data_dir + patient
    df = pd.read_csv(patient_dir +'/Inspiration/' + patient +'CenterlineMetrics_peakInhalation.csv')
    dfs_in.append(df)

#%% plot csa - Inhalation
fig = go.Figure()

for count, (condition, df) in enumerate(zip(patient_names, dfs_in)):
    fig.add_trace(go.Scatter(x = df['Centerline Distance (m)'], y = df['Cross Sectional Area (mm2)'], name=condition, marker_color=colors[count]))

fig.update_layout(template='plotly_white')
fig.show()

#%% plot resistance - Inhalation
fig = go.Figure()

for count, (condition, df) in enumerate(zip(patient_names, dfs_in)):
    fig.add_trace(go.Scatter(x = df['Centerline Distance (m)'], y = df['Average Pressure (Pa)'], name=condition, marker_color=colors[count]))

fig.update_layout(template='plotly_white')
fig.show()


#%% plot csa - Exhalation
fig = go.Figure()

for count, (condition, df) in enumerate(zip(patient_names, dfs_ex)):
    fig.add_trace(go.Scatter(x = df['Centerline Distance (m)'], y = df['Cross Sectional Area (mm2)'], name=condition, marker_color=colors[count]))

fig.update_layout(template='plotly_white')
fig.show()


#%% create bar charts


oneDrive_data = '/Users/mbarbour/OneDrive - UW/RobinSequence/Data/'
df = pd.read_excel(oneDrive_data + 'RS Simulation Parameters.xlsx', sheet_name='Simulation Results')

fig_dir = '/Users/mbarbour/OneDrive - UW/RS Flow Modeling/Manuscripts/Engineering_RS_Baseline/Figures/'
#%% re-calculate the centerline data for Subj01 -using smoothing cl

patient_dir = '/Users/mbarbour/OneDrive - UW/RobinSequence/Data/CFD/Subj01/Inspiration/'


cl = pv.read(glob.glob(patient_dir + "*CL2.vtp")[0])
cl.scale(1e-3, inplace=True)

vol, wall, inlet = load_casedata_to_poly(patient_dir + "CaseFile_correctFlow/")

# you might want to flip the orientation of your centerlines file
flip_centerline = True

if flip_centerline:
    temp = cl.points
    cl_new = np.flip(temp, axis=0)
    cl = pv.PolyData(cl_new)

    # confirm the change
    p1 = pv.PolyData(cl.points[0])
    pend = pv.PolyData(cl.points[-1])
    p = pv.Plotter()
    p.add_mesh(wall, opacity=0.5, scalars='Pressure')
    p.add_mesh(cl)
    p.add_mesh(p1, point_size=20, color='red', label='centerline Start')
    p.add_mesh(pend, point_size=20, color='white', label='centerline End')
    p.add_text("End Point: {:d}".format(len(cl.points)))
    p.add_legend()

    p.show()
    
df = centerline_metric_comp(cl, vol)
df.to_csv(patient_dir +  'Subj01CenterlineMetrics_peakInspiration_v2.csv')
#%%

df["Minimum CSA (mm2)"] = df["Minimum CSA (m2)"]*1e6

inhalation_index = [0,2,4,6,8,10,12,14]
exhalation_index = [1,3,5,7,9,11,13,15]

df_in = df.loc[inhalation_index]
df_ex = df.loc[exhalation_index]


#%%
inhalation_index = [6, 2, 8, 14, 4, 10, 12, 0]
exhalation_index = [7, 3, 9, 15, 5, 11, 13, 1]



plot_variables = ['Resistance (Pa/(mL/s))', 'Total Dissipation (mW)', 'Peak Velocity (m/s)', 'Minimum CSA (mm2)']
plot_variables_latex = ['Resistance <br> (Pa/(mL/s))', 'Energy Dissipation <br> (mW)', 'Peak Velocity <br> (m/s)', 'Minimum CSA <br> (mm2)']
# patient_list = ["Pt" + str(int(a)) for a in df.Patient[inhalation_index].values]

patient_letters = ["E", "C", "A", "D", "B"] # 109, 106, 101, 107, 103
patient_letters = ["A", "B", "C", "D", "E", "F", "G", "H"] # 109, 106, 101, 107, 103

patient_list_manuscript = ["Patient " + a for a in patient_letters]

fig = go.Figure()

ncols=2
nrows=2

fig = make_subplots(rows=2, cols=2)

for count, variable in enumerate(plot_variables):
    
    row = int(count/ncols) + 1
    col = count%(ncols) +1
    legendflag = False
    if count == 0:
        legendflag = True
        
    fig.add_trace(go.Bar(y=df.loc[inhalation_index,variable], marker_color=colors[0], showlegend=legendflag, name='Inhalation'), row=row, col=col)
    fig.add_trace(go.Bar(y=df.loc[exhalation_index,variable], marker_color=colors[1], showlegend=legendflag, name='Exhalation'), row=row, col=col)
    fig.update_yaxes(title=plot_variables_latex[count], row=row, col=col)
    fig.update_xaxes(title='Patient', ticktext = patient_letters, tickvals=[0,1,2,3,4,5,6,7])
    fig.update_layout(template='plotly_white')

fig.update_layout(legend=dict(
   orientation='h',
   yanchor='top',
   y=1,
   xanchor='right',
   x=1),
    font_family="arial",
    font_size=20,
    font_color='black')
fig.write_image(oneDrive_data + 'RS_CFDMetrics_071123.png', width=1300, height=700, scale=2)
fig.show()



#%% create flow vis for Subj02
patient_dir = '/Users/mbarbour/OneDrive - UW/RobinSequence/Data/CFD/Subj02/Inspiration/'


cl = pv.read(glob.glob(patient_dir + "*CL.vtp")[0])
cl.scale(1e-3, inplace=True)

vol, wall, inlet = load_casedata_to_poly(patient_dir + "CaseFiles/")

# you might want to flip the orientation of your centerlines file
flip_centerline = True

if flip_centerline:
    temp = cl.points
    cl_new = np.flip(temp, axis=0)
    cl = pv.PolyData(cl_new)

    # confirm the change
    p1 = pv.PolyData(cl.points[0])
    pend = pv.PolyData(cl.points[-1])
    p = pv.Plotter()
    p.add_mesh(wall, opacity=0.5, scalars='Pressure')
    p.add_mesh(cl)
    p.add_mesh(p1, point_size=20, color='red', label='centerline Start')
    p.add_mesh(pend, point_size=20, color='white', label='centerline End')
    p.add_text("End Point: {:d}".format(len(cl.points)))
    p.add_legend()

    p.show()
    


# extract slices and plot velocity in pyvista

p = pv.Plotter()
p.add_mesh(wall, opacity=0.2, color='white')
p.background_color = 'w'
pv.global_theme.font.color = 'black'

#pmin = pv.PolyData(cl.points[187]) 
sargs = dict(interactive=True)
n_points = len(cl.points)
start_point = 70 # starting centerline integer
delta_points = 25 # centerline locations to skip
for point_id in range(start_point, n_points - 50, delta_points):  
    print(point_id)

    cut_slice, n, x0 = plane_extract(vol, cl, point_id)
    connected = cut_slice.connectivity()
    surface_count = len(np.unique(connected['RegionId']))

    print('N surfaces: {:d}'.format(surface_count))
    if surface_count >= 2:
        cut_slice = split_surface(connected, cl.points[point_id])

    # p.add_mesh(cut_slice, scalars='VelocityMagnitude',clim=[0,20], cmap='cet_linear_kbgyw_5_98_c62')
    # p.add_mesh(cut_slice, scalars='VelocityMagnitude',clim=[0,20], cmap='cet_bgy', scalar_bar_args=sargs)
    p.add_mesh(cut_slice, scalars='VelocityMagnitude',clim=[0,12], cmap='cet_bmy', scalar_bar_args=sargs)
    # p.add_mesh(cut_slice, scalars='VelocityMagnitude',clim=[0,25], cmap='viridis', smooth_shading=True)
#p.add_mesh(pmin, render_points_as_spheres=True, color='cyan', point_size=45)
p.show()


#%% centerline values 
patient_flowrate = 130
df_subj02 = pd.read_csv(patient_dir + 'Subj02CenterlineMetrics_peakInhalation.csv')
df_subj02['Local Resistance (Pa/m/(mL/s)'] = np.gradient(df_subj02['Average Pressure (Pa)'], df_subj02['Centerline Distance (m)']) / patient_flowrate
df_subj02['Energy Loss'] = abs(df_subj02['Energy Flux (W)'] - df_subj02['Energy Flux (W)'][0])

resistance = (df_subj02['Average Pressure (Pa)'][0] - df_subj02['Average Pressure (Pa)']) / patient_flowrate
df_subj02['Resistance'] = resistance

fig = make_subplots(rows=2, cols=1, specs=[
                    [{"secondary_y": True}], [{"secondary_y": True}]])
ticktext = [0.04, 0.06, 0.08, 0.1, ]
tickvals = [0.04, 0.06, 0.08, 0.1, 0.12, 0.14]


# Bottom subplot
plot_row = 1
fig.add_trace(go.Scatter(x=df_subj02['Centerline Distance (m)'],
              y=df_subj02['Cross Sectional Area (mm2)']*1e6, marker_color=colors[2], line_width=3), row=plot_row, col=1)
fig.add_trace(go.Scatter(x=df_subj02['Centerline Distance (m)'],
              y=df_subj02['Local Resistance (Pa/m/(mL/s)']*-1, marker_color=colors[3], line_width=4), row=plot_row, col=1, secondary_y=True)

fig.update_xaxes(range=[0.06, 0.13], row=plot_row, col=1, tickangle=-90,
                 gridwidth=1, showline=True, gridcolor='lightgrey', linecolor='lightgrey',
                 ticktext=ticktext, tickvals=tickvals)
fig.update_yaxes(range=[0, 100], title='Area (mm2)', row=plot_row,
                 col=1, secondary_y=False, tickangle=-90, gridwidth=1, showline=True, gridcolor='lightgrey', linecolor='lightgrey')

fig.update_yaxes(title='Local Resistance <br> (Pa/m/(mL/s))',
                 row=plot_row, col=1, color=colors[3], secondary_y=True, tickangle=-90, showgrid=False, showline=True, linecolor='lightgrey')


# Top Subplot
plot_row = 2
fig.add_trace(go.Scatter(x=df_subj02['Centerline Distance (m)'], y=df_subj02['Energy Loss']*1e3,
                         marker_color=colors[0],
                         line_width=4), row=plot_row, col=1, secondary_y=True)


fig.add_trace(go.Scatter(x=df_subj02['Centerline Distance (m)'], y=(df_subj02['Resistance']),
                         marker_color=colors[1],
                         line_width=4), row=plot_row, col=1)

fig.update_yaxes(title='Resistance <br> (Pa/(mL/s))',
                 color=colors[1], row=plot_row, col=1, secondary_y=False, tickangle=-90, gridwidth=1, showline=True, gridcolor='lightgrey')
fig.update_yaxes(title='Accumulated <br> Energy Loss (mW)',
                 color=colors[0], row=plot_row, col=1, secondary_y=True, tickangle=-90, showgrid=False, showline=True)

fig.update_xaxes(range=[0.06, 0.13], row=plot_row, col=1, tickangle=-90,
                 gridwidth=1, gridcolor='lightgrey', showline=True, linecolor='lightgrey',
                 ticktext=ticktext, tickvals=tickvals)



fig.update_layout(template='plotly_white', showlegend=False, font_size=16)
fig.update_xaxes(title='Airway Centerline Distance (m)', row=2, col=1)


fig.write_image(patient_dir + 'subj02_local_metrics_v3.png', scale=2)

fig.show()

#%% create flow vis for Subj01
patient_dir = '/Users/mbarbour/OneDrive - UW/RobinSequence/Data/CFD/Subj01/Inspiration/'


cl = pv.read(glob.glob(patient_dir + "*CL.vtp")[0])
# cl.scale(1e-3, inplace=True)

vol, wall, inlet = load_casedata_to_poly(patient_dir + "CaseFiles_correctFlow/")

# you might want to flip the orientation of your centerlines file
flip_centerline = True

if flip_centerline:
    temp = cl.points
    cl_new = np.flip(temp, axis=0)
    cl = pv.PolyData(cl_new)

    # confirm the change
    p1 = pv.PolyData(cl.points[0])
    pend = pv.PolyData(cl.points[-1])
    p = pv.Plotter()
    p.add_mesh(wall, opacity=0.5, scalars='Pressure')
    p.add_mesh(cl)
    p.add_mesh(p1, point_size=20, color='red', label='centerline Start')
    p.add_mesh(pend, point_size=20, color='white', label='centerline End')
    p.add_text("End Point: {:d}".format(len(cl.points)))
    p.add_legend()

    p.show()
    

#%%
# extract slices and plot velocity in pyvista

p = pv.Plotter()
p.add_mesh(wall, opacity=0.2, color='white')
p.background_color = 'w'
pv.global_theme.font.color = 'black'

#pmin = pv.PolyData(cl.points[187]) 
sargs = dict(interactive=True)
n_points = len(cl.points)
start_point = 50 # starting centerline integer
delta_points = 12 # centerline locations to skip
for point_id in range(start_point, n_points - 50, delta_points):  
    print(point_id)

    cut_slice, n, x0 = plane_extract(vol, cl, point_id)
    connected = cut_slice.connectivity()
    surface_count = len(np.unique(connected['RegionId']))

    print('N surfaces: {:d}'.format(surface_count))
    if surface_count >= 2:
        cut_slice = split_surface(connected, cl.points[point_id])

    # p.add_mesh(cut_slice, scalars='VelocityMagnitude',clim=[0,20], cmap='cet_linear_kbgyw_5_98_c62')
    # p.add_mesh(cut_slice, scalars='VelocityMagnitude',clim=[0,20], cmap='cet_bgy', scalar_bar_args=sargs)
    p.add_mesh(cut_slice, scalars='VelocityMagnitude',clim=[0,4], cmap='cet_bmy', scalar_bar_args=sargs)
    # p.add_mesh(cut_slice, scalars='VelocityMagnitude',clim=[0,25], cmap='viridis', smooth_shading=True)
#p.add_mesh(pmin, render_points_as_spheres=True, color='cyan', point_size=45)
p.show()

#%% centerline values - subj01
patient_flowrate = 78
df_subj01 = pd.read_csv(patient_dir + 'Subj01CenterlineMetrics_peakInspiration_v2.csv')
df_subj01['Local Resistance (Pa/m/(mL/s)'] = np.gradient(df_subj01['Average Pressure (Pa)'], df_subj01['Centerline Distance (m)']) / patient_flowrate
df_subj01['Energy Loss'] = abs(df_subj01['Energy Flux (W)'] - df_subj01['Energy Flux (W)'][0])

resistance = (df_subj01['Average Pressure (Pa)'][0] - df_subj01['Average Pressure (Pa)']) / patient_flowrate
df_subj01['Resistance'] = resistance

fig = make_subplots(rows=2, cols=1, specs=[
                    [{"secondary_y": True}], [{"secondary_y": True}]])
ticktext = [0.04, 0.06, 0.08, 0.1, ]
tickvals = [0.04, 0.06, 0.08, 0.1, 0.12, 0.14]


kernel_size = 10
kernel = np.ones(kernel_size) / kernel_size
local_resistance_smooth = np.convolve(df_subj01['Local Resistance (Pa/m/(mL/s)'], kernel, mode='same')


# Bottom subplot
plot_row = 1
fig.add_trace(go.Scatter(x=df_subj01['Centerline Distance (m)'],
              y=df_subj01['Cross Sectional Area (mm2)']*1e6, marker_color=colors[2], line_width=3), row=plot_row, col=1)
fig.add_trace(go.Scatter(x=df_subj01['Centerline Distance (m)'],
              y=local_resistance_smooth*-1, marker_color=colors[3], line_width=4), row=plot_row, col=1, secondary_y=True)

fig.update_xaxes(range=[0.07, 0.14], row=plot_row, col=1, tickangle=-90,
                 gridwidth=1, showline=True, gridcolor='lightgrey', linecolor='lightgrey',
                 ticktext=ticktext, tickvals=tickvals)
fig.update_yaxes(range=[0, 250], title='Area (mm2)', row=plot_row,
                 col=1, secondary_y=False, tickangle=-90, gridwidth=1, showline=True, gridcolor='lightgrey', linecolor='lightgrey')

fig.update_yaxes(title='Local Resistance <br> (Pa/m/(mL/s))',
                 row=plot_row, col=1, color=colors[3], secondary_y=True, tickangle=-90, showgrid=False, showline=True, linecolor='lightgrey')


# Top Subplot
plot_row = 2
fig.add_trace(go.Scatter(x=df_subj01['Centerline Distance (m)'], y=df_subj01['Energy Loss']*1e3,
                         marker_color=colors[0],
                         line_width=4), row=plot_row, col=1, secondary_y=True)


fig.add_trace(go.Scatter(x=df_subj01['Centerline Distance (m)'], y=(df_subj01['Resistance']),
                         marker_color=colors[1],
                         line_width=4), row=plot_row, col=1)

fig.update_yaxes(title='Resistance <br> (Pa/(mL/s))',
                 color=colors[1], row=plot_row, col=1, secondary_y=False, tickangle=-90, gridwidth=1, showline=True, gridcolor='lightgrey')
fig.update_yaxes(title='Accumulated <br> Energy Loss (mW)',
                 color=colors[0], row=plot_row, col=1, secondary_y=True, tickangle=-90, showgrid=False, showline=True)

fig.update_xaxes(range=[0.07, 0.14], row=plot_row, col=1, tickangle=-90,
                 gridwidth=1, gridcolor='lightgrey', showline=True, linecolor='lightgrey',
                 ticktext=ticktext, tickvals=tickvals)



fig.update_layout(template='plotly_white', showlegend=False, font_size=16)
fig.update_xaxes(title='Airway Centerline Distance (m)', row=2, col=1)


fig.write_image(patient_dir + 'subj01_local_metrics_v3.png', scale=2)

fig.show()

#%% plot three peak resistance values

patient_flowrates = [78, 130, 70]
fig = make_subplots(rows=3, cols=1)
max_indices = []

for count, df in enumerate(dfs_in):
    df['Local Resistance (Pa/m/(mL/s)'] = np.gradient(df['Average Pressure (Pa)'], df['Centerline Distance (m)']) / patient_flowrates[count]

    kernel_size = 10
    kernel = np.ones(kernel_size) / kernel_size
    local_resistance_smooth = -1*np.convolve(df['Local Resistance (Pa/m/(mL/s)'], kernel, mode='same')
    
    max_res_index = np.argmax(abs(local_resistance_smooth))

    fig.add_trace(go.Scatter(x=df["Centerline Distance (m)"], y=local_resistance_smooth, name=names[count], marker_color=colors[count]), row=count+1, col=1)
    fig.add_trace(go.Scatter(x=[df["Centerline Distance (m)"][max_res_index]], y=[local_resistance_smooth[max_res_index]], mode='markers', marker_color='red', marker_size=15, showlegend=False), row=count+1, col=1)
    
    max_indices.append(max_res_index)
    fig.update_xaxes(title='Centerline Distance (m)')
    fig.update_yaxes(title='Local Resistance <br> (Pa/m/(mL/s))')
fig.update_layout(template='plotly_white')

fig.show()


#%% surfaces

names = ["Subj01, H", "Subj02, B", "Subj03, E"]

p = pv.Plotter(shape=(1,3))

for count, surf in enumerate(walls):
    p.subplot(0,count)
    p.add_mesh(surf, opacity=0.5, color='white')

    point = [dfs_in[count].x0[max_indices[count]], dfs_in[count].y0[max_indices[count]], dfs_in[count].z0[max_indices[count]]]

    p.add_mesh(pv.PolyData(point), render_points_as_spheres=True, point_size=20, color='red')
    p.add_text(names[count])
p.show()
                

#%% repeat the same analysis but for exhalation


dfs_ex = []
for patient in patient_names:
    patient_dir = data_dir + patient
    df = pd.read_csv(patient_dir +'/Expiration/' + patient +'CenterlineMetrics_peakExpiration.csv')
    dfs_ex.append(df)




data_dir = '/Users/mbarbour/OneDrive - UW/RobinSequence/Data/CFD/'


patient_names = ['Subj01', 'Subj02', 'Subj03']

scale = [False, True, True]
regions = []
walls = []
inlets = []


for count, patient in enumerate(patient_names):
    case_dir = glob.glob(data_dir + patient + '/Expiration/CaseFiles/')
    vol_, wall_, inlet_ = load_casedata_to_poly(case_dir[0])


    
    regions.append(vol_)
    walls.append(wall_)
    inlets.append(inlet_)




#%%

max_indices = []
patient_flowrates = [78, 130, 87]
fig = make_subplots(rows=3, cols=1)


for count, df in enumerate(dfs_ex):
    df['Local Resistance (Pa/m/(mL/s)'] = np.gradient(df['Average Pressure (Pa)'], df['Centerline Distance (m)']) / patient_flowrates[count]

    kernel_size = 10
    kernel = np.ones(kernel_size) / kernel_size
    local_resistance_smooth = -1*np.convolve(df['Local Resistance (Pa/m/(mL/s)'], kernel, mode='same')
    
    max_res_index = np.argmax(abs(local_resistance_smooth))

    # fig.add_trace(go.Scatter(x=df["Centerline Distance (m)"], y=local_resistance_smooth, name=names[count], marker_color=colors[count]), row=count+1, col=1)
    # fig.add_trace(go.Scatter(x=[df["Centerline Distance (m)"][max_res_index]], y=[local_resistance_smooth[max_res_index]], mode='markers', marker_color='red', marker_size=15, showlegend=False), row=count+1, col=1)
    
    fig.add_trace(go.Scatter( y=local_resistance_smooth, name=names[count], marker_color=colors[count]), row=count+1, col=1)
    fig.add_trace(go.Scatter(x=[max_indices[count]], y=[local_resistance_smooth[max_indices[count]]], mode='markers', marker_color='red', marker_size=15, showlegend=False), row=count+1, col=1)
    
    
    
    max_indices.append(max_res_index)
    fig.update_xaxes(title='Centerline Distance (m)')
    fig.update_yaxes(title='Local Resistance <br> (Pa/m/(mL/s))')
fig.update_layout(template='plotly_white')

fig.show()


names = ["Subj01, H", "Subj02, B", "Subj03, E"]

p = pv.Plotter(shape=(1,3))

for count, surf in enumerate(walls):
    p.subplot(0,count)
    p.add_mesh(surf, opacity=0.5)

    point = [dfs_ex[count].x0[max_indices[count]], dfs_ex[count].y0[max_indices[count]], dfs_ex[count].z0[max_indices[count]]]

    p.add_mesh(pv.PolyData(point), render_points_as_spheres=True, point_size=20, color='red')
    p.add_text(names[count])
p.show()
           




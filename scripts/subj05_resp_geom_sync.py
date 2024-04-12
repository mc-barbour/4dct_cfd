# -*- coding: utf-8 -*-
"""
Created on Mon May  8 10:42:23 2023

@author: barbourm
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import pandas as pd
import glob
from natsort import natsorted
import pyvista as pv
from scipy import signal

import plotly.io as pio
from plotly.subplots import make_subplots

from source.paths import *
from source.respiratory_functions import *

pio.renderers.default = "browser"


pio.renderers.default = "browser"
colors = [px.colors.qualitative.Dark24[23],  px.colors.qualitative.Dark24[22],
          px.colors.qualitative.Dark2[7], px.colors.qualitative.Vivid[0]]



#%%  Load segmentation files
data_dir = onedrive_dir / 'Area and Volume Data/4DCT_analysis/Subj05/'

seg_files = natsorted( (data_dir / "STL_Sequence/").glob("*.stl"))


#%% compute defined region for airway volume slice

surf = pv.read(seg_files[0])

voxels = pv.voxelize(surf, density=surf.length / 200)

y_slice_location = -1.65e2
p = pv.Plotter()
p.add_mesh(voxels.clip(normal="y", origin=(0, y_slice_location,0), invert=False), color=True, show_edges=True, opacity=0.5)
p.add_mesh(surf, color="red", opacity=0.5)
p.add_axes()
p.show()


#%% compute airway volume for all segmentations

volume = []
for surf_file in seg_files:
    surf = pv.read(surf_file)
    voxels = pv.voxelize(surf, density=surf.length / 200)
    clip = voxels.clip(normal="y", origin=(0, y_slice_location, 0), invert=False)
    volume.append(abs(clip.volume))
    print(clip.volume)
    
    
fig = go.Figure()
fig.add_trace(go.Scatter(y=volume))
fig.show()

np.savetxt(data_dir / "STL_airway_volumes.txt", volume)
    

#%% load in the respiratory phase data


df = pd.read_excel(data_dir / "S05_VideoRespMotionData.xlsx")

chest_motion = df['Chest Motion Curve (pixels, zero mean, increasing means chest rise)']
vid_time = df["Time (s)"]
ct_time = df['CT Frame Times (s)']

#%% need to calculate binned ct data
# aquisition_time = .25

# ct_time = df["CT Frame Times (s)"][0:38]

# binned_chest_motion = []

# for time in ct_time:
#     start = time - aquisition_time / 2
#     end = time + aquisition_time / 2
    
#     start_idx = np.argmin(abs(df["Time (s)"] - start))
#     end_idx = np.argmin(abs(df["Time (s)"] - end))
    
#     avg_chest_motion = np.mean(df["Chest Motion Curve (pixels, zero mean, increasing means chest rise)"][start_idx:end_idx])
#     print(df["Time (s)"][end_idx] - df["Time (s)"][start_idx], avg_chest_motion)
#     binned_chest_motion.append(avg_chest_motion)


#%% smooth the chestmotion data and compute gradient

smooth_chest = signal.savgol_filter(chest_motion, window_length=43, polyorder=3)
grad_chest = np.gradient(smooth_chest)


#%% plot chest motion and volume curves

scale = 10
fig = make_subplots(rows=2, cols=1)

fig.add_trace(go.Scatter(x=vid_time, y=chest_motion, name = "Raw chest motion"), row=1, col=1)
fig.add_trace(go.Scatter(x=vid_time, y=smooth_chest, name = "Raw chest motion - smooth"), row=1, col=1)
fig.add_trace(go.Scatter(x=vid_time, y=grad_chest*scale, name = "chest motion gradient"), row=1, col=1)

fig.update_yaxes(title='Chest Motion (pixel)', row=1, col=1)
fig.update_xaxes(title='Time (s)', range=[0,4.1], row=1, col=1)


fig.add_trace(go.Scatter(x=ct_time, y=volume), row=2, col=1)
fig.update_yaxes(title='Volume (mm3)', row=2, col=1)
fig.update_xaxes(title='Time (s)', range=[0,4.1], row=2, col=1)

fig.show()






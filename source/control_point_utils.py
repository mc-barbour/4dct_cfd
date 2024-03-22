#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 12:27:26 2022

@author: mbarbour
"""

import numpy as np
import pyvista as pv
import pandas as pd
import json
import os
import csv
import vtk
import glob
from matplotlib import cm
from scipy.spatial import cKDTree as KDTree



colors = cm.jet(np.linspace(0, 1, 6))


"""
Functions for defining control point placement and modifying control point tracks


"""


def intersecting_cells(vox, surf):
    
    cell_index = 2
    single_cell = vox.extract_cells(cell_index)
    cell_length = abs(single_cell.bounds[0] - single_cell.bounds[1])
    test_length = np.sqrt(3)*cell_length/2
    
    voxel_centers = vox.cell_centers()
    voxel_centers.compute_implicit_distance(surf, inplace=True)
    threshold = voxel_centers.threshold(test_length)

    distance = voxel_centers['implicit_distance']
    threshold_bool = [(abs(x) < test_length) for x in distance]
    keep_ids = np.where(threshold_bool)[0]

    return vox.extract_cells(keep_ids)


def unstructured_to_poly_convert(grid):
    """
    Convert unstructured grid to polydata
    """
    filterGeom = vtk.vtkGeometryFilter()
    filterGeom.SetInputData(grid)

    filterGeom.Update()

    return pv.wrap(filterGeom.GetOutput())

def getSurfacePoints(vox, surf):
    

    surf_points = []
    for cell_id in range(vox.n_cells):

        point_index = surf.find_closest_point(vox.cell_centers().points[cell_id])
        surf_points.append(point_index)

    return surf_points

def pointPickerRemove(init_points, surf):
    
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


def point_distance(p1,p2):
    return abs(np.linalg.norm(p1 - p2))




def compute_displacement_vec(cp1, cp0):
    """
    computed the displacement between two sets of points
    input: dataframes
    outpt: vector (numpy), and magnitude
    """

    vec = np.vstack([
        cp1.X - cp0.X,
        cp1.Y - cp0.Y,
        cp1.Z - cp0.Z,
        ]).T
    mag = np.linalg.norm(vec, axis = 1)
    return vec, mag

def computeDisplacementAngle(pc1, pc0):
    """
    computed the angle between two displaceemnt vectors
    """

    vdot = np.array([np.dot(a,b) for (a,b) in zip(pc1['disp_vec'], pc0['disp_vec'])])

    angle = np.arccos( vdot / (pc1['mag'] * pc0['mag']))
    

    return angle

def loadControlPoints(folder):
    """
    Load control points of a provided test folder
    returns: list of control point dataframes (xyz)
                list of control point polydata structures with displacement vectors
    """
    
    cp_files = sorted(glob.glob(folder + "/*.fcsv"))
    print(cp_files)
    if len(cp_files) == 0:
        raise ValueError("No control point files found at ", folder)
    
    control_points = []
    point_clouds = []
    
    for cp_file in cp_files:
        
        cp_df = pd.read_csv(cp_file, skiprows=3, usecols=[1,2,3], names=['X', 'Y', 'Z'])
    
        print(len(cp_df))
        control_points.append(cp_df)

        point_cloud = pv.PolyData(cp_df.values)
        point_clouds.append(point_cloud)
        
    for count in range(len(control_points)-1):


        vec, mag = compute_displacement_vec(control_points[count+1], control_points[count])

        control_points[count]['dx'] = vec[:,0]
        control_points[count]['dy'] = vec[:,1]
        control_points[count]['dz'] = vec[:,2]
        control_points[count]['mag'] = mag

        point_clouds[count]['dx'] = vec[:,0]
        point_clouds[count]['dy'] = vec[:,1]
        point_clouds[count]['dz'] = vec[:,2]
        point_clouds[count]['mag'] = mag
        point_clouds[count]['disp_vec'] = vec
        point_clouds[count].set_active_vectors('disp_vec')

    return control_points, point_clouds



def loadSurfaces(folder):
    """
    Returns a list of surfaces
    format: polydata
    
    """

    surf_files = sorted(glob.glob(folder + "/*.ply"))
    print(surf_files)
    if len(surf_files) == 0:
        raise ValueError("No control point files found at ", folder)

    surfaces = []
    for file in surf_files:
        surf = pv.read(file)
        surfaces.append(surf)
        
    return surfaces
        
        
def createPointCloudUniform(surface, box_density=42, show_init_plot=True):
    """
    create point cloud from voxelized mesh. 
    """
    
    #create two voxelized grids - using surf0, bc cp_surf is not closed
    vox_coarse = pv.voxelize(surface, density = surface.length / box_density)
    
    vox_coarse = intersecting_cells(vox_coarse, surface)
    
    coarse_points = getSurfacePoints(vox_coarse, surface)

    if show_init_plot:
        p = pv.Plotter()
        p.add_mesh(surface, opacity = 0.5, color='white')
        p.add_mesh(surface.points[coarse_points], render_points_as_spheres=True, color='blue',point_size=15)
        p.add_mesh(vox_coarse, render_points_as_spheres=True, color='red',point_size=15)
        p.add_text("total points: {:d}".format(len(coarse_points)))
        p.show()

    coarse_pc = pv.PolyData(surface.points[coarse_points])
    remove_points, remove_ids = pointPickerRemove(coarse_pc, surface)

    # add_points, add_ids = pointPickerAdd(coarse_pc, surface)
    # df_add = pd.DataFrame(data = add_points, columns=["l", "p", "s"])
    # print(df_add)


    cp_df = pd.DataFrame(data = coarse_pc.points, columns=["l", "p", "s"])
    cp_df.drop(remove_ids, inplace=True)
    # cp_df_out = cp_df.append(df_add, ignore_index=True)
    
    return cp_df, remove_ids

def createPointCloudRefinement(surface, density_coarse=3, density_fine=1.5, show_init_plot=True):
    """
    create point cloud from voxelized mesh. Refine control points on surfaces with higher curvature
    """
    shift_x = 0
    shift_y = 45
    shift_z = 15
    center = [surface.center[0] + shift_x,
              surface.center[1] + shift_y,
              surface.center[2] + shift_z]
    sphere = pv.Sphere(radius = 44, center = center)
    
    surface['curvature'] = surface.curvature(curv_type='Mean')
    
    high_curvature = surface.threshold(scalars='curvature', value = 0.4)
    low_curvature = surface.threshold(scalars='curvature', value = 0.4, invert=True)
    
    #create two voxelized grids - using surf0, bc cp_surf is not closed
    vox_coarse = pv.voxelize(surface, density = density_coarse)
    vox_fine = pv.voxelize(surface, density = density_fine)
    
    vox_coarse = intersecting_cells(vox_coarse, unstructured_to_poly_convert(low_curvature))
    vox_fine = intersecting_cells(vox_fine, unstructured_to_poly_convert(high_curvature))
    
    coarse_points = getSurfacePoints(vox_coarse, surface)
    fine_points = getSurfacePoints(vox_fine, surface)
    
    points_in, points_out = pointsWithinSphere(pv.PolyData(surface.points[fine_points]), sphere, surface)
    fine_points_in = np.array(fine_points)[points_in]
    
    print(len(fine_points_in))
    print(len(points_out))
    print(len(coarse_points))

    if show_init_plot:
        p = pv.Plotter()
        p.add_mesh(surface, opacity = 0.5, color='white')
        p.add_mesh(surface.points[coarse_points], render_points_as_spheres=True, color='blue',point_size=15)
        p.add_mesh(surface.points[fine_points_in], render_points_as_spheres=True, color='red',point_size=15)
        p.add_text("total points: {:d}".format(len(coarse_points) + len(fine_points_in)))
        p.show()
        
    all_points = np.append(coarse_points, fine_points_in)

    pc = pv.PolyData(surface.points[all_points])
    remove_points, remove_ids = pointPickerRemove(pc, surface)

    add_points, add_ids = pointPickerAdd(pc, surface)
    df_add = pd.DataFrame(data = add_points, columns=["l", "p", "s"])


    cp_df = pd.DataFrame(data = pc.points, columns=["l", "p", "s"])
    cp_df.drop(remove_ids, inplace=True)
    cp_df_out = cp_df.append(df_add, ignore_index=True)
    
    return cp_df_out



        
def savePointsFCSV(file_path, df):
    """
    save a dataframe of fiducials in the correct format for slicer3D
    """

    with open(file_path, 'w') as f:
         f.write("# Markups fiducial file version = 4.11\n")
         f.write("# CoordinateSystem = LPS\n")
         f.write('# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n')

    df.to_csv(file_path, header=False, mode="a")
    
    
    
def badTracksPicker(point_clouds, surfaces):
    """
    Identify bad tracks by selecting the last point of the track on the surface. 
    The last point is represented by a blue dot.
    
    """

    
    picked_points = []
    picked_ids = []

    p = pv.Plotter(notebook=0)
    p.add_mesh(surfaces[0], pickable=False,opacity = 0.4, color='white')
    p.add_mesh(surfaces[-1], pickable=False,opacity = 0.4, color='white')

    for count in range(len(point_clouds)-1):

        point_clouds[count].set_active_vectors('disp_vec')
        p.add_mesh(point_clouds[count].arrows, color=colors[count,0:3], pickable=False)

    p.add_mesh(point_clouds[-1], render_points_as_spheres=True, color='blue', point_size=15, pickable=True)
    p.add_text("\n Select Tracks to Remove")
    
    def point_callback(mesh, picked_id):

        point = mesh.points[picked_id]
        picked_points.append(point)
        picked_ids.append(picked_id)
        p.add_point_labels([point,], ['You picked me!'], point_color='red', point_size=15)

    p.enable_point_picking(point_callback, use_mesh=True)
    p.show()
    
    return picked_points, picked_ids
    
    
def addNewTrack(point_clouds, surfaces):
    
    """
    Add a new track. Steps through all available surfaces and the user can select a single point for the new track on each surface
    Returns the point locations of the new track (numpy ndarray)
    """

    picked_points = []
    picked_ids = []

    def point_callback(mesh, picked_id):

        point = mesh.points[picked_id]
        picked_points.append(point)
        picked_ids.append(picked_id)
        p.add_point_labels([point, ], ['You picked me!'],
                           point_color='red', point_size=15)

    n_surfaces = len(surfaces)

    for count, surface in enumerate(surfaces):
        print(count)
        p = pv.Plotter()
        p.add_mesh(surfaces[0], opacity=0.2, color='white', pickable=False)
        p.add_mesh(surface, pickable=True, opacity=1.0, color='white')

        for inner_count in range(len(point_clouds)-1):

            point_clouds[inner_count].set_active_vectors('disp_vec')
            p.add_mesh(point_clouds[inner_count].arrows,
                       color='red', pickable=False)

        if count == 0:

            point_clouds[count].set_active_vectors('disp_vec')
            p.add_mesh(point_clouds[count].arrows, color=colors[count, 0:3], pickable=False)
            p.add_mesh(point_clouds[count], render_points_as_spheres=True, point_size=15,
                       color='red', pickable=False)
            p.add_text("\n Add first point of new control point track")

        elif count == n_surfaces - 1:
            p.add_mesh(point_clouds[count-1].arrows, color=colors[count,
                       0:3], pickable=False)
            p.add_mesh(pv.PolyData(picked_points), render_points_as_spheres=True,
                       color='blue', point_size=25, pickable=False)
            p.add_mesh(point_clouds[count], render_points_as_spheres=True, point_size=15,
                       color='red', pickable=False)
            p.add_text("\n Add the last point of new control point track")

        else:
            point_clouds[count].set_active_vectors('disp_vec')
            p.add_mesh(point_clouds[count].arrows, color=colors[count, 0:3], pickable=False)
            p.add_mesh(pv.PolyData(picked_points), render_points_as_spheres=True,
                       color='blue', point_size=25, pickable=False)
            p.add_mesh(point_clouds[count], render_points_as_spheres=True, point_size=15,
                       color='red', pickable=False)
            p.add_text("\n Add the next point of new track")

        p.enable_point_picking(point_callback, use_mesh=True)
        p.show()

        print(picked_points)
        
def create_newTrack(point_clouds, surfaces):
    
    """
    Add a new track. Steps through all available surfaces and the user can select a single point for the new track on each surface
    Returns the point locations of the new track (numpy ndarray)
    """

    picked_points = []
    picked_ids = []

    def point_callback(mesh, picked_id):

        point = mesh.points[picked_id]
        picked_points.append(point)
        picked_ids.append(picked_id)
        p.add_point_labels([point, ], ['You picked me!'],
                           point_color='red', point_size=15)

    n_surfaces = len(surfaces)
    
    
    for count, surface in enumerate(surfaces):
        
        p = pv.Plotter(shape=(1,2), window_size=(1500,1200))
        p.subplot(0,1)
        p.add_mesh(surfaces[0], opacity=1.0, color='white', pickable=False)
        
        for inner_count in range(len(point_clouds)-1):
            
            point_clouds[inner_count].set_active_vectors('disp_vec')
            p.add_mesh(surfaces[inner_count], pickable=True, opacity=0.2, color='white')

            p.add_mesh(point_clouds[inner_count].arrows, color=colors[inner_count, 0:3], pickable=False)
            
            if count > 0:
                p.add_mesh(pv.PolyData(picked_points), render_points_as_spheres=True,
                           color='blue', point_size=25, pickable=False)
        
        p.subplot(0,0)
        p.add_mesh(surface, pickable=True, opacity=1.0, color='white')

        for inner_count in range(len(point_clouds)-1):

            point_clouds[inner_count].set_active_vectors('disp_vec')
            p.add_mesh(point_clouds[inner_count].arrows,
                       color='red', pickable=False)

        if count == 0:

            point_clouds[count].set_active_vectors('disp_vec')
            p.add_mesh(point_clouds[count].arrows, color=colors[count, 0:3], pickable=False)
            p.add_mesh(point_clouds[count], render_points_as_spheres=True, point_size=15,
                       color='red', pickable=False)
            p.add_text("\n Add first point of new control point track")

        elif count == n_surfaces - 1:
            p.add_mesh(point_clouds[count-1].arrows, color=colors[count,
                       0:3], pickable=False)
            p.add_mesh(pv.PolyData(picked_points), render_points_as_spheres=True,
                       color='blue', point_size=25, pickable=False)
            p.add_mesh(point_clouds[count], render_points_as_spheres=True, point_size=15,
                       color='red', pickable=False)
            p.add_text("\n Add the last point of new control point track")

        else:
            point_clouds[count].set_active_vectors('disp_vec')
            p.add_mesh(point_clouds[count].arrows, color=colors[count, 0:3], pickable=False)
            p.add_mesh(pv.PolyData(picked_points), render_points_as_spheres=True,
                       color='blue', point_size=25, pickable=False)
            p.add_mesh(point_clouds[count], render_points_as_spheres=True, point_size=15,
                       color='red', pickable=False)
            p.add_text("\n Add the next point of new track")

        p.link_views()
        p.enable_point_picking(point_callback, use_mesh=True, pickable_window=False)
        p.show()

        print(picked_points)
    
    # Visualize the new track
    p = pv.Plotter()
    for count in range(len(point_clouds)-1):
        point_clouds[count].set_active_vectors('disp_vec')
        p.add_mesh(point_clouds[count].arrows, color=colors[count,0:3] )
        p.add_mesh(surfaces[count], color='white', opacity=0.5)
    p.add_mesh(pv.PolyData(picked_points), render_points_as_spheres=True,
               color='red', point_size=25, pickable=False)
    p.show()

    return picked_points, picked_ids


def addNewTrackSimple(point_clouds, surfaces):

    picked_points = []
    picked_ids = []

    def point_callback(mesh, picked_id):

        point = mesh.points[picked_id]
        picked_points.append(point)
        picked_ids.append(picked_id)
        p.add_point_labels([point, ], ['You picked me!'],
                           point_color='red', point_size=15)

    for count, surface in enumerate(surfaces):
        p = pv.Plotter()

        p.add_mesh(surface, pickable=True, opacity=0.6, color='white')

        p.enable_point_picking(point_callback, use_mesh=True)
        p.show()

    return picked_points, picked_ids


def addTrack(new_points, point_clouds):
    """
    Add a new track to the list of point clouds (single track)
    """

    if len(new_points) != len(point_clouds):
        raise ValueError(
            "Number of new points and point clouds to do not match")

    point_clouds_new = []
    control_points_df = []
    for (point, point_cloud) in zip(new_points, point_clouds):

        pc = pv.PolyData(np.vstack((point_cloud.points, point)))
        point_clouds_new.append(pc)
        df = pd.DataFrame(data=pc.points, columns=["X", "Y", "Z"])
        control_points_df.append(df)

    for count in range(len(point_clouds_new)-1):

        vec, mag = compute_displacement_vec(
            control_points_df[count+1], control_points_df[count])

        point_clouds_new[count]['dx'] = vec[:, 0]
        point_clouds_new[count]['dy'] = vec[:, 1]
        point_clouds_new[count]['dz'] = vec[:, 2]

        point_clouds_new[count]['disp_vec'] = vec
        point_clouds_new[count].set_active_vectors('disp_vec')

    return point_clouds_new, control_points_df


def removeTracks(track_ids, point_clouds):
    """
    Delete tracks from point clouds.
    Returns new list of point clouds and dfs with points removed
    
    """

    cleaned_point_clouds = []
    control_points_df = []

    for pc in point_clouds:
        points = pc.points
        cleaned_points = np.delete(points, track_ids, axis=0)
        cleaned_point_clouds.append(pv.PolyData(cleaned_points))

        df = pd.DataFrame(data=cleaned_points, columns=["X", "Y", "Z"])
        control_points_df.append(df)

    for count in range(len(cleaned_point_clouds)-1):

        vec, mag = compute_displacement_vec(
            control_points_df[count+1], control_points_df[count])

        cleaned_point_clouds[count]['dx'] = vec[:, 0]
        cleaned_point_clouds[count]['dy'] = vec[:, 1]
        cleaned_point_clouds[count]['dz'] = vec[:, 2]

        cleaned_point_clouds[count]['disp_vec'] = vec
        cleaned_point_clouds[count].set_active_vectors('disp_vec')

    return cleaned_point_clouds, control_points_df


def starControlPointsCSV(dfs, dt=0.1):
    """
    Format a list of control point dataframes into a single dataframe ready for reading by startccm. 
    """

    df_save = dfs[0]

    for count, df in enumerate(dfs):
        time = count * dt
        column_x = "X[t={:1.3f}s]".format(time)
        column_y = "Y[t={:1.3f}s]".format(time)
        column_z = "Z[t={:1.3f}s]".format(time)

        df_save[column_x] = df.X
        df_save[column_y] = df.Y
        df_save[column_z] = df.Z
        
    return df_save


def pointsWithinSphere(pc, sphere, surface, debug=True):
    """
    Identify points that reside inside and outside a sphere
    """
    
    
    pc_sphere = pc.select_enclosed_points(sphere)
    inside_point_ids = np.where(pc_sphere["SelectedPoints"] == 1)[0]
    outside_point_ids = np.where(pc_sphere["SelectedPoints"] == 0)[0]


    if debug:
        p = pv.Plotter()
        p.add_mesh(sphere, opacity = 0.5)
        p.add_mesh(surface)
        p.add_mesh(pc_sphere, render_points_as_spheres=True, scalars='SelectedPoints')
        p.show()
    
    return inside_point_ids, outside_point_ids

def save_new_point_clouds(point_clouds, export_dir, file_prefix):
    """
    
    Export modified set of point clouds as dfs
    
    """
    
    if not os.path.exists(export_dir):
        os.makedirs(export_dir)
        
    for count,pc in enumerate(point_clouds):
        df = pd.DataFrame(data = pc.points, columns=["X", "Y", "Z"])
        save_name = export_dir + "/" + file_prefix + "_" + str(count) + '.fcsv'
        print(save_name)
        df.to_csv(save_name)


        
    

def pointPicker_Normal(point_clouds, surfaces):
    """
    Select multiple points for integraton with the normal vector function
    
    """    
    picked_points = []
    picked_ids = []

    p = pv.Plotter(shape = (1,2), notebook=0)
    
    p.subplot(0,1)
    for count in range(len(point_clouds) - 1):
        point_clouds[count].set_active_vectors('disp_vec')
        p.add_mesh(point_clouds[count].arrows, color=colors[count,0:3] )
        p.add_mesh(surfaces[count], color='white', opacity=0.5)
    
    
    p.subplot(0,0)
    p.add_mesh(surfaces[0], pickable=True, opacity=1, color='white')
    p.add_mesh(point_clouds[0], render_points_as_spheres=True, color='blue', point_size=15, pickable=False)
    p.add_text("/n Select Points to Add")


    def point_callback(mesh, picked_id):

        point = mesh.points[picked_id]
        picked_points.append(point)
        picked_ids.append(picked_id)
        p.add_point_labels([point,], ['Adding Point!'], point_color='red', point_size=15)

    p.enable_point_picking(point_callback, use_mesh=True)
    p.link_views()
    p.show()
    
    return picked_points, picked_ids



def add_normal_tracks(point_clouds, surfaces, normal_scale=5, show_tracks=True):
    
    """
    Define a new set of tracks that expand based on the surface normal
    
    """

    start_points, starting_point_ids = pointPicker_Normal(point_clouds, surfaces)
    n_tracks = len(start_points)
    n_surfaces = len(surfaces)
    
    point_ids = starting_point_ids

    new_tracks = np.zeros((n_tracks, 3, n_surfaces))
    new_tracks[:,:,0] = start_points
   
    for surf_count in range(len(surfaces) - 1):
        
        surfaces[surf_count].compute_normals(inplace=True)
        intersection_points = []
        
        for point_id in point_ids:
            
            extracted = surfaces[surf_count].extract_points(point_id)
            normal = np.mean(extracted["Normals"], axis=0)
            
            #create line
            p0 = surfaces[surf_count].points[point_id]
            p1 = p0 + normal * normal_scale
            
            intersection_point, intersection_cell = surfaces[surf_count + 1].ray_trace(p0, p1, plot=False)
            
            if len(intersection_point) == 0:
                
                print("Warning: no intersection, try inverting normal")
                p1 = p0 + normal * normal_scale * -1
                intersection_point, intersection_cell = surfaces[surf_count + 1].ray_trace(p0, p1, plot=False)
            
                if len(intersection_point) == 0:
                    raise Exception("No intersection Found. Debug")
                    
            intersection_points.append(intersection_point[0])
            
        #update the new new points ids - find closest in the next surface
        tree = KDTree(surfaces[surf_count + 1].points)
        cdist, point_ids = tree.query(intersection_points, k=1)
        
        new_tracks[:,:,surf_count + 1] = intersection_points
    
    point_clouds_updated = point_clouds.copy()
    for count in range(n_tracks):
        print("Adding track, ", count)
        point_clouds_updated, control_points_updated = addTrack(new_tracks[count, :, :].T, point_clouds_updated)
        
    
    if show_tracks:
        p = pv.Plotter()
        p.add_mesh(pv.PolyData(start_points), render_points_as_spheres=True, color='red', point_size=15)
        for surf_count in range(len(surfaces) - 1):        
            p.add_mesh(surfaces[surf_count], opacity = 0.4)
            p.add_mesh(new_tracks[:,:,surf_count], render_points_as_spheres=True, color='blue', point_size=15)
        p.show()
            
    
    return point_clouds_updated, control_points_updated
        
        
        
        
        
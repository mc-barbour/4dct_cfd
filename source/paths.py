# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:35:25 2024

@author: mbarb1
"""

from sys import platform
from pathlib import Path


if platform == "win32":
    data_dir_sch = Path("Y:\\dahl_j\\4D CT")
    onedrive_dir = Path("D:\Barbour\OneDrive - UW\RobinSequence\Data")
    
elif platform == "darwin":
    data_dir_sch = Path("/Volumes/Active/dahl_j/4D CT")
    onedrive_dir = Path("/Users/mbarb1/OneDrive - UW/RobinSequence/Data/")
    print("need to define mac paths")

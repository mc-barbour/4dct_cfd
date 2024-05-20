# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 12:02:48 2024

@author: mbarb1
"""
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def breathing_waveform_sinusoid(tidal_volume, period, period_insp, D_inlet, D_min):
    
    dt = 0.1
    nu = 1.5e-5
    
    period_exp = period - period_insp
    print(period_exp)
    time_insp = np.linspace(0,period_insp,100)
    time_exp = np.linspace(0,period_exp,100)

    # Inspiration phase
    Q_insp = -np.sin((np.pi/period_insp)*time_insp)
    scale_insp = tidal_volume / np.trapz(Q_insp,x=time_insp)
    Q_insp = scale_insp * np.sin((np.pi/period_insp)*time_insp)
    print("Inspiration tidal volume: ", np.trapz(Q_insp, x=time_insp))
    
    # Expiration phase
    Q_exp = np.sin((np.pi/period_exp)*time_exp)
    scale_exp = tidal_volume / np.trapz(Q_exp, x=time_exp)
    Q_exp = scale_exp*np.sin((np.pi/period_exp)*time_exp)
    print("Expiration tidal volume: ", np.trapz(Q_exp, x=time_exp))

    # combine
    Q = np.append(Q_insp,Q_exp)
    time = np.append(time_insp, time_exp + time_insp[-1])
    
    fig = make_subplots(rows=2, cols=1)
    fig.add_trace(go.Scatter(x=time,y=Q), row=1, col=1)
    
    fig.update_layout(template='plotly_white')
    fig.update_yaxes(title='Flow rate (mL/s)', row=1, col=1)
    
    fig.add_trace(go.Scatter(x=time, y= Q*1e-6*4/(np.pi * D_inlet * nu), name='Inlet Re'), row=2, col=1)
    fig.add_trace(go.Scatter(x=time, y= Q*1e-6*4/(np.pi * D_min * nu), name='constriction Re'), row=2, col=1)

    fig.update_layout(template='plotly_white')
    fig.update_yaxes(title='Flow rate (mL/s)', row=1, col=1)
    fig.update_yaxes(title='Reynolds Number', row=2, col=1)

    fig.update_xaxes(title='time (s)')
    fig.show()
    
    print("Peak Inspiratory Flow: {:.4f} ml/s".format(min(Q)))
    print("Peak Inspiratory Re inlet: {:.4f} ml/s".format(min(Q*1e-6*4/(np.pi * D_inlet * nu))))
    print("Peak Inspiratory Re narrowing: {:.4f} ml/s".format(min(Q*1e-6*4/(np.pi * D_min * nu))))
    
    print("Peak Expiratory Flow: {:.4f} ml/s".format(max(Q)))
    print("Peak Expiratory Re inlet: {:.4f} ml/s".format(max(Q*1e-6*4/(np.pi * D_inlet * nu))))
    print("Peak Expiratory Re narrowing: {:.4f} ml/s".format(max(Q*1e-6*4/(np.pi * D_min * nu))))
    
    return Q, time
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 11:39:33 2020

@author: 1
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes
import time
from datetime import datetime
import json
from mpl_toolkits.basemap import Basemap


def ploth(lon, lat, h_bg, h_an, h_tr, va_value, runtime, curr_time, colorrange, level_in, Numlevel):

    timep1=time.time()

    runtime = np.array(runtime)
    va_value = np.array(va_value)

    fig = plt.figure(figsize=(20, 15))

    basemap = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,resolution='c')

    ax1 = plt.axes([0, 0.55, 0.3, 0.3])
    pcl = plt.imshow(h_bg, cmap='RdBu', vmin=colorrange[0,0], vmax=colorrange[0,1], extent=[min(lon), max(lon), min(lat), max(lat)], interpolation='nearest', origin='lower', aspect='auto')
    cbar = plt.colorbar(pcl, orientation='horizontal', pad=0.1, aspect=20)
#    cbar.set_ticks(range(75000, 92000, 5000))
#    cbar.set_ticklabels(range(75000, 92000, 5000))
    cbar.ax.tick_params(labelsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    ax1.xaxis.set_major_locator(MultipleLocator(60))
    ax1.yaxis.set_major_locator(MultipleLocator(30))
    #ax1.set_yticklabels('') #设定y轴标签文字
    plt.title('h-$\mathregular{X_b}$', fontsize=20)
    #timep3=time.time()
    #print("draw:h_bg_other", timep3-timep2)

    ax2 = plt.axes([0.35, 0.55, 0.3, 0.3])
#    pcl = plt.pcolormesh(lon, lat, h_an, vmin=colorrange[0,0], vmax=colorrange[0,1])

    pcl = plt.imshow(h_an, cmap='RdBu', vmin=colorrange[0,0], vmax=colorrange[0,1], extent=[min(lon), max(lon), min(lat), max(lat)], interpolation='nearest', origin='lower', aspect='auto')
    cbar = plt.colorbar(pcl, orientation='horizontal', pad=0.1, aspect=20)
#    cbar.set_ticks(range(75000, 92000, 5000))
#    cbar.set_ticklabels(range(75000, 92000, 5000))
    cbar.ax.tick_params(labelsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    ax2.xaxis.set_major_locator(MultipleLocator(60))
    ax2.yaxis.set_major_locator(MultipleLocator(30))
    #ax2.set_yticklabels('') #设定y轴标签文字
    plt.title('h-$\mathregular{X_a}$', fontsize=20)

    ax3 = plt.axes([0.7, 0.55, 0.3, 0.3])
#    pcl = plt.pcolormesh(lon, lat, h_tr, vmin=colorrange[0,0], vmax=colorrange[0,1])
    pcl = plt.imshow(h_tr, cmap='RdBu', vmin=colorrange[0,0], vmax=colorrange[0,1], extent=[min(lon), max(lon), min(lat), max(lat)], interpolation='nearest', origin='lower', aspect='auto')
    cbar = plt.colorbar(pcl, orientation='horizontal', pad=0.1, aspect=20)
#    cbar.set_ticks(range(75000, 92000, 5000))
#    cbar.set_ticklabels(range(75000, 92000, 5000))
    cbar.ax.tick_params(labelsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    ax3.xaxis.set_major_locator(MultipleLocator(60))
    ax3.yaxis.set_major_locator(MultipleLocator(30))
    #ax3.set_yticklabels('')
    plt.title('h-$\mathregular{X_t}$', fontsize=20)

    ax4 = plt.axes([0, 0.15, 0.3, 0.3])
#    pcl = plt.pcolormesh(lon, lat, h_an-h_bg, cmap=plt.cm.bwr, vmin=colorrange[3,0], vmax=colorrange[3,1])
    pcl = plt.imshow(h_an-h_bg, cmap='RdBu', vmin=colorrange[1,0], vmax=colorrange[1,1], extent=[min(lon), max(lon), min(lat), max(lat)], interpolation='nearest', origin='lower', aspect='auto')
    cbar = plt.colorbar(pcl, orientation='horizontal', pad=0.1, aspect=20)
#    cbar.set_ticks(range(-100,101,20))
#    cbar.set_ticklabels(range(-100,101,20))
    cbar.ax.tick_params(labelsize=20)
    plt.tick_params(axis='both',which='major',labelsize=20)
    ax4.xaxis.set_major_locator(MultipleLocator(60))
    ax4.yaxis.set_major_locator(MultipleLocator(30))
    plt.title('h-($\mathregular{X_a}$-$\mathregular{X_b}$)',fontsize=20)

    ax5 = HostAxes(fig, [0.38, 0.20, 0.6, 0.20])

#    ax_MeanErr = ParasiteAxes(ax5, sharex=ax5) #寄生的y轴共用ax5的x轴
#    ax_VarRlt = ParasiteAxes(ax5, sharex=ax5) 

#    ax5.parasites.append(ax_MeanErr) #ax_MeanErr拥有一条y轴
#    ax5.parasites.append(ax_VarRlt) #ax_VarRlt拥有一条y轴

    ax5.axis['right'].set_visible(False) #主图右侧的轴看不到
    ax5.axis['top'].set_visible(False) #主图上侧的轴看不到

    fig.add_axes(ax5)

    ax5.set_xlabel('nt', fontsize = 20)#x轴标签名称为nt
    

    timep_rmse0=time.time()
    
    ax5.plot(runtime, va_value[:,0], 'ko-', color='red', label="RMSE-6090N", markersize=5)
    ax5.plot(runtime, va_value[:,1], 'ko-', color='orange', label="RMSE-3060N", markersize=5)
    ax5.plot(runtime, va_value[:,2], 'ko-', color='yellow', label="RMSE-0030N", markersize=5)
    ax5.plot(runtime, va_value[:,3], 'ko-', color='green', label="RMSE-0030S", markersize=5)
    ax5.plot(runtime, va_value[:,4], 'ko-', color='blue', label="RMSE-3060S", markersize=5)
    ax5.plot(runtime, va_value[:,5], 'ko-', color='black', label="RMSE-6090S", markersize=5)
    ax5.legend()
 

    timep_rmse1=time.time()
    print("draw:rmse", timep_rmse1-timep_rmse0)
    
    path = './pics-letkf1.0/h/'+str(level_in)+'/'
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(path+'h_'+curr_time+".png", bbox_inches = 'tight') #bbox_inches尝试找到封装图窗口中所有内容的最紧密的边界>框
    plt.close()


    timep_save=time.time()
    print("save:", timep_save-timep_rmse1)



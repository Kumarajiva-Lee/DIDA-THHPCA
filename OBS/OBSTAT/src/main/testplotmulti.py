# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 11:39:33 2020

@author: 1
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes
import time
from datetime import datetime
from mpl_toolkits.basemap import Basemap
import configparser
from matplotlib.dates import (YEARLY, DateFormatter, rrulewrapper, RRuleLocator, drange)

#import matplotlib
#matplotlib.use('Agg')

#picpath='/home/users/yuz/job/vv/integrated/pic'

def plotlev(bg, an, tr, va_theta, pic_param):

    input_path=pic_param[8]

    cf = configparser.ConfigParser()
    cf.read(input_path)
    m_g=eval(cf.get("obstat_ini", "mpas_or_gmcore"))

    cc = configparser.ConfigParser()
    if m_g == 0:
        cc.read("cc_ini_4mpas")
    else:
        cc.read("cc_ini_4gmcore")


    tp0=time.time()

    runtime = np.array(pic_param[0])
    curr_time = pic_param[1]
    level_in_key = pic_param[2]
    level_in = pic_param[3]
    variable = pic_param[4]
    picpath = pic_param[5]
    cc_case = pic_param[6]
    casename = pic_param[7]
    #va = np.array(va)
    va_theta = np.array(va_theta)

    cc_b_min = float(cc.get(cc_case, cc_case+'_'+variable+'_'+str(level_in_key)+'_min'))
    cc_b_max = float(cc.get(cc_case, cc_case+'_'+variable+'_'+str(level_in_key)+'_max'))
    cc_amb_min = float(cc.get(cc_case, cc_case+'_'+variable+'_'+str(level_in_key)+'_amb_min'))
    cc_amb_max = float(cc.get(cc_case, cc_case+'_'+variable+'_'+str(level_in_key)+'_amb_max'))

    tp1=time.time()
    print('get_param', tp1-tp0)


    fig = plt.figure(figsize=(20, 15))
    ax1 = plt.axes([0.025, 0.55, 0.3, 0.3])
    basemap = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,resolution='c')
    pcl = basemap.imshow(bg, cmap='rainbow', vmin=cc_b_min, vmax=cc_b_max, extent=[0, 360, -90, 90], interpolation='nearest', origin='lower', aspect='auto')
    basemap.drawcoastlines()
    basemap.drawparallels(np.arange(-90,91,30),labels=[1,0,0,0])
    basemap.drawmeridians(np.arange(-180,181,60),labels=[0,0,0,1])

    cbar = plt.colorbar(pcl, orientation='horizontal', pad=0.1, aspect=20) #色标轴，水平分布，间隔0.1，20个颜色
    cbar.ax.tick_params(labelsize=20) #色标大小
    plt.title(variable+'-$\mathregular{X_b}$', fontsize=20) #图题目，大小

    tp2=time.time()
    print('imshow1', tp2-tp1)

#####################################################################################################


    ax2 = plt.axes([0.35, 0.55, 0.3, 0.3])
    pcl = basemap.imshow(an, cmap='rainbow', vmin=cc_b_min, vmax=cc_b_max, extent=[0, 360, -90, 90], interpolation='nearest', origin='lower', aspect='auto')
    basemap.drawcoastlines()
    basemap.drawparallels(np.arange(-90,91,30),labels=[1,0,0,0])
    basemap.drawmeridians(np.arange(-180,181,60),labels=[0,0,0,1])
    cbar = plt.colorbar(pcl, orientation='horizontal', pad=0.1, aspect=20) #色标轴，水平分布，间隔0.1，20个颜色
    cbar.ax.tick_params(labelsize=20) #色标大小
    plt.title(variable+'-$\mathregular{X_a}$', fontsize=20) #图题目，大小

#####################################################################################################

    ax3 = plt.axes([0.675, 0.55, 0.3, 0.3])
    pcl = basemap.imshow(tr, cmap='rainbow', vmin=cc_b_min, vmax=cc_b_max, extent=[0, 360, -90, 90], interpolation='nearest', origin='lower', aspect='auto')
    basemap.drawcoastlines()
    basemap.drawparallels(np.arange(-90,91,30),labels=[1,0,0,0])
    basemap.drawmeridians(np.arange(-180,181,60),labels=[0,0,0,1])
    cbar = plt.colorbar(pcl, orientation='horizontal', pad=0.1, aspect=20) #色标轴，水平分布，间隔0.1，20个颜色
    cbar.ax.tick_params(labelsize=20) #色标大小
    plt.title(variable+'-$\mathregular{X_t}$', fontsize=20) #图题目，大小

######################################################################################################


    ax4 = plt.axes([0.025, 0.15, 0.3, 0.3])
    #pcl = basemap.imshow(an-bg, cmap='RdBu_r', vmin=cc_amb_min, vmax=cc_amb_max, extent=[0, 360, -90, 90], interpolation='nearest', origin='lower', aspect='auto')
    pcl = basemap.imshow(an-bg, cmap='RdBu_r', vmin=-abs((an-bg).max()), vmax=abs((an-bg).max()), extent=[0, 360, -90, 90], interpolation='nearest', origin='lower', aspect='auto')
    basemap.drawcoastlines()
    basemap.drawparallels(np.arange(-90,91,30),labels=[1,0,0,0])
    basemap.drawmeridians(np.arange(-180,181,60),labels=[0,0,0,1])
    cbar = plt.colorbar(pcl, orientation='horizontal', pad=0.1, aspect=20) #色标轴，水平分布，间隔0.1，20个颜色
    cbar.ax.tick_params(labelsize=20) #色标大小
    plt.title(variable+'-$\mathregular{X_a}$-$\mathregular{X_b}$', fontsize=20) #图题目，大小

############################################################################################

    #ax5 = HostAxes(fig, [0.38, 0.35, 0.6, 0.15])
    #ax5.axis['right'].set_visible(False) #主图右侧的轴看不到
    #ax5.axis['top'].set_visible(False) #主图上侧的轴看不到
    #fig.add_axes(ax5)
    #ax5.set_xlabel('nt', fontsize = 20)#x轴标签名称为nt

    
    #ax5.plot(runtime, va[:,0], color='red', label="RMSE-6090N", lw=2, ls=':')
    #ax5.plot(runtime, va[:,1], color='orange', label="RMSE-3060N", lw=2, ls=':')
    #ax5.plot(runtime, va[:,2], color='deepskyblue', label="RMSE-0030N", lw=2, ls=':')
    #ax5.plot(runtime, va[:,3], color='green', label="RMSE-0030S", lw=2, ls=':')
    #ax5.plot(runtime, va[:,4], color='blue', label="RMSE-3060S", lw=2, ls=':')
    #ax5.plot(runtime, va[:,5], color='purple', label="RMSE-6090S", lw=2, ls=':')
    #ax5.plot(runtime, va[:,6], color='black', label="RMSE-level", lw=3, ls='-')
    #ax5.plot(runtime, va[:,7], color='brown', label="RMSE-global", lw=3, ls='-')
    #ax5.set_ylabel(variable+'-RMSE')
    #ax5.legend(loc='upper right')


    ax6 = HostAxes(fig, [0.36, 0.16, 0.6, 0.3])
    ax6.axis['right'].set_visible(False) #主图右侧的轴看不到
    ax6.axis['top'].set_visible(False) #主图上侧的轴看不到
    fig.add_axes(ax6)
    #ax6.set_xlabel('nt', fontsize = 20)#x轴标签名称为nt

    
    ax6.plot(runtime, va_theta[:,0], color='red', label="RMSE-6090N", lw=2, ls=':', marker="o", markersize=5)
    ax6.plot(runtime, va_theta[:,1], color='orange', label="RMSE-3060N", lw=2, ls=':', marker="o", markersize=5)
    ax6.plot(runtime, va_theta[:,2], color='deepskyblue', label="RMSE-0030N", lw=2, ls=':', marker="o", markersize=5)
    ax6.plot(runtime, va_theta[:,3], color='green', label="RMSE-0030S", lw=2, ls=':', marker="o", markersize=5)
    ax6.plot(runtime, va_theta[:,4], color='blue', label="RMSE-3060S", lw=2, ls=':', marker="o", markersize=5)
    ax6.plot(runtime, va_theta[:,5], color='purple', label="RMSE-6090S", lw=2, ls=':', marker="o", markersize=5)
    ax6.plot(runtime, va_theta[:,6], color='black', label="RMSE-level", lw=3, ls='-', marker="o", markersize=5)
    ax6.plot(runtime, va_theta[:,7], color='brown', label="RMSE-global", lw=3, ls='-', marker="o", markersize=5)
    ax6.set_ylim([0, np.max(va_theta)])
    if len(runtime) != 1:
        ax6.set_xlim([runtime[0],runtime[-1]])
    ax6.set_ylabel(variable+'-RMSE-theta')
    ax6.legend(loc='center', bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=15)


    path = picpath+'/pics_level/'+variable+'/'+str(level_in_key)+'/'
    if not os.path.exists(path):
        os.makedirs(path)

    tp3=time.time()

    #plt.savefig(path+variable+'_'+curr_time+".png", bbox_inches = 'tight') #bbox_inches尝试找到封装图窗口中所有内容的最紧密的边界>框
    plt.savefig(path+casename+'_'+variable+'_'+curr_time+".png") #bbox_inches尝试找到封装图窗口中所有内容的最紧密的边界>框

    tp4=time.time()
    print('savefig', tp4-tp3)

    plt.close()





def plots(lev_lat_xx, lev_lat_yy, lev_lat, lev_lon_xx, lev_lon_yy, lev_lon, lat_lon, pic_param):


        curr_time = pic_param[0]
        variable_in_key = pic_param[1]
        level_in = pic_param[2]
        level_in_key = pic_param[3]
        picpath = pic_param[4]
        cc_case = pic_param[5]

        cc_amb_min = float(cc.get(cc_case, cc_case+'_'+variable_in_key+'_'+str(level_in_key)+'_amb_min'))
        cc_amb_max = float(cc.get(cc_case, cc_case+'_'+variable_in_key+'_'+str(level_in_key)+'_amb_max'))

    #colorrange = Colorrange(variable)
    #deltacolorrange = Deltacolorrange(variable)


        #fig, ax = plt.subplots()
        fig = plt.figure(figsize=(20, 15))
        cm = plt.cm.get_cmap('RdBu_r')
        ax1 = plt.axes([0.5, 0.55, 0.5, 0.17])
        sc1 = ax1.scatter(lev_lat_xx,lev_lat_yy,c=lev_lat, vmin=cc_amb_min, vmax=cc_amb_max, cmap=cm)
        plt.tick_params(axis='both', which='both', labelsize=15)
        ax1.xaxis.set_major_locator(MultipleLocator(30))
        ax1.yaxis.set_major_locator(MultipleLocator(5))
        plt.title('vertical-$\mathregular{lev&lat}$', fontsize=15)

        #ax1.set_xlim(-90,90)
        #ax1.set_ylim(0,30)
        #cbar1 = fig.colorbar(sc1)
        ax2 = plt.axes([0.5, 0.30, 0.5, 0.17])
        sc2 = ax2.scatter(lev_lon_xx,lev_lon_yy,c=lev_lon, vmin=cc_amb_min, vmax=cc_amb_max, cmap=cm)
        plt.tick_params(axis='both', which='both', labelsize=15)
        ax2.xaxis.set_major_locator(MultipleLocator(30))
        ax2.yaxis.set_major_locator(MultipleLocator(5))
        plt.title('vertical-$\mathregular{lev&lat}$', fontsize=15)

        cbar = fig.colorbar(sc2,ax=[ax2,ax1], orientation='vertical')
        cbar.ax.tick_params(labelsize=15)

        ax3 = plt.axes([0.1, 0.35, 0.35, 0.35])
        basemap = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,resolution='c')
        pcl = basemap.imshow(lat_lon, cmap='RdBu_r', vmin=cc_amb_min, vmax=cc_amb_max, extent=[0, 360, -90, 90], interpolation='nearest', origin='lower', aspect='auto')
        basemap.drawcoastlines()
        basemap.drawparallels(np.arange(-90,91,30),labels=[1,0,0,0],fontsize=15)
        basemap.drawmeridians(np.arange(-180,181,60),labels=[0,0,0,1],fontsize=15)
        cbar = plt.colorbar(pcl, orientation='horizontal', pad=0.1, aspect=20)
        cbar.ax.tick_params(labelsize=15)
        plt.title('horizontal-$\mathregular{lat&lon}$', fontsize=15)


        path = picpath+'/pics_single/'+variable_in_key+'/'
        if not os.path.exists(path):
               os.makedirs(path)
        plt.savefig(path+variable_in_key+'_'+curr_time+".png", bbox_inches = 'tight') #bbox_inches尝试找到封装图窗口中所有内容的最紧密的边界>框
        plt.close()



def ploti(lev_lat, lev_lon, lon, lat, lat_lon, pic_param):


        curr_time = pic_param[0]
        variable_in_key = pic_param[1]
        level_in = pic_param[2]
        level_in_key = pic_param[3]
        extent_top = pic_param[4]
        picpath = pic_param[5]
        cc_case = pic_param[6]

        cc_amb_min = float(cc.get(cc_case, cc_case+'_'+variable_in_key+'_'+str(level_in_key)+'_amb_min'))
        cc_amb_max = float(cc.get(cc_case, cc_case+'_'+variable_in_key+'_'+str(level_in_key)+'_amb_max'))


        if variable_in_key == 'ps':
            bottom = -0.5
            top = 0.5
        else:
            bottom = 0
            top = extent_top

        #fig, ax = plt.subplo
        fig = plt.figure(figsize=(20, 15))
        cm = plt.cm.get_cmap('RdBu_r')
        ax1 = plt.axes([0.5, 0.55, 0.5, 0.17])
        sc1 = plt.imshow(lev_lat, cmap='RdBu_r', vmin=cc_amb_min, vmax=cc_amb_max, extent=[-90, 90, bottom, top], interpolation='nearest', origin='lower', aspect='auto')
        plt.tick_params(axis='both', which='both', labelsize=15)
        ax1.xaxis.set_major_locator(MultipleLocator(30))
        ax1.yaxis.set_major_locator(MultipleLocator(5))
        plt.title('vertical-$\mathregular{lev&lat}$', fontsize=15)

        #ax1.set_xlim(-90,90)
        #ax1.set_ylim(0,30)
        #cbar1 = fig.colorbar(sc1)
        ax2 = plt.axes([0.5, 0.30, 0.5, 0.17])
        sc2 = plt.imshow(lev_lon, cmap='RdBu_r', vmin=cc_amb_min, vmax=cc_amb_max, extent=[0, 360, bottom, top],interpolation='nearest', origin='lower', aspect='auto')
        plt.tick_params(axis='both', which='both', labelsize=15)
        ax2.xaxis.set_major_locator(MultipleLocator(30))
        ax2.yaxis.set_major_locator(MultipleLocator(5))
        plt.title('vertical-$\mathregular{lev&lon}$', fontsize=15)

        cbar = fig.colorbar(sc2,ax=[ax2,ax1], orientation='vertical')
        cbar.ax.tick_params(labelsize=15)

        ax3 = plt.axes([0.1, 0.35, 0.35, 0.35])
        basemap = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,resolution='c')
        pcl = basemap.imshow(lat_lon, cmap='RdBu_r', vmin=cc_amb_min, vmax=cc_amb_max, extent=[0, 360, -90, 90], interpolation='nearest', origin='lower', aspect='auto')
        basemap.drawcoastlines()
        basemap.drawparallels(np.arange(-90,91,30),labels=[1,0,0,0],fontsize=15)
        basemap.drawmeridians(np.arange(-180,181,60),labels=[0,0,0,1],fontsize=15)
        plt.title('horizontal-$\mathregular{lat&lon}$', fontsize=15)
        cbar = plt.colorbar(pcl, orientation='horizontal', pad=0.1, aspect=20)
        cbar.ax.tick_params(labelsize=15)


        path = picpath+'/pics_single/'+variable_in_key+'/'
        if not os.path.exists(path):
               os.makedirs(path)
        plt.savefig(path+variable_in_key+'_'+curr_time+".png", bbox_inches = 'tight') #bbox_inches尝试找到封装图窗口中所有内容的最紧密的边界>框
        plt.close()



def ploth(va_value, runtime, curr_time, variable, obs_open, picpath):

        runtime = np.array(runtime)
        va_value = np.array(va_value)

        fig = plt.figure(figsize=(20, 15))

        if obs_open != 1:

            #ax5 = HostAxes(fig, [0.2, 0.45, 0.6, 0.15])
            ax5 = plt.axes([0.5, 0.30, 0.5, 0.17])
            ax5.axis['right'].set_visible(False) #主图右侧的轴看不到
            ax5.axis['top'].set_visible(False) #主图上侧的轴看不到
            fig.add_axes(ax5)
            #ax5.set_xlabel('nt', fontsize = 20)#x轴标签名称为nt
            #ax5.plot(runtime, va_value[:,0], 'ko-', color='red', label="grmse", markersize=7)
            ax5.plot(runtime, va_value[:,1], color='black', label="grmse_theta", lw=2, ls=':')
            ax5.set_ylim(bottom=0)
            #ax5.set_xlim(left=runtime[0])
            #ax5.legend(loc='center', bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=15)

        else:

            #ax5 = HostAxes(fig, [0.2, 0.45, 0.6, 0.15])
            ax5 = plt.axes([0.2, 0.75, 0.6, 0.15])
            #ax5=plt.subplots()
            ax5.spines['right'].set_visible(False) #主图右侧的轴看不到
            ax5.spines['top'].set_visible(False) #主图上侧的轴看不到
            #fig.add_axes(ax5)
            #ax5.set_xlabel('nt', fontsize = 20)#x轴标签名称为nt
            ax5.plot(runtime, va_value[:,5], color='purple', label="ormse",  lw=2, ls=':')
            ax5.plot(runtime, va_value[:,6], color='deepskyblue', label="brmse",  lw=2, ls=':')
            ax5.plot(runtime, va_value[:,7], color='red', label="armse", lw=2, ls=':')
            ax5.plot(runtime, va_value[:,9], color='orange', label="grmse_theta", lw=2, ls=':')
            ax5.set_ylim(bottom=0)
            ax5.set_xlim([runtime[0],runtime[-1]])
            ax5.legend(loc='upper right', bbox_to_anchor=(1.0, 1.15), ncol=4, fontsize=10)
            #plt.tick_params(axis='both', which='major', labelsize=10)
            #formatter = DateFormatter('%m/%d/%y')
            #ax5.xaxis.set_major_formatter(formatter)
            #ax5.xaxis.set_tick_params(rotation=30, labelsize=10)
            #plt.xticks(rotation=30)
            ax5.tick_params(axis='x', labelrotation=90, labelsize=10)
            interval=len(runtime)//8
            ax5.xaxis.set_major_locator(MultipleLocator(5))
            #for xtick in ax5.get_xticklabels():
                #xtick.set_rotation(90)

            #ax6 = HostAxes(fig, [0.2, 0.20, 0.6, 0.15])
            ax6 = plt.axes([0.2, 0.45, 0.6, 0.15])
            ax6.spines['right'].set_visible(False) #主图右侧的轴看不到
            ax6.spines['top'].set_visible(False) #主图上侧的轴看不到
            #fig.add_axes(ax6)
            #ax6.set_xlabel('nt', fontsize = 20)#x轴标签名称为nt
            ax6.plot(runtime, va_value[:,3], color='green', label="asd_osd", lw=2, ls='-')
            ax6.plot(runtime, va_value[:,4], color='blue', label="asd_bsd", lw=2, ls='-')
            ax6.set_ylim([0, 90])
            ax6.set_xlim([runtime[0],runtime[-1]])
            #ax6.set_xlim(left=runtime[0])
            ax6.tick_params(axis='x', labelrotation=90, labelsize=10)
            ax6.legend(loc='upper right', bbox_to_anchor=(1.0, 1.15), ncol=4, fontsize=10)



        path = picpath+'/pics_global/'+variable+'/'
        if not os.path.exists(path):
            os.makedirs(path)
        plt.savefig(path+variable+'_'+curr_time+".png", bbox_inches = 'tight') #bbox_inches尝试找到封装图窗口中所有内容的最紧密的边界>框
        plt.close()


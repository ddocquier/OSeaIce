#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 11: Bar plots of ice growth/melt terms, computed via fig12-13.py
PROGRAMMER
    D. Docquier
LAST UPDATE
    17/11/2020
'''

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt

# Working directory
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/OSeaIce_Paper/'

# Parameters
save_fig = True

# Load vfxice
filename = dir_output + 'vfxice_siconc0.npy'
vfxice_D000,vfxice_D012,vfxice_D013,vfxice_D014,vfxice_D015,vfxice_D016,vfxice_D017,vfxice_D018,vfxice_D019,vfxice_D020,vfxice_D021,vfxice_D022,vfxice_D023,vfxice_D024,vfxice_D025,vfxice_D027,vfxice_D028,vfxice_D029,vfxice_D030 = np.load(filename)

# Load vfxbog
filename = dir_output + 'vfxbog_siconc0.npy'
vfxbog_D000,vfxbog_D012,vfxbog_D013,vfxbog_D014,vfxbog_D015,vfxbog_D016,vfxbog_D017,vfxbog_D018,vfxbog_D019,vfxbog_D020,vfxbog_D021,vfxbog_D022,vfxbog_D023,vfxbog_D024,vfxbog_D025,vfxbog_D027,vfxbog_D028,vfxbog_D029,vfxbog_D030 = np.load(filename)

# Load vfxopw
filename = dir_output + 'vfxopw_siconc0.npy'
vfxopw_D000,vfxopw_D012,vfxopw_D013,vfxopw_D014,vfxopw_D015,vfxopw_D016,vfxopw_D017,vfxopw_D018,vfxopw_D019,vfxopw_D020,vfxopw_D021,vfxopw_D022,vfxopw_D023,vfxopw_D024,vfxopw_D025,vfxopw_D027,vfxopw_D028,vfxopw_D029,vfxopw_D030 = np.load(filename)

# Load vfxdyn
filename = dir_output + 'vfxdyn_siconc0.npy'    
vfxdyn_D000,vfxdyn_D012,vfxdyn_D013,vfxdyn_D014,vfxdyn_D015,vfxdyn_D016,vfxdyn_D017,vfxdyn_D018,vfxdyn_D019,vfxdyn_D020,vfxdyn_D021,vfxdyn_D022,vfxdyn_D023,vfxdyn_D024,vfxdyn_D025,vfxdyn_D027,vfxdyn_D028,vfxdyn_D029,vfxdyn_D030 = np.load(filename)

# Load vfxsni
filename = dir_output + 'vfxsni_siconc0.npy'
vfxsni_D000,vfxsni_D012,vfxsni_D013,vfxsni_D014,vfxsni_D015,vfxsni_D016,vfxsni_D017,vfxsni_D018,vfxsni_D019,vfxsni_D020,vfxsni_D021,vfxsni_D022,vfxsni_D023,vfxsni_D024,vfxsni_D025,vfxsni_D027,vfxsni_D028,vfxsni_D029,vfxsni_D030 = np.load(filename)

# Load vfxbom
filename = dir_output + 'vfxbom_siconc0.npy'
vfxbom_D000,vfxbom_D012,vfxbom_D013,vfxbom_D014,vfxbom_D015,vfxbom_D016,vfxbom_D017,vfxbom_D018,vfxbom_D019,vfxbom_D020,vfxbom_D021,vfxbom_D022,vfxbom_D023,vfxbom_D024,vfxbom_D025,vfxbom_D027,vfxbom_D028,vfxbom_D029,vfxbom_D030 = np.load(filename)

# Load vfxsum
filename = dir_output + 'vfxsum_siconc0.npy'
vfxsum_D000,vfxsum_D012,vfxsum_D013,vfxsum_D014,vfxsum_D015,vfxsum_D016,vfxsum_D017,vfxsum_D018,vfxsum_D019,vfxsum_D020,vfxsum_D021,vfxsum_D022,vfxsum_D023,vfxsum_D024,vfxsum_D025,vfxsum_D027,vfxsum_D028,vfxsum_D029,vfxsum_D030 = np.load(filename)

print((vfxice_D012-vfxice_D000)*100.)
print((vfxice_D015-vfxice_D000)*100.)
print((vfxice_D018-vfxice_D000)*100.)
print((vfxice_D021-vfxice_D000)*100.)
print((vfxice_D022-vfxice_D000)*100.)
print((vfxice_D023-vfxice_D000)*100.)


# Maps of difference in mean var between the SST restoring experiments and CTRL (except for vfxice: absolute values)
fig,ax=plt.subplots(3,3,figsize=(18,18))
fig.subplots_adjust(left=0.1,bottom=0.05,right=0.95,top=0.95,wspace=0.3,hspace=0.2)
index = np.arange(7)
bar_width = 1
name_xticks = ['','','','','','','']
yrange = np.arange(-120, 121, 30)
if siconc_threshold == 0.:
    if lat_threshold == 80.:
        yrange2 = np.arange(-40, 40.1, 10)
    else:
        yrange = np.arange(-100, 101, 25)
        yrange2 = np.arange(-5, 5.1, 1)
elif siconc_threshold == 15.:
    yrange2 = np.arange(-9, 9.1, 3)
elif siconc_threshold == 30.:
    yrange2 = np.arange(-12, 12.1, 4)

# D000
ax[0,0].set_title('CTRL',fontsize=26)
ax[0,0].set_title('a',loc='left',fontsize=26,fontweight='bold')
ax[0,0].bar(index[0],vfxbog_D000*100.,bar_width,color='b',label='Basal growth')
ax[0,0].bar(index[1],vfxopw_D000*100.,bar_width,color='g',label='Open-water growth')
ax[0,0].bar(index[2],vfxdyn_D000*100.,bar_width,color='orange',label='Dynamic growth')
ax[0,0].bar(index[3],vfxsni_D000*100.,bar_width,color='gray',label='Snow-ice formation')
ax[0,0].bar(index[4],vfxbom_D000*100.,bar_width,color='r',label='Basal melt')
ax[0,0].bar(index[5],vfxsum_D000*100.,bar_width,color='purple',label='Surface melt')
ax[0,0].bar(index[6],vfxice_D000*100.,bar_width,color='k',label='Net growth')
ax[0,0].set_ylabel('Ice growth (cm year$^{-1}$)',fontsize=24)
ax[0,0].set_yticks(yrange)
ax[0,0].set_xticks(index)
ax[0,0].set_xticklabels(name_xticks)
ax[0,0].tick_params(axis='both',labelsize=18)
ax[0,0].axhline(c='k')
ax[0,0].legend(shadow=True,frameon=False,fontsize=24,bbox_to_anchor=(1.2,1))

# Delete axes
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

# D012
ax[1,0].set_title('ATL1+3$^{\circ}$C',fontsize=26)
ax[1,0].set_title('b',loc='left',fontsize=26,fontweight='bold')
ax[1,0].bar(index[0],(vfxbog_D012-vfxbog_D000)*100.,bar_width,color='b')
ax[1,0].bar(index[1],(vfxopw_D012-vfxopw_D000)*100.,bar_width,color='g')
ax[1,0].bar(index[2],(vfxdyn_D012-vfxdyn_D000)*100.,bar_width,color='orange')
ax[1,0].bar(index[3],(vfxsni_D012-vfxsni_D000)*100.,bar_width,color='gray')
ax[1,0].bar(index[4],(vfxbom_D012-vfxbom_D000)*100.,bar_width,color='r')
ax[1,0].bar(index[5],(vfxsum_D012-vfxsum_D000)*100.,bar_width,color='purple')
ax[1,0].bar(index[6],(vfxice_D012-vfxice_D000)*100.,bar_width,color='k')
ax[1,0].set_ylabel('Ice growth PERT-CTRL \n (cm year$^{-1}$)',fontsize=24)
ax[1,0].set_yticks(yrange2)
ax[1,0].set_xticks(index)
ax[1,0].set_xticklabels(name_xticks)
ax[1,0].tick_params(axis='both',labelsize=18)
ax[1,0].axhline(c='k')

# D015
ax[1,1].set_title('ATL2+3$^{\circ}$C',fontsize=26)
ax[1,1].set_title('c',loc='left',fontsize=26,fontweight='bold')
ax[1,1].bar(index[0],(vfxbog_D015-vfxbog_D000)*100.,bar_width,color='b')
ax[1,1].bar(index[1],(vfxopw_D015-vfxopw_D000)*100.,bar_width,color='g')
ax[1,1].bar(index[2],(vfxdyn_D015-vfxdyn_D000)*100.,bar_width,color='orange')
ax[1,1].bar(index[3],(vfxsni_D015-vfxsni_D000)*100.,bar_width,color='gray')
ax[1,1].bar(index[4],(vfxbom_D015-vfxbom_D000)*100.,bar_width,color='r')
ax[1,1].bar(index[5],(vfxsum_D015-vfxsum_D000)*100.,bar_width,color='purple')
ax[1,1].bar(index[6],(vfxice_D015-vfxice_D000)*100.,bar_width,color='k')
ax[1,1].set_yticks(yrange2)
ax[1,1].set_xticks(index)
ax[1,1].set_xticklabels(name_xticks)
ax[1,1].tick_params(axis='both',labelsize=18)
ax[1,1].axhline(c='k')

# D018
ax[1,2].set_title('ATL3+3$^{\circ}$C',fontsize=26)
ax[1,2].set_title('d',loc='left',fontsize=26,fontweight='bold')
ax[1,2].bar(index[0],(vfxbog_D018-vfxbog_D000)*100.,bar_width,color='b')
ax[1,2].bar(index[1],(vfxopw_D018-vfxopw_D000)*100.,bar_width,color='g')
ax[1,2].bar(index[2],(vfxdyn_D018-vfxdyn_D000)*100.,bar_width,color='orange')
ax[1,2].bar(index[3],(vfxsni_D018-vfxsni_D000)*100.,bar_width,color='gray')
ax[1,2].bar(index[4],(vfxbom_D018-vfxbom_D000)*100.,bar_width,color='r')
ax[1,2].bar(index[5],(vfxsum_D018-vfxsum_D000)*100.,bar_width,color='purple')
ax[1,2].bar(index[6],(vfxice_D018-vfxice_D000)*100.,bar_width,color='k')
ax[1,2].set_yticks(yrange2)
ax[1,2].set_xticks(index)
ax[1,2].set_xticklabels(name_xticks)
ax[1,2].tick_params(axis='both',labelsize=18)
ax[1,2].axhline(c='k')

# D021
ax[2,0].set_title('PAC1+3$^{\circ}$C',fontsize=26)
ax[2,0].set_title('e',loc='left',fontsize=26,fontweight='bold')
ax[2,0].bar(index[0],(vfxbog_D021-vfxbog_D000)*100.,bar_width,color='b')
ax[2,0].bar(index[1],(vfxopw_D021-vfxopw_D000)*100.,bar_width,color='g')
ax[2,0].bar(index[2],(vfxdyn_D021-vfxdyn_D000)*100.,bar_width,color='orange')
ax[2,0].bar(index[3],(vfxsni_D021-vfxsni_D000)*100.,bar_width,color='gray')
ax[2,0].bar(index[4],(vfxbom_D021-vfxbom_D000)*100.,bar_width,color='r')
ax[2,0].bar(index[5],(vfxsum_D021-vfxsum_D000)*100.,bar_width,color='purple')
ax[2,0].bar(index[6],(vfxice_D021-vfxice_D000)*100.,bar_width,color='k')
ax[2,0].set_ylabel('Ice growth PERT-CTRL \n (cm year$^{-1}$)',fontsize=24)
ax[2,0].set_yticks(yrange2)
ax[2,0].set_xticks(index)
ax[2,0].set_xticklabels(name_xticks)
ax[2,0].tick_params(axis='both',labelsize=18)
ax[2,0].axhline(c='k')

# D022
ax[2,1].set_title('PAC2+3$^{\circ}$C',fontsize=26)
ax[2,1].set_title('f',loc='left',fontsize=26,fontweight='bold')
ax[2,1].bar(index[0],(vfxbog_D022-vfxbog_D000)*100.,bar_width,color='b')
ax[2,1].bar(index[1],(vfxopw_D022-vfxopw_D000)*100.,bar_width,color='g')
ax[2,1].bar(index[2],(vfxdyn_D022-vfxdyn_D000)*100.,bar_width,color='orange')
ax[2,1].bar(index[3],(vfxsni_D022-vfxsni_D000)*100., bar_width,color='gray')
ax[2,1].bar(index[4],(vfxbom_D022-vfxbom_D000)*100.,bar_width,color='r')
ax[2,1].bar(index[5],(vfxsum_D022-vfxsum_D000)*100.,bar_width,color='purple')
ax[2,1].bar(index[6],(vfxice_D022-vfxice_D000)*100.,bar_width,color='k')
ax[2,1].set_yticks(yrange2)
ax[2,1].set_xticks(index)
ax[2,1].set_xticklabels(name_xticks)
ax[2,1].tick_params(axis='both',labelsize=18)
ax[2,1].axhline(c='k')

# D023
ax[2,2].set_title('PAC3+3$^{\circ}$C',fontsize=26)
ax[2,2].set_title('g',loc='left',fontsize=26,fontweight='bold')
ax[2,2].bar(index[0],(vfxbog_D023-vfxbog_D000)*100.,bar_width,color='b')
ax[2,2].bar(index[1],(vfxopw_D023-vfxopw_D000)*100.,bar_width,color='g')
ax[2,2].bar(index[2],(vfxdyn_D023-vfxdyn_D000)*100.,bar_width,color='orange')
ax[2,2].bar(index[3],(vfxsni_D023-vfxsni_D000)*100.,bar_width,color='gray')
ax[2,2].bar(index[4],(vfxbom_D023-vfxbom_D000)*100.,bar_width,color='r')
ax[2,2].bar(index[5],(vfxsum_D023-vfxsum_D000)*100.,bar_width,color='purple')
ax[2,2].bar(index[6],(vfxice_D023-vfxice_D000)*100.,bar_width,color='k')
ax[2,2].set_yticks(yrange2)
ax[2,2].set_xticks(index)
ax[2,2].set_xticklabels(name_xticks)
ax[2,2].tick_params(axis='both',labelsize=18)
ax[2,2].axhline(c='k')

# Save figure
if save_fig == True:
    fig.savefig(dir_output+'fig11.png')
    fig.savefig(dir_output+'fig11.eps',dpi=300)

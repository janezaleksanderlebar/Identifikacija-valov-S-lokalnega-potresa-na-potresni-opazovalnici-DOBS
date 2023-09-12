# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 09:22:52 2023

@author: LebarJ
"""

from obspy import read
import os
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import scipy.signal
import numpy as np
from scipy import signal
import statistics 
from os import listdir
from os.path import isfile, join


def list_vseh_potresov(x):
    onlyfiles = [f for f in listdir(r"{a}".format(a=x)) if isfile(join(r"{a}".format(a=x), f))]
    return onlyfiles

def list_potresov_in_branje_v_st(x):
    list_potresa=list_vseh_potresov(x)
    st=read(x+"\\"+list_potresa[0])
    st.clear()
    for a in list_potresa:
        st+=read(x+"\\"+a)
    return st

 
potresi_files_2=[]


rootdir = r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC"
for file in os.listdir(rootdir):
    d = os.path.join(rootdir, file)
    potresi_files_2.append(file)

potresi_files=[]
for a in potresi_files_2:
    potresi_files.append(rootdir+"\\"+a)


vsi_potresi=[]
vsi_potresi_imena=[]


for potresi in potresi_files:
    vsi_potresi.append(list_potresov_in_branje_v_st(potresi))
    t=os.path.basename(potresi)



channel_input = "HHE"
input_2 = "SOS_FILTER_PSD"


#------------------------------------------------------------------------------
# PSD + FILTER ZA E KOMPONENTO
    

# DOBS NO FILTER
if input_2 == "SOS_FILTER_PSD":
    povprecje_array=[]
    for i in range(len(vsi_potresi)):
        for a in vsi_potresi[i]:
            if a.stats.channel == channel_input:
                if a.stats.station == 'DOBS':
                    tr=a
                    
                    # DETREND
                    detrended=scipy.signal.detrend(tr.data, axis=-1, type='linear')
                    
                    # PSD
                    y=detrended
                    f,Pxy=scipy.signal.csd(y,y, fs=200, nperseg=tr.stats.npts)
                    
                    # KONVOLUCIJA
                    w=np.ones(5)
                    w=w/sum(w)
                    convoluted=np.convolve(Pxy, w, mode='full')
                    y=convoluted
                    
                    # POVPREČJE
                    y=y.tolist()
                    povprecje_array.append(y)

                                 
if input_2 == "SOS_FILTER_PSD":

    povprecje=np.average(povprecje_array, axis=0)
    
    
    # GRAF
    y=povprecje*10**(-9)
    
    f=f.tolist()
    y=y.tolist()
    while len(y) > len(f):
        y.pop(-1)
    
    while len(f) > len(y):
        f.pop(-1)
    f_E = np.array(f)
    y_E = np.array(y)

    
    

channel_input = "HHN"

# DOBS NO FILTER HHN
if input_2 == "SOS_FILTER_PSD":
    povprecje_array=[]
    for i in range(len(vsi_potresi)):
        for a in vsi_potresi[i]:
            if a.stats.channel == channel_input:
                if a.stats.station == 'DOBS':
                    tr=a
                    
                    # DETREND
                    detrended=scipy.signal.detrend(tr.data, axis=-1, type='linear')
                    
                    # PSD
                    y=detrended
                    f,Pxy=scipy.signal.csd(y,y, fs=200, nperseg=tr.stats.npts)
                    
                    # KONVOLUCIJA
                    w=np.ones(5)
                    w=w/sum(w)
                    convoluted=np.convolve(Pxy, w, mode='full')
                    y=convoluted
                    
                    # POVPREČJE
                    y=y.tolist()
                    povprecje_array.append(y)

                                 
if input_2 == "SOS_FILTER_PSD":

    povprecje=np.average(povprecje_array, axis=0)
    
    
    # GRAF
    y=povprecje*10**(-9)
    
    f=f.tolist()
    y=y.tolist()
    while len(y) > len(f):
        y.pop(-1)
    
    while len(f) > len(y):
        f.pop(-1)
    f_N = np.array(f)
    y_N = np.array(y)
   


channel_input = "HHZ"

# DOBS NO FILTER HHZ
if input_2 == "SOS_FILTER_PSD":
    povprecje_array=[]
    for i in range(len(vsi_potresi)):
        for a in vsi_potresi[i]:
            if a.stats.channel == channel_input:
                if a.stats.station == 'DOBS':
                    tr=a
                    
                    # DETREND
                    detrended=scipy.signal.detrend(tr.data, axis=-1, type='linear')
                    
                    # PSD
                    y=detrended
                    f,Pxy=scipy.signal.csd(y,y, fs=200, nperseg=tr.stats.npts)
                    
                    # KONVOLUCIJA
                    w=np.ones(5)
                    w=w/sum(w)
                    convoluted=np.convolve(Pxy, w, mode='full')
                    y=convoluted
                    
                    # POVPREČJE
                    y=y.tolist()
                    povprecje_array.append(y)

                                 
if input_2 == "SOS_FILTER_PSD":

    povprecje=np.average(povprecje_array, axis=0)
    
    
    # GRAF
    y=povprecje*10**(-9)
    
    f=f.tolist()
    y=y.tolist()
    while len(y) > len(f):
        y.pop(-1)
    
    while len(f) > len(y):
        f.pop(-1)
    f_Z = np.array(f)
    y_Z = np.array(y)
    



    plt.figure(100)
    # plt.subplot(2, 1, 2,sharex=ax1)
    plt.semilogx(f_E, 10*np.log10(y_E*f_E**2*4*np.pi**2),linewidth=1,\
              label='{x}'.format(x="PSD DOBS HHE"))
    
    plt.figure(100)
    plt.semilogx(f_N, 10*np.log10(y_N*f_N**2*4*np.pi**2),linewidth=1,\
              label='{x}'.format(x="PSD DOBS HHN"))
        
    plt.figure(100)
    plt.semilogx(f_Z, 10*np.log10(y_Z*f_Z**2*4*np.pi**2),linewidth=1,\
              label='{x}'.format(x="PSD DOBS HHZ"))

    plt.xlabel("Frekvenca[Hz]",fontsize="25")
    plt.ylabel("Amplituda [Db]",fontsize="25")
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.plot([6.8], [-0.488], 'rx',markersize=20)
    plt.plot([6.837], [0.561], 'rx',markersize=20)
    plt.plot([6.7999], [-7.2], 'rx',markersize=20)
    plt.plot([12.783], [-8.636], 'rx',markersize=20)
    
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=20)    
    plt.legend(fontsize="30",loc="lower left")    
    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.tight_layout()
    plt.savefig('PSD')      

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

    plt.figure(200)
    # plt.subplot(2, 1, 2,sharex=ax1)
    plt.semilogx(f_E, 10*np.log10(y_E*f_E**2*4*np.pi**2),linewidth=1,\
              label='{x}'.format(x="PSD DOBS HHE"))
    
    plt.figure(200)
    plt.semilogx(f_N, 10*np.log10(y_N*f_N**2*4*np.pi**2),linewidth=1,\
              label='{x}'.format(x="PSD DOBS HHN"))
        
    plt.figure(200)
    plt.semilogx(f_Z, 10*np.log10(y_Z*f_Z**2*4*np.pi**2),linewidth=1,\
              label='{x}'.format(x="PSD DOBS HHZ"))

    plt.xlabel("Frekvenca[Hz]",fontsize="25")
    plt.ylabel("Amplituda [Db]",fontsize="25")
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.plot([6.8], [-0.488], 'rx',markersize=20)
    plt.plot([6.837], [0.561], 'rx',markersize=20)
    plt.plot([6.7999], [-7.2], 'rx',markersize=20)
    plt.plot([12.783], [-8.636], 'rx',markersize=20)
    ax = plt.gca()
    ax.set_xlim([5, 18])
    ax.set_ylim([-30, 2])
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=20)    
    plt.legend(fontsize="30",loc="upper right")    
    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.tight_layout()
    plt.savefig('PSD_ZOOMED')   

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

    plt.figure(300)
    ax1=plt.subplot(1, 3, 1)
    plt.semilogx(f_E, 10*np.log10(y_E*f_E**2*4*np.pi**2),linewidth=1,\
               label='{x}'.format(x="PSD DOBS HHE"))
    plt.plot([6.8], [-0.488], 'rx',markersize=20)
    plt.legend(fontsize="20")    
    ax = plt.gca()
    ax.set_xlim([5, 18])
    ax.set_ylim([-20, 3])
    plt.ylabel("Amplituda [Db]",fontsize="25")

    
    
    plt.figure(300)
    plt.subplot(1, 3, 2,sharex=ax1)
    plt.semilogx(f_N, 10*np.log10(y_N*f_N**2*4*np.pi**2),"orange",linewidth=1,\
              label='{x}'.format(x="PSD DOBS HHN"))
    plt.plot([6.837], [0.561], 'rx',markersize=20)
    plt.legend(fontsize="20")    
    ax = plt.gca()
    ax.set_xlim([5, 18])
    ax.set_ylim([-20, 3])
    plt.xlabel("Frekvenca[Hz]",fontsize="25")


        
    plt.figure(300)
    plt.subplot(1, 3, 3,sharex=ax1)
    plt.semilogx(f_Z, 10*np.log10(y_Z*f_Z**2*4*np.pi**2),"green",linewidth=1,\
              label='{x}'.format(x="PSD DOBS HHZ"))
    plt.plot([6.7999], [-7.2], 'rx',markersize=20)
    plt.plot([12.783], [-8.636], 'rx',markersize=20)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax = plt.gca()
    ax.set_xlim([5, 18])
    ax.set_ylim([-30, 3])
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=20)    
    plt.legend(fontsize="20")    
    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.tight_layout()
    plt.savefig('PSD_ZOOMED_POSEBEJ')  
    
    
    
    
    
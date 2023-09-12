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

if input_2 == "SOS_FILTER_PSD":

# 2 VRH
    fc=6.2
    r=0.945

    a_0=1
    a_1=-2*r*np.cos(2*np.pi*fc/200)
    a_2=r**2
    
    denum=[a_0,a_1,a_2]
    num=[0.1]

    sos_1=signal.tf2sos(num, denum, analog=False)

# 3 VRH
    fc=6.9
    r=0.9999

    a_0=12.5
    a_1=-2*r*np.cos(2*np.pi*fc/200)
    a_2=r**2
    
    denum=[a_0,a_1,a_2]
    num=[0.4]

    sos_3=signal.tf2sos(num, denum, analog=False)

# 1 VRH
    fc=8
    r=0.955
    
    a_0=1
    a_1=-2*r*np.cos(2*np.pi*fc/200)
    a_2=r**2

    denum=[a_0,a_1,a_2]
    num=[0.49]
    sos_2=signal.tf2sos(num, denum, analog=False)
  
    
# LOW PASS FILTER
    num, denum = scipy.signal.iirfilter(2, Wn=8.5, fs=200, btype="low", ftype="butter")
    sos_low=signal.tf2sos(denum, num, analog=False)
    
    num, denum = scipy.signal.iirfilter(2, Wn=9, fs=200, btype="low", ftype="butter")
    sos_low_2=signal.tf2sos(denum, num, analog=False)
    
    num, denum = scipy.signal.iirfilter(2, Wn=11, fs=200, btype="low", ftype="cheby1",rp=3,rs=8)
    sos_low_3=signal.tf2sos(denum, num, analog=False)
    
    num, denum = scipy.signal.iirfilter(2, Wn=10, fs=200, btype="low", ftype="bessel")
    sos_low_4=signal.tf2sos(denum, num, analog=False)
    
    sos_skp=np.concatenate((sos_1,sos_2,sos_3,sos_low,sos_low_2,sos_low_3,sos_low_4),axis=0)
    # sos_skp=np.concatenate((sos_1,sos_2,sos_3),axis=0)
    # sos_skp=np.concatenate((sos_1,sos_2,sos_3,sos_band),axis=0)
    # sos_skp=np.concatenate((sos_1,sos_2,sos_3),axis=0)
    
    
    num,denum=signal.sos2tf(sos_skp)
    sos_skp=signal.tf2sos(denum, num, analog=False)


    povprecje_array=[]
    for i in range(len(vsi_potresi)):
        for a in vsi_potresi[i]:
            if a.stats.channel == channel_input:
                if a.stats.station == 'DOBS':
                    tr=a
                    
                    # DETREND
                    y=scipy.signal.detrend(tr.data, axis=-1, type='linear')
                    
                    # FILTER
                    y=scipy.signal.sosfilt(sos_skp,y,axis=-1)
                    # y= scipy.signal.lfilter(denum, num, y, axis=-1)
                    
                    # PSD
                    y=scipy.signal.detrend(y, axis=-1, type='linear')
                    f,y=scipy.signal.csd(y,y, fs=200, nperseg=tr.stats.npts)
                    
                    # KONVOLUCIJA
                    w=np.ones(5)
                    w=w/sum(w)
                    y=np.convolve(y, w, mode='full')
                    
                    # POVPREČJE
                    y=y.tolist()
                    povprecje_array.append(y)
                    
                                 
if input_2 == "SOS_FILTER_PSD":

    povprecje=np.mean(povprecje_array, axis=0)
    
    # GRAF
    y=povprecje*10**(-9)
    
    
    f=f.tolist()
    y=y.tolist()
    while len(y) > len(f):
        y.pop(-1)
    
    while len(f) > len(y):
        f.pop(-1)
    f = np.array(f)
    y = np.array(y)
    
    plt.figure(100)
    plt.semilogx(f, 10*np.log10(y*f**2*4*np.pi**2),linewidth=1.5,\
              label='{x}'.format(x="povprecje DOBS"))
    plt.legend()
    
    # RISANJE FILTRA
    plt.figure(200)
    w,h=scipy.signal.freqz(denum,num, worN=1024, whole=False, plot=None, fs=200, include_nyquist=True)
    plt.xlabel("Frekvenca[Hz]",fontsize="25")
    plt.ylabel("Amplituda [Db]",fontsize="25")
    plt.legend(fontsize="25")
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=20)
    plt.semilogx(w,  20 * np.log10(abs(h)),linewidth=1.5,label='{x}'.format(x="Resonančni filter"))
    plt.legend(fontsize="20")


# OSTALI PODATKI NE DOBS
if input_2 == "SOS_FILTER_PSD":
    povprecje_array=[]
    for i in range(len(vsi_potresi)):
        for a in vsi_potresi[i]:
            if a.stats.channel == channel_input:
                if a.stats.station != 'DOBS':
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
    f = np.array(f)
    y = np.array(y)
    
    plt.figure(100)
    plt.semilogx(f, 10*np.log10(y*f**2*4*np.pi**2),linewidth=0.5,\
              label='{x}'.format(x="povprecje_vseh_ostalih"))
    plt.legend()    
    

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
    f = np.array(f)
    y = np.array(y)
    
    plt.figure(100)
    plt.semilogx(f, 10*np.log10(y*f**2*4*np.pi**2),linewidth=0.5,\
              label='{x}'.format(x="normalen DOBS"))
    plt.legend()    
    

    
    
    
    
    
    
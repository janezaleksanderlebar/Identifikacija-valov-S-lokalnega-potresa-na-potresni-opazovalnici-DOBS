from obspy import read
import os
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import scipy.signal
import numpy as np
from scipy import signal
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


potresi_files=[
    
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data 2.0\data 2.0\S20211221144400_NovoMesto_Mlv3.2",
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data 2.0\data 2.0\S20211225073000_Ljutomer_Mlv2.0",
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data 2.0\data 2.0\S20220115201200_Zagorje_2.6",
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data 2.0\data 2.0\S20220308215300_DolPriLitiji_2.0"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20150210073731"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20150531231013"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20190122180210"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20150529055815"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20170308012354"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220510122540"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220521155147"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20170308012354"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20170706165833"          
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20170711220254" 
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220813194201" #okish
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220526211437"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220521155147"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20221222110630"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220813194201"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220526211437" #okish
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220521155147"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220510122540"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220504113556"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220502191414"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220308215358"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220302023259"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220219061434"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220203170258"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220124220049" #DOBRA ZA POKAZAT BOLJŠI FILTER Z
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20220115201233"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20211124145420"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20210809160826" #ok
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20210223103633"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20210111163148"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20201019063846"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20200828194457"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20200803202335"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20200729141120"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20200621192209"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20200621181215"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20200414052720" #DOBRA Z
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20191013075614"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20190427093132" #mogoče
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20190122180210"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20190427093132"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20190122180210"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20181225182541"
# r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20181116141202"
r"C:\Users\lebarj\Desktop\Janez\Diploma\data SAC\S20181104041446" #mogoče

                ]




vsi_potresi=[]
vsi_potresi_imena=[]


for potresi in potresi_files:
    vsi_potresi.append(list_potresov_in_branje_v_st(potresi))
    t=os.path.basename(potresi)
    

channel_input = "HHZ"
input_2="RISANJE GRAFA FILTRIRANEGA Z SOS FILTROM SAMO ZA Z"


if input_2 == "RISANJE GRAFA FILTRIRANEGA Z SOS FILTROM SAMO ZA Z":

# 2 VRH
    fc=12.6
    r=0.98

    a_0=1
    a_1=-2*r*np.cos(2*np.pi*fc/200)
    a_2=r**2
    
    denum=[a_0,a_1,a_2]
    num=[0.25]

    sos_1=signal.tf2sos(num, denum, analog=False)

# 3 VRH
    fc=12.5
    r=0.99

    a_0=12.5
    a_1=-2*r*np.cos(2*np.pi*fc/200)
    a_2=r**2
    
    denum=[a_0,a_1,a_2]
    num=[0.45]

    sos_3=signal.tf2sos(num, denum, analog=False)

# 1 VRH
    fc=6.8
    r=0.98
    
    a_0=1
    a_1=-2*r*np.cos(2*np.pi*fc/200)
    a_2=r**2

    denum=[a_0,a_1,a_2]
    num=[0.2]
    sos_2=signal.tf2sos(num, denum, analog=False)


# 1 VRH
    fc=15
    r=0.9
    
    a_0=1
    a_1=-2*r*np.cos(2*np.pi*fc/200)
    a_2=r**2

    denum=[a_0,a_1,a_2]
    num=[0.52]
    sos_4=signal.tf2sos(num, denum, analog=False)
    
    
    

    num, denum = scipy.signal.iirfilter(2, Wn=13, fs=200, btype="low", ftype="butter")
    sos_low=signal.tf2sos(denum, num, analog=False)
    
    num, denum = scipy.signal.iirfilter(2, Wn=14, fs=200, btype="low", ftype="butter")
    sos_low_2=signal.tf2sos(denum, num, analog=False)
    
    num, denum = scipy.signal.iirfilter(2, Wn=15, fs=200, btype="low", ftype="cheby1",rp=3,rs=10)
    sos_low_3=signal.tf2sos(denum, num, analog=False)
    
    # num, denum = scipy.signal.iirfilter(2, Wn=12.8, fs=200, btype="low", ftype="bessel")
    # sos_low_4=signal.tf2sos(denum, num, analog=False)
    
    
    sos_skp=np.concatenate((sos_1,sos_2,sos_3,sos_low,sos_low_2,sos_low_3,sos_4),axis=0)
    # sos_skp=np.concatenate((sos_1,sos_2,sos_3,sos_band),axis=0)
    # sos_skp=np.concatenate((sos_1,sos_2,sos_3),axis=0)
    
    
    num,denum=signal.sos2tf(sos_skp)    
    sos_skp=signal.tf2sos(denum, num, analog=False)

# RISANJE SAMO DOBS FILTRIRAN -------------------------------------------------
    lil_count=0
    premik=[]
    st_grafov=0
    if input_2 == "RISANJE GRAFA FILTRIRANEGA Z SOS FILTROM SAMO ZA Z":
        for i in range(len(vsi_potresi)):
            for a in vsi_potresi[i]:
                if a.stats.channel == channel_input:
                    if a.stats.station == "DOBS":
                        premik.append(max(a.data))
                        tr=a
                        y=scipy.signal.detrend(tr.data, axis=-1, type='linear')

                        y=scipy.signal.sosfilt(sos_skp,y,axis=-1)
                        # y = signal.lfilter(denum, num, tr.data)

                        # y=y+50000*st_grafov #AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                        st_grafov+=1
                        
                        y_izb_filt=y*10**(-9)
                        



# RISANJE DOBS ----------------------------------------------------------------     
    if input_2 == "RISANJE GRAFA FILTRIRANEGA Z SOS FILTROM SAMO ZA Z":
        for i in range(len(vsi_potresi)):
            for a in vsi_potresi[i]:
                if a.stats.channel == channel_input:
                    if a.stats.station == "DOBS":
                        tr=a
                        y=tr.data
                        y=scipy.signal.detrend(y, axis=-1, type='linear')

                        # y=y+20000 #AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                        y_no_filt=y*10**(-9)
                        
                        

            
                    
    
            
    
# 3 VRH
    fc=12.5
    r=0.9

    a_0=12.5
    a_1=-2*r*np.cos(2*np.pi*fc/200)
    a_2=r**2
    
    denum=[a_0,a_1,a_2]
    num=[0.7725]

    sos_3=signal.tf2sos(num, denum, analog=False)

# 1 VRH
    fc=6.8
    r=0.9
    
    a_0=1
    a_1=-2*r*np.cos(2*np.pi*fc/200)
    a_2=r**2

    denum=[a_0,a_1,a_2]
    num=[0.7725]
    sos_2=signal.tf2sos(num, denum, analog=False)
    
    
    sos_skp=np.concatenate((sos_3,sos_2),axis=0)
    num,denum=signal.sos2tf(sos_skp)    
    sos_skp=signal.tf2sos(denum, num, analog=False)

    lil_count=0
    premik=[]
    st_grafov=0
    if input_2 == "RISANJE GRAFA FILTRIRANEGA Z SOS FILTROM SAMO ZA Z":
        for i in range(len(vsi_potresi)):
            for a in vsi_potresi[i]:
                if a.stats.channel == channel_input:
                    if a.stats.station == "DOBS":
                        premik.append(max(a.data))
                        tr=a
                        y=scipy.signal.detrend(tr.data, axis=-1, type='linear')

                        y=scipy.signal.sosfilt(sos_skp,y,axis=-1)
                        # y = signal.lfilter(denum, num, tr.data)

                        # y=y+60000 #AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                        st_grafov+=1
                        
                        y_osnvn_filt=y*10**(-9)
                        

                        
        # if input_2 == "RISANJE GRAFA FILTRIRANEGA Z SOS FILTROM SAMO ZA Z":
        #         for i in range(len(vsi_potresi)):
        #             for a in vsi_potresi[i]:
        #                 if a.stats.channel == channel_input:
        #                     if a.stats.station != "DOBS":
        #                         tr=a
                                
        #                         y=tr.data
        #                         y=scipy.signal.detrend(y, axis=-1, type='linear')                        
                                
        #                         y=y-20000*st_grafov #AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        #                         y=y*10**(-9)
                                
        #                         st_grafov+=1

        #                         plt.figure(100)
        #                         plt.semilogx(a.times(), y,linewidth=1,
        #                                   label='{x}'.format(x=tr.stats.station))
                                
        #                         # plt.xlim([4, 10])
        #                         plt.xlabel("Frekvenca[Hz]",fontsize="25")
        #                         plt.ylabel("Amplituda [Db]",fontsize="25")
        #                         plt.legend(fontsize="25")
        #                         plt.xticks(fontsize=20)
        #                         plt.yticks(fontsize=20)
        #                         plt.tick_params(axis='both', which='major', labelsize=20)
        #                         plt.tick_params(axis='both', which='minor', labelsize=20)    
        #                         plt.legend(fontsize="20")  


plt.rc('font', **{'size':'20'})

risanje_test_1= False
risanje_test_2= False
risanje_test_3= False
risanje_test_4= False
risanje_test_5= False 
risanje_test_6= False  
risanje_test_7= False  

       
# Risanje grafov samodejno shranjevanje grafov:
# 1=Brez filtra, Normalen filter, Izboljšan filter
# 2=Brez filtra, Normalen filter
# 3=normalen filter graf
# 4=Izboljšan filter graf
# 5=ZOOMED IN  Brez filtra, Normalen filter, Izboljšan filter
# 6=ZOOMED IN  Brez filtra, Normalen filter
# 7=končna analiza

risanje_test_1= True
risanje_test_2= True
risanje_test_3= True
risanje_test_4= True
risanje_test_5= True
risanje_test_6= True
risanje_test_7=True


# Risanje grafov:

plt.figure(100)
ax1=plt.subplot(3, 1, 3)
plt.plot(a.times(), y_izb_filt,linewidth=1.5,
             label='{x}'.format(x="Filter z več resonančnimi vrhovi"))
plt.axvline(x=12.45, color='r')
plt.xlabel("Čas[s]",fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)    
plt.legend(fontsize="25",loc="lower right")

                        
plt.figure(100)
plt.subplot(3, 1, 1,sharex=ax1)
plt.plot(a.times(), y_no_filt,linewidth=1.5,
             label='{x}'.format(x="Brez filtra"))
plt.axvline(x=12.45, color='r',linestyle = 'dashed')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.legend(fontsize="25",loc="lower right")

plt.figure(100)
plt.subplot(3, 1, 2,sharex=ax1)
plt.plot(a.times(), y_osnvn_filt,linewidth=1.5,
             label='{x}'.format(x="Enojni resonančni filter"))
plt.axvline(x=12.45, color='r')
plt.ylabel("Premik[nm]",fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.legend(fontsize="25",loc="lower right")
manager = plt.get_current_fig_manager()
manager.full_screen_toggle()
plt.tight_layout()

if risanje_test_1 == True:    
    plt.savefig('Z_Normal_Filter_izbfilt.png')


# Risanje grafov brez izboljšanega
plt.figure(200)
plt.subplot(2, 1, 1,sharex=ax1)
plt.plot(a.times(), y_no_filt,linewidth=1.5,
             label='{x}'.format(x="Brez filtra"))
plt.axvline(x=12.45, color='r',linestyle = 'dashed')
plt.ylabel("Premik[nm]",fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.legend(fontsize="30",loc="upper right")

plt.figure(200)
plt.subplot(2, 1, 2,sharex=ax1)
plt.plot(a.times(), y_osnvn_filt,linewidth=1.5,
             label='{x}'.format(x="Enojni resonančni filter"))
plt.xlabel("Čas[s]",fontsize="25")
plt.axvline(x=12.45, color='r')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.legend(fontsize="30")
manager = plt.get_current_fig_manager()
manager.full_screen_toggle()
plt.tight_layout()

if risanje_test_2 == True:    
    plt.savefig('Z_Normal_Filter.png')      


# Risanje ZOOMED IN 
plt.figure(500)
ax1=plt.subplot(3, 1, 3)
plt.plot(a.times(), y_izb_filt,linewidth=1.5,
             label='{x}'.format(x="Filter z več resonančnimi vrhovi"))
plt.axvline(x=12.45, color='r')
plt.xlabel("Čas[s]",fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)    
plt.legend(fontsize="25",loc="lower left")

                        
plt.figure(500)
plt.subplot(3, 1, 1,sharex=ax1)
plt.plot(a.times(), y_no_filt,linewidth=1.5,
             label='{x}'.format(x="Brez filtra"))
plt.axvline(x=12.45, color='r',linestyle = 'dashed')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.legend(fontsize="25",loc="lower left")

plt.figure(500)
plt.subplot(3, 1, 2,sharex=ax1)
plt.plot(a.times(), y_osnvn_filt,linewidth=1.5,
             label='{x}'.format(x="Enojni resonančni filter"))
plt.axvline(x=12.45, color='r')
plt.ylabel("Premik[nm]",fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.legend(fontsize="25",loc="lower left")
ax = plt.gca()
ax.set_xlim([5, 15])
manager = plt.get_current_fig_manager()
manager.full_screen_toggle()
plt.tight_layout()

if risanje_test_5 == True:    
    plt.savefig('Z_ZOOMED_Normal_Filter_izbfilt.png')

# Risanje ZOOMED IN grafov brez izboljšanega
plt.figure(600)
plt.subplot(2, 1, 1,sharex=ax1)
plt.plot(a.times(), y_no_filt,linewidth=1.5,
             label='{x}'.format(x="Brez filtra"))
plt.axvline(x=12.45, color='r',linestyle = 'dashed')
plt.ylabel("Premik[nm]",fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.legend(fontsize="25",loc="lower left")

plt.figure(600)
plt.subplot(2, 1, 2,sharex=ax1)
plt.plot(a.times(), y_osnvn_filt,linewidth=1.5,
             label='{x}'.format(x="Enojni resonančni filter"))
plt.xlabel("Čas[s]",fontsize="25")
plt.axvline(x=12.45, color='r')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.legend(fontsize="25",loc="lower left")
ax = plt.gca()
ax.set_xlim([5, 15])
manager = plt.get_current_fig_manager()
manager.full_screen_toggle()
plt.tight_layout()

if risanje_test_6 == True:    
    plt.savefig('Z_ZOOMED_Normal_Filter.png')      


#Risanje filtra osnovnega 

# 3 VRH
fc=12.5
r=0.9

a_0=12.5
a_1=-2*r*np.cos(2*np.pi*fc/200)
a_2=r**2

denum=[a_0,a_1,a_2]
num=[0.7725]

sos_3=signal.tf2sos(num, denum, analog=False)

# 1 VRH
fc=6.8
r=0.9

a_0=1
a_1=-2*r*np.cos(2*np.pi*fc/200)
a_2=r**2

denum=[a_0,a_1,a_2]
num=[0.7725]
sos_2=signal.tf2sos(num, denum, analog=False)


sos_skp=np.concatenate((sos_3,sos_2),axis=0)
num,denum=signal.sos2tf(sos_skp)
sos_skp=signal.tf2sos(denum, num, analog=False)

plt.figure(300)
ax1=plt.subplot(2, 1, 1)
w,h=scipy.signal.freqz(num,denum, worN=1024, whole=False, plot=None, fs=200, include_nyquist=True)
plt.semilogx(w,  20 * np.log10(abs(h)),linewidth=2.5,
             label='{x}'.format(x="Resonančni filter"))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel("Amplituda [Db]",fontsize="25")
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)    
plt.legend(fontsize="30",loc='lower left') 
ax = plt.gca()
ax.set_xlim([0.1, 25])
ax.set_ylim([-50, 15])


plt.figure(300)
plt.subplot(2, 1, 2,sharex=ax1)
w,h=scipy.signal.freqz(denum,num, worN=1024, whole=False, plot=None, fs=200, include_nyquist=True)
plt.semilogx(w,  20 * np.log10(abs(h)),linewidth=2.5,
             label='{x}'.format(x="Inverzni resonančni filter"))
ax = plt.gca()
ax.set_xlim([0.1, 25])
ax.set_ylim([-15, 50])
plt.xlabel("Frekvenca[Hz]",fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)    
plt.legend(fontsize="30",loc='upper left')  
manager = plt.get_current_fig_manager()
manager.full_screen_toggle()
plt.tight_layout()

if risanje_test_3 == True:    
    plt.savefig('Z_normalfilt')      



# Risanje izboljšanega filtra:
    # 2 VRH
fc=12.6
r=0.98

a_0=1
a_1=-2*r*np.cos(2*np.pi*fc/200)
a_2=r**2

denum=[a_0,a_1,a_2]
num=[0.25]

sos_1=signal.tf2sos(num, denum, analog=False)

    # 3 VRH
fc=12.5
r=0.99

a_0=12.5
a_1=-2*r*np.cos(2*np.pi*fc/200)
a_2=r**2

denum=[a_0,a_1,a_2]
num=[0.45]

sos_3=signal.tf2sos(num, denum, analog=False)

    # 1 VRH
fc=6.8
r=0.98

a_0=1
a_1=-2*r*np.cos(2*np.pi*fc/200)
a_2=r**2

denum=[a_0,a_1,a_2]
num=[0.2]
sos_2=signal.tf2sos(num, denum, analog=False)


    # 1 VRH
fc=15
r=0.9

a_0=1
a_1=-2*r*np.cos(2*np.pi*fc/200)
a_2=r**2

denum=[a_0,a_1,a_2]
num=[0.52]
sos_4=signal.tf2sos(num, denum, analog=False)




num, denum = scipy.signal.iirfilter(2, Wn=13, fs=200, btype="low", ftype="butter")
sos_low=signal.tf2sos(denum, num, analog=False)

num, denum = scipy.signal.iirfilter(2, Wn=14, fs=200, btype="low", ftype="butter")
sos_low_2=signal.tf2sos(denum, num, analog=False)

num, denum = scipy.signal.iirfilter(2, Wn=15, fs=200, btype="low", ftype="cheby1",rp=3,rs=10)
sos_low_3=signal.tf2sos(denum, num, analog=False)

# num, denum = scipy.signal.iirfilter(2, Wn=12.8, fs=200, btype="low", ftype="bessel")
# sos_low_4=signal.tf2sos(denum, num, analog=False)


sos_skp=np.concatenate((sos_1,sos_2,sos_3,sos_low,sos_low_2,sos_low_3,sos_4),axis=0)
# sos_skp=np.concatenate((sos_1,sos_2,sos_3,sos_band),axis=0)
# sos_skp=np.concatenate((sos_1,sos_2,sos_3),axis=0)


num,denum=signal.sos2tf(sos_skp)    
sos_skp=signal.tf2sos(denum, num, analog=False)


plt.figure(400)
ax1=plt.subplot(2, 1, 1)
w,h=scipy.signal.freqz(num,denum, worN=1024, whole=False, plot=None, fs=200, include_nyquist=True)
plt.semilogx(w,  20 * np.log10(abs(h)),linewidth=2.5,
             label='{x}'.format(x="Resonančni filter"))
plt.ylabel("Amplituda [Db]",fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)    
plt.legend(fontsize="30",loc='upper left') 
ax = plt.gca()
ax.set_xlim([0.1, 75])
ax.set_ylim([-50, 85])


plt.figure(400)
plt.subplot(2, 1, 2,sharex=ax1)
w,h=scipy.signal.freqz(denum,num, worN=1024, whole=False, plot=None, fs=200, include_nyquist=True)
plt.semilogx(w,  20 * np.log10(abs(h)),linewidth=2.5,
             label='{x}'.format(x="Inverzni resonančni filter"))
ax = plt.gca()
ax.set_xlim([0.1, 75])
ax.set_ylim([-85, 50])
plt.xlabel("Frekvenca[Hz]",fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)    
plt.legend(fontsize="30",loc='upper left')  
manager = plt.get_current_fig_manager()
manager.full_screen_toggle()
plt.tight_layout()

if risanje_test_4 == True:    
    plt.savefig('Z_izblfilt')      


plt.figure(700)
ax1=plt.subplot(3, 1, 3)
plt.plot(a.times(), y,linewidth=1.5,
             label='{x}'.format(x="Filter z več resonančnimi vrhovi"))
# plt.xlabel("Čas[s]",fontsize="25")
# plt.ylabel("Premik[nm]",fontsize="25")
plt.axvline(x=12.46, color='r')
plt.xlabel("Čas[s]",fontsize="25")
ax = plt.gca()
ax.set_xlim([10, 14])
# ax.set_ylim([-85, 50])
plt.legend(fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.legend(fontsize="25",loc="lower left")


plt.figure(700)
# ax2 = plt.subplot(2,1,2, )
plt.subplot(3, 1, 1,sharex=ax1)
plt.plot(a.times(), y_no_filt,linewidth=1.5,
             label='{x}'.format(x="Brez filtra"))
# plt.xlabel("Čas[s]",fontsize="25")
plt.axvline(x=12.4, color='r',linestyle = 'dashed',linewidth=70,alpha=0.2)
plt.legend(fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)    
plt.legend(fontsize="25",loc="lower left")

 
plt.figure(700)
plt.subplot(3, 1, 2,sharex=ax1)
plt.plot(a.times(),  y_osnvn_filt,linewidth=1.5,
          label='{x}'.format(x="Enojni resonančni filter"))
plt.axvline(x=12.45, color='r')
plt.ylabel("Premik[nm]",fontsize="25")
plt.legend(fontsize="25")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)    
plt.legend(fontsize="25",loc="lower left")
manager = plt.get_current_fig_manager()
manager.full_screen_toggle()
plt.tight_layout()

if risanje_test_7 == True:    
    plt.savefig('Z_koncna_analiza')  




fig_close_1=plt.figure(100)
plt.close(fig_close_1)

fig_close_2=plt.figure(200)
plt.close(fig_close_2)

fig_close_3=plt.figure(300)
plt.close(fig_close_3)

fig_close_6=plt.figure(400)
plt.close(fig_close_6)

fig_close_4=plt.figure(500)
plt.close(fig_close_4)

fig_close_5=plt.figure(600)
plt.close(fig_close_5)

fig_close_7=plt.figure(700)
plt.close(fig_close_7)












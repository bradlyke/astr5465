import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('text', usetex=True)
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord
import progressBar as pb
from astropy.io import fits
from matplotlib.lines import Line2D

#NOTE: Requires gcdata.fits data file to run properly. Found in my user
#folder.

def rat_func(A,r,wp):
    if wp == '1':
        return 10**(0.2*A*r)
    elif wp == '2':
        return 0.2*A*r

def q1_plot(pchoice):
    x = np.linspace(0,20,1000)
    Av = np.array([0.1,0.5,1.0,1.5])

    y0 = rat_func(Av[0],x,pchoice)
    y1 = rat_func(Av[1],x,pchoice)
    y2 = rat_func(Av[2],x,pchoice)
    y3 = rat_func(Av[3],x,pchoice)

    z0 = y0*(3)
    z1 = y1*(3)
    z2 = y2*(3)
    z3 = y3*(3)

    x0 = y0*(-3)
    x1 = y1*(-3)
    x2 = y2*(-3)
    x3 = y3*(-3)

    fig1,ax1 = plt.subplots(figsize=(8,6))
    ax1.plot(x,y0,linestyle='-',color='black',label='$A_v = 0.1$')
    ax1.plot(x,y1,linestyle='--',color='blue',label='$A_v = 0.5$')
    ax1.plot(x,y2,linestyle='-.',color='orange',label='$A_v = 1.0$')
    ax1.plot(x,y3,linestyle=':',color='green',label='$A_v = 1.5$')
    ax1.tick_params(axis='both',direction='in')
    ax1.set_xticklabels(ax1.get_xticks(),fontsize=15)
    ax1.set_yticklabels(ax1.get_yticks(),fontsize=15)
    axbox = ax1.get_position()
    #plt.legend(fontsize=15,loc=(axbox.x0 + 0.62,axbox.y0 + 0))
    plt.legend(fontsize=15,loc='upper left')
    ax1.set_xlabel('Distance [kpc]',fontsize=15)
    ax1.set_ylabel('$\log[d/r]$',fontsize=15)

    fig2,ax2 = plt.subplots(figsize=(8,6))
    ax2.plot(x,z0,linestyle='-',color='black',label='$A_v = 0.1$')
    ax2.plot(x,z1,linestyle='--',color='blue',label='$A_v = 0.5$')
    ax2.plot(x,z2,linestyle='-.',color='orange',label='$A_v = 1.0$')
    ax2.plot(x,z3,linestyle=':',color='green',label='$A_v = 1.5$')
    ax2.tick_params(axis='both',direction='in')
    ax2.set_xticklabels(ax2.get_xticks(),fontsize=15)
    ax2.set_yticklabels(ax2.get_yticks(),fontsize=15)
    axbox = ax2.get_position()
    #plt.legend(fontsize=15,loc=(axbox.x0 + 0.62,axbox.y0 + 0))
    plt.legend(fontsize=15,loc='upper left')
    ax2.set_xlabel('Distance [kpc]',fontsize=15)
    ax2.set_ylabel('$\log[(d/r)^3]$',fontsize=15)

    fig3,ax3 = plt.subplots(figsize=(8,6))
    ax3.plot(x,x0,linestyle='-',color='black',label='$A_v = 0.1$')
    ax3.plot(x,x1,linestyle='--',color='blue',label='$A_v = 0.5$')
    ax3.plot(x,x2,linestyle='-.',color='orange',label='$A_v = 1.0$')
    ax3.plot(x,x3,linestyle=':',color='green',label='$A_v = 1.5$')
    ax3.tick_params(axis='both',direction='in')
    ax3.set_xticklabels(ax3.get_xticks(),fontsize=15)
    ax3.set_yticklabels(ax3.get_yticks(),fontsize=15)
    axbox = ax3.get_position()
    #plt.legend(fontsize=15,loc=(axbox.x0 + 0.62,axbox.y0 + 0))
    plt.legend(fontsize=15,loc='lower left')
    ax3.set_xlabel('Distance [kpc]',fontsize=15)
    ax3.set_ylabel('$\log[(r/d)^3]$',fontsize=15)

    plt.show()

def q2_plot():

    #Above is the old way to load data.
    gcords = fits.open('gcdata.fits')[1].data
    fig0, ax0 = plt.subplots(figsize=(8,6),subplot_kw=dict(polar=True))
    ax0.scatter(gcords['l'],gcords['rsun'],marker='o',linewidths=1,s=10)
    ax0.set_rmax(150)
    ax0.set_rticks([30,60,90,120,150])
    ax0.set_xticks(np.pi/180 * np.linspace(0,360,4,endpoint=False))
    ax0.set_yticklabels(['30','60','90','120','150$~$kpc'],fontsize=15)
    ax0.set_xticklabels(['0$^\circ$','90$^\circ$','180$^\circ$','270$^\circ$'],fontsize=15)
    ax0.grid(True)

    fig1,ax1 = plt.subplots(figsize=(8,6))
    ax1.hist(gcords['l'],bins=20)
    ax1.tick_params(axis='both',direction='in')
    ax1.set_xticklabels(ax1.get_xticks(),fontsize=15)
    ax1.set_yticklabels(ax1.get_yticks(),fontsize=15)
    #axbox = ax1.get_position()
    #plt.legend(fontsize=15,loc=(axbox.x0 + 0.62,axbox.y0 + 0))
    #plt.legend(fontsize=15,loc='upper left')
    ax1.set_xlabel('Galatic Longitude (deg)',fontsize=15)
    ax1.set_ylabel('Number of Globular Clusters',fontsize=15)

    proj_dist = gcords['rsun']*np.cos(np.deg2rad(gcords['b']))
    fig2,ax2 = plt.subplots(figsize=(8,6))
    ax2.scatter(gcords['b'],proj_dist,marker='o',color='black',s=10)
    ax2.set_xticklabels(ax2.get_xticks(),fontsize=15)
    ax2.set_yticklabels(ax2.get_yticks(),fontsize=15)
    ax2.tick_params(axis='both',direction='in')
    ax2.set_ylabel('Projected Distance (kpc)',fontsize=15)
    ax2.set_xlabel('Galactic Latitude (deg)',fontsize=15)

    plt.show()

def main_menu():
    #while True:
    print('\n')
    print('ASTR 5465 Homework 1 Plots')
    print('Which question do you want to generate plots for')
    print('\n[1] Q1\n[2] Q2\n[3] Quit')
    qq = input('Selection: ')
    if qq == '1':
        print('[1] Linear\n[2] Logarithmic:')
        wloglin = input('Selection: ')
        q1_plot(wloglin)
    elif qq == '2':
        q2_plot()
main_menu()

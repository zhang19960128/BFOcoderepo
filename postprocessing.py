#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 14:50:37 2020

@author: jiahaoz
"""
import numpy as np
import math
import matplotlib.pyplot as plt
class Eenthalpy():
    def __init__(self,path,scfin,highsymm):
        self.scffile=path+'/'+scfin;
        self.path=path;
        self.highfile=path+'/'+highsymm;
        self.obtainvol();
    def obtainvol(self):
        f=open(self.scffile,'r');
        lines=f.readlines();
        length=len(lines);
        self.axis=np.zeros([3,3]);
        atoau=0.529;
        for i in range(length):
            if lines[i].find("CELL_PARAMETERS")!=-1:
                for j in range(3):
                    for k in range(3):
                        self.axis[j][k]=float(lines[i+j+1].split()[k])/atoau;
        f.close();
    def obtaindipole(self,file):
        f=open(file,'r');
        lines=f.readlines();
        length=len(lines);
        edipole=np.zeros(3);
        idipole=np.zeros(3);
        dipoles=np.zeros([2,3]);
        for i in range(length):
            if(lines[i].find("Electronic Dipole on Cartesian axes")!=-1):
                for j in range(3):
                    edipole[j]=float(lines[i+j+1].split()[1]);
            if(lines[i].find("Ionic Dipole on Cartesian axes")!=-1):
                for j in range(3):
                    idipole[j]=float(lines[i+j+1].split()[1]);
        dipoles[0][0:3]=np.copy(edipole*1.0/math.sqrt(2));
        dipoles[1][0:3]=np.copy(idipole*1.0/math.sqrt(2));
        return dipoles;
    def scan(self,start,end,step):
        estep=np.arange(start,end+step,step);
        length=len(estep);
        edipole=np.zeros([length,3]);
        idipole=np.zeros([length,3]);
        totaldipole=np.zeros([length,3]);
        totalpolar=np.zeros([length,3]);
        symdipole=self.obtaindipole(self.highfile);
        for i in range(length):
            dip=self.obtaindipole(self.path+'/'+'bfoE'+"{:3.1f}".format(estep[i])+'.out');
            edipole[i,0:3]=dip[0][0:3]-symdipole[0][0:3];
            idipole[i,0:3]=dip[1][0:3]-symdipole[1][0:3];
            totaldipole[i,0:3]=edipole[i,0:3]+idipole[i,0:3];
            totalpolar[i,0:3]=totaldipole[i,0:3]/np.linalg.det(self.axis)*57.137;
        plt.plot(estep,totalpolar[:,0],marker='o',markersize=8,Linewidth=2)
        plt.plot(estep,totalpolar[:,1],marker='s',markersize=8,Linewidth=2)
        plt.legend(["px","py"])
        plt.xlabel("Efield(mv/cm)")
        plt.ylabel('Polarization(c/m^2)')
en=Eenthalpy('.','bfo_bulk.in','bfo_high_symm.out');
en.scan(-6.0,6.0,0.4)

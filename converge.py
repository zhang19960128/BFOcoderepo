#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 14:50:37 2020

@author: jiahaoz
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
class Eenthalpy():
    def __init__(self,natoms,path,scfin,base):
        self.scffile=path+'/'+scfin;
        self.path=path;
        self.basefile=path+'/'+base;
        self.obtainvol();
        self.natoms=natoms;
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
    def convergeP(self,filelist):
        length=len(filelist);
        edipolelist=np.zeros([length,3]);
        idipolelist=np.zeros([length,3]);
        tpdipolelist=np.zeros([length,3]);
        referdip=self.obtaindipole(self.basefile);
        for i in range(length):
            dip=self.obtaindipole(filelist[i]);
            edipolelist[i][0:3]=dip[0,0:3]-referdip[0,0:3];
            idipolelist[i][0:3]=dip[1,0:3]-referdip[1,0:3];
            tpdipolelist[i][0:3]=(edipolelist[i,0:3]+idipolelist[i,0:3])/np.linalg.det(self.axis)*57.137;
        matplotlib.rc('font',size=20)
        plt.plot(range(length),tpdipolelist[:,2],marker='o',markersize=5,Linewidth=2);
        plt.xlabel("iteration #");
        plt.ylabel("$\Delta P_{z}$");
    def obtainforce(self,file): 
        # post-processing the scf.out
        scfout=open(file,'r');
        lines=scfout.readlines();
        force=np.zeros([self.natoms,3]);
        for i in range(len(lines)):
            if lines[i].find("Forces acting on atoms (cartesian axes, Ry/au):")!=-1:
                for j in range(self.natoms):
                    force[j][0]=float(lines[i+j+2].split()[6]);
                    force[j][1]=float(lines[i+j+2].split()[7]);
                    force[j][2]=float(lines[i+j+2].split()[8]);
        scfout.close();
        return force;
    def convergeF(self,filelist):
        length=len(filelist);
        norm=np.zeros(length);
        for i in range(length):
            force=self.obtainforce(filelist[i]);
            norm[i]=0.0;
            for j in range(self.natoms):
                for k in range(3):
                    if norm[i] < np.sqrt(force[j][k]*force[j][k]):
                        norm[i]=np.sqrt(force[j][k]*force[j][k]);
                    else:
                        norm[i]=norm[i];
            norm[i]=np.sqrt(norm[i]/self.natoms/3);
        matplotlib.rc('font',size=20)
        plt.plot(range(length),norm,marker='o',markersize=5,Linewidth=2);
        plt.xlabel("iteration #");
        plt.ylabel("$|F_{max}|(Ry/au)$");        
    def scan(self,estep):
        length=len(estep);
        edipole=np.zeros([length,3]);
        idipole=np.zeros([length,3]);
        totaldipole=np.zeros([length,3]);
        totalpolar=np.zeros([length,3]);
        symdipole=self.obtaindipole(self.path+'/'+'bfoE'+'0.0'+'.out');
        for i in range(length):
            dip=self.obtaindipole(self.path+'/'+'bfoE'+"{:3.1f}".format(estep[i])+'.out');
            edipole[i,0:3]=dip[0][0:3]-symdipole[0][0:3];
            idipole[i,0:3]=dip[1][0:3]-symdipole[1][0:3];
            totaldipole[i,0:3]=edipole[i,0:3]+idipole[i,0:3];
            totalpolar[i,0:3]=totaldipole[i,0:3]/np.linalg.det(self.axis)*57.137;
        plt.plot(estep,totalpolar[:,0],marker='o',markersize=4,Linewidth=2)
        plt.plot(estep,totalpolar[:,1],marker='s',markersize=4,Linewidth=2)
        plt.plot(estep,totalpolar[:,2],marker='1',markersize=4,Linewidth=2)
        plt.legend(["px","py",'pz'])
        plt.xlabel("Ez(mv/cm)")
        plt.xlim([-2.0,2.0])
        plt.ylabel('Change of Polarization(c/m^2)')
        matplotlib.rc('font',size=20)
en=Eenthalpy(5,'.','./PWOUT/bto.in','./PWOUT/ite.out0');
flist=["./PWOUT/ite.out"+str(i) for i in range(7)];
#en.convergeP(flist)
en.convergeF(flist)
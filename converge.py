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
        self.atommass=np.zeros(natoms);
        self.edie=np.zeros([3,3]);
        self.atomcharge=np.zeros([natoms,3,3]);
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
                        self.axis[j][k]=float(lines[i+j+1].split()[k]);
        f.close();
    def obtainph(self,file):
        phout=open(file,'r');
        #post-processing ph.out
        lines=phout.readlines();
        for i in range(len(lines)):
            if lines[i].find("site n.  atom      mass           positions (alat units)")!=-1:
                for j in range(self.natoms):
                    for k in range(3):
                        self.atommass[j]=float(lines[i+j+1].split()[2]);
            if lines[i].find("Dielectric constant in cartesian axis")!=-1:
                for j in range(3):
                    for k in range(3):
                        self.edie[j][k]=float(lines[i+j+2].split()[k+1]);
            if lines[i].find("Effective charges (d P / du) in cartesian axis")!=-1:
                for j in range(self.natoms):
                    for k in range(3):
                        for m in range(3):
                            self.atomcharge[j][k][m]=lines[i+2+4*j+k+1].split()[m+2];
        phout.close();
    def obtaindipolescftotal(self,file):
        berryout=open(file,'r');
        lines=berryout.readlines();
        epolar=np.zeros(3);
        apolar=np.zeros(3);
        total=np.zeros(3);
        Atoau=0.52917721067121;
        for i in range(len(lines)):
            if(lines[i].find("Electronic Dipole on Cartesian axes")!=-1):
                for j in range(3):
                    epolar[j]=float(lines[i+1+j].split()[1])/math.sqrt(2);
            elif(lines[i].find("Ionic Dipole on Cartesian axes")!=-1):
                for j in range(3):
                    apolar[j]=float(lines[i+1+j].split()[1])/math.sqrt(2);
        return[epolar,apolar];
    def readposition(self,scfin,natoms):
        scfinput=open(scfin,'r');
        lines=scfinput.readlines();
        length=len(lines);
        atomp=np.zeros([natoms,3]);
        for i in range(length):
            if lines[i].find("ATOMIC_POSITIONS")!=-1 and lines[i].find("(angstrom)")!=-1:
                for j in range(natoms):
                    for k in range(3):
                        atomp[j,k]=float(lines[i+j+1].split()[k+1]);
        return atomp;
    def obtaindipolediffperiod(self,pscfin,pscfout,ascfin,ascfout):
        # this function is the function used to calibrate the periodical effects of berry phase:
        print("====================================================");
        Atoau=0.52917721067121;
        previous=self.obtaindipolescftotal(pscfout);
        now=self.obtaindipolescftotal(ascfout);
        period=np.zeros(3);
        for i in range(3):
            period[i]=self.axis[i][i]/Atoau;
        deltae=np.zeros(3);
        deltai=np.zeros(3);
        for i in range(3):
            deltae[i]=now[0][i]-previous[0][i];
            deltai[i]=now[1][i]-previous[1][i];
        f=open(pscfin,'r');
        lines=f.readlines();
        f.close();
        ebefore=np.zeros(3);
        eafter=np.zeros(3);
        for i in range(len(lines)):
            for j in range(3):
                if lines[i].find("efield_cart("+str(j+1)+")")!=-1:
                    ebefore[j]=float(lines[i].split("=")[-1]);
        f=open(ascfin,'r');
        lines=f.readlines();
        f.close();
        for i in range(len(lines)):
            for j in range(3):
                if lines[i].find("efield_cart("+str(j+1)+")")!=-1:
                    eafter[j]=float(lines[i].split("=")[-1]);
        diffe=eafter-ebefore;
        # estimate the electrical polarization change;
        epsilzero=8.8541878*10**-12;
        angstrom=10**-10;
        efieldqe=36.3609*10**10;
        echarge=1.60217662*10**-19;
        au=0.529*10**-10;
        dipoleqe=epsilzero*angstrom**3*efieldqe/echarge/au;
        polarest=np.matmul(self.edie,diffe)*np.linalg.det(self.axis)*dipoleqe;# units e*au;
        print("estimate electric dipole change: by suseptbility ",polarest,'QE output: ',deltae);
        for i in range(3):
            deltae[i]=deltae[i]-round(deltae[i]/period[i])*period[i];
            temp=period[i]*(round(polarest[i]/period[i])-round(deltae[i]/period[i]));
            deltae[i]=deltae[i]+temp;
        print("QE output: ",deltae,"Period: ",period);
        #estimate the ioninc polarization change;
        atombefore=self.readposition(pscfin,self.natoms);
        atomafter=self.readposition(ascfin,self.natoms);
        diffp=atomafter-atombefore;
        dipolestp=np.zeros(3);
        for i in range(self.natoms):
            dipolestp=dipolestp+np.matmul(self.atomcharge[i,0:3,0:3],diffp[i])/0.529;
        print("estimate ionic dipole change: by Born effective charge: ",dipolestp,'QE ouput',deltai);
        for i in range(3):
            deltai[i]=deltai[i]-round(deltai[i]/period[i])*period[i];
            temp=period[i]*(round(dipolestp[i]/period[i])-round(deltai[i]/period[i]));
            deltai[i]=deltai[i]+temp;
        print("QE output: ",deltai,"Period: ",period);
        print("====================================================");
        return [deltae/np.linalg.det(self.axis/Atoau),deltai/np.linalg.det(self.axis/Atoau)];
    def convergeP(self,scfinlist,scfoutlist,phlist):
        length=len(scfinlist);
        delpolare=np.zeros([length,3]);
        delpolari=np.zeros([length,3]);
        total=np.zeros([length,3]);
        for i in range(1,length):
            self.obtainph(phlist[i]);
            re=self.obtaindipolediffperiod(scfinlist[i-1],scfoutlist[i-1],scfinlist[i],scfoutlist[i]);
            delpolare[i,0:3]=np.copy(re[0]);
            delpolari[i,0:3]=np.copy(re[1]);
            for j in range(i+1):
                total[i,0:3]=delpolare[j,0:3]+delpolari[j,0:3]+total[i,0:3];
            total[i,0:3]=57.2148*total[i,0:3];
        print(total)
        matplotlib.rc('font',size=20);
        plt.plot(range(length),total[:,0],marker=7,markersize=10);
        plt.plot(range(length),total[:,1],marker='X',markersize=10);
        plt.plot(range(length),total[:,2],marker=4,markersize=10);
        plt.ylim([-0.10,0.05])
        plt.legend(["Px","Py","Pz"])
        plt.xlabel("iteration #")
        plt.ylabel("Polarization (C/$m^2$)")
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
en=Eenthalpy(5,'.','./PWOUT/ite0','./PWOUT/ite.out0');
iten=6;
scfinlist=["./PWOUT/ite.out"+str(i) for i in range(iten)];
scfoutlist=["./PWOUT/ite.out"+str(i) for i in range(iten)];
phlist=["./PWOUT/ph.out"+str(i) for i in range(iten)];
en.convergeP(scfinlist,scfoutlist,phlist)
#en.convergeF(flist)
import numpy as np
import re
import math
import sys
class pwout:
    def __init__(self,filepath,phout,dynout,zero,scfin):
        self.path=filepath;
        self.phfile=self.path+"/"+phout;
        self.dynfile=self.path+'/'+dynout;
        self.zerofile=self.path+'/'+zero;
        self.scfin=self.path+'/'+scfin;
    def obtain(self,natom):
        au=0.52917721067121;
        self.natoms=natom;
        self.axis=np.zeros([3,3]);
        self.atomp=np.zeros([self.natoms,3]);
        self.atommass=np.zeros([self.natoms,1]);
        self.atomcharge=np.zeros([self.natoms,3,3]);
        self.atomdis=np.arange(self.natoms*3);
        self.dynmatrix=np.zeros([self.natoms,self.natoms,3,3]);
        self.edie=np.zeros([3,3]);
        self.force=np.zeros([self.natoms,3]);
        self.polar=np.zeros(3);
        self.readscf();
        self.obtainph(self.phfile);
        self.obtaindyn(self.dynfile);
        self.obtainforce(self.zerofile);
        self.vol=np.linalg.det(pw.axis)/au/au/au;
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
    def obtaindyn(self,file):
        dyn=open(file);
        # post-processing BFO.dyn
        lines=dyn.readlines();
        for i in range(len(lines)):
            if lines[i].find("Dynamical  Matrix in cartesian axes")!=-1:
                for j in range(self.natoms):
                    for k in range(self.natoms):
                        for m in range(3):
                            self.dynmatrix[j][k][0][m]=float(lines[i+4+4*(j*self.natoms+k)+1].split()[m*2]);
                            self.dynmatrix[j][k][1][m]=float(lines[i+4+4*(j*self.natoms+k)+2].split()[m*2]);                    
                            self.dynmatrix[j][k][2][m]=float(lines[i+4+4*(j*self.natoms+k)+3].split()[m*2]);
        for i in range(self.natoms):
            for j in range(self.natoms):
                #be careful about w and f, those are circle frequence and normal frequency.
                self.dynmatrix[i][j]=self.dynmatrix[i][j];
        dyn.close();
    def obtainforce(self,file): 
        # post-processing the scf.out
        scfout=open(file,'r');
        lines=scfout.readlines();
        for i in range(len(lines)):
            if lines[i].find("Forces acting on atoms (cartesian axes, Ry/au):")!=-1:
                for j in range(self.natoms):
                    self.force[j][0]=float(lines[i+j+2].split()[6]);
                    self.force[j][1]=float(lines[i+j+2].split()[7]);
                    self.force[j][2]=float(lines[i+j+2].split()[8]);
        scfout.close();
    def obtaindipolescf(self,file):
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
        for i in range(3):
            total[i]=epolar[i]+apolar[i];
        # now the units of dipole is e * bohr   
        print('the polarization now is: ',total/np.linalg.det(self.axis/Atoau))
        return total/np.linalg.det(self.axis/Atoau);
    def readscf(self):
        scfinput=open(self.scfin,'r');
        lines=scfinput.readlines();
        length=len(lines);
        for i in range(length):
            if lines[i].find("CELL_PARAMETERS")!=-1 and lines[i].find("(angstrom)")!=-1:
                for j in range(3):
                    for k in range(3):
                        self.axis[j][k]=float(lines[i+j+1].split()[k]);
            if lines[i].find("ATOMIC_POSITIONS")!=-1 and lines[i].find("(angstrom)")!=-1:
                for j in range(self.natoms):
                    for k in range(3):
                        self.atomp[j,k]=float(lines[i+j+1].split()[k+1]);
    def writenewscf(self,Efield,atomposition,filename):
        # the Efield units now is MV/cm
        # the atomposition units is the angstrom;
        changeunits=Efield*10**6/10**(-2)/(36.3509*10**10);
        scffiles=open(self.scfin,'r');
        newfilename=open(self.path+'/'+filename,'w');
        lines=scffiles.readlines();
        tick1=0;
        tick2=0;
        for i in range(len(lines)):
            if lines[i].find("efield_cart")!=-1:
                tick1=tick1+1;
                if tick1>1:
                    continue;
                else:
                    for j in range(3):
                        newfilename.write("efield_cart"+"("+str(j+1)+")="+str(changeunits[j])+"\n");
            elif lines[i].find("ATOMIC_POSITIONS (angstrom)")!=-1:
                newfilename.write("ATOMIC_POSITIONS (angstrom)\n");
                tick2=i;
                for j in range(self.natoms):
                    temp="";
                    for t in range(3):
                        temp=temp+" "+'{:12.8f}'.format(atomposition[j,t]);
                    newfilename.write(lines[j+i+1].split()[0]+" "+temp+"\n");
            elif tick2+self.natoms >= i and tick2 >1:
                continue;
            elif lines[i].find("verbosity")!=-1:
                continue;
            elif lines[i].find("wf_collect")!=-1:
                continue;
            else:
                newfilename.write(lines[i]);
        newfilename.close();
        scffiles.close();
    def writenewscfnoe(self,atomposition,filename):
        scffiles=open(self.scfin,'r');
        newfilename=open(self.path+'/'+filename,'w');
        lines=scffiles.readlines();
        tick1=0;
        tick2=0;
        for i in range(len(lines)):
            if lines[i].find("efield_cart")!=-1:
                continue;
            elif lines[i].find("lelfield")!=-1:
                continue;
            elif lines[i].find("ATOMIC_POSITIONS (angstrom)")!=-1:
                newfilename.write("ATOMIC_POSITIONS (angstrom)\n");
                tick2=i;
                for j in range(self.natoms):
                    temp="";
                    for t in range(3):
                        temp=temp+'{:12.8f}'.format(atomposition[j,t]);
                    newfilename.write(lines[j+i+1].split()[0]+" "+temp+'\n');
            elif tick2+self.natoms >= i and tick2 >1:
                continue;
            else:
                newfilename.write(lines[i]);
        newfilename.close();
        scffiles.close();        
    def iterate(self,times,step):
        startP=self.obtaindipolescf(self.zerofile); #units e/bohr^2;
        endP=startP+np.array([0,0,1])*step/57.137;#units e/bohr^2;step units is C/m^2
        deltaP=startP-endP;
        print('the aim is: ',endP,'Please also be notifying that forces should also be zero');
        self.obtainph(self.phfile);
        self.obtaindyn(self.dynfile);
        # starting the first guess: think the polarization increase partly come from the electrical contribution
        efield=np.matmul(-1*np.linalg.inv(self.edie*self.vol),deltaP*self.vol)/0.0000154961;
        efield=efield*1.0;#only count as 10% from the electrical contribution, the first step.
        print('electrical field=',efield)
        au=0.52917721067121;
        self.writenewscf(efield,self.atomp,"ite1");
        self.writenewscfnoe(self.atomp,"itenoe1");
        for i in range(1,times):
            self.obtainforce(self.path+'/'+'ite.out'+str(i));
            Pnow=self.obtaindipolescf(self.path+'/'+'ite.out'+str(i));
            self.obtainph(self.path+'/'+'ph.out'+str(i));
            self.obtaindyn(self.path+'/'+'dyn.out'+str(i));
            deltaP=Pnow-endP;
            deltax=self.solve(self.force,deltaP);
            posit=au*np.reshape(deltax[0:self.natoms*3],[self.natoms,3]);
            posit=posit+self.atomp;
            self.writenewscf(deltax[self.natoms*3:(self.natoms+1)*3],posit,'ite'+str(i+1));
            self.writenewscfnoe(posit, 'itenoe'+str(i+1));
    def iteratetwo(self,times,step):
        startP=self.obtaindipolescf(self.zerofile); #units e/bohr^2;
        endP=startP+np.array([0,0,1])*step/57.137;#units e/bohr^2;step units is C/m^2
        deltaP=startP-endP;
        print('the aim is: ',endP,'Please also be notifying that forces should also be zero');
        self.obtainph(self.phfile);
        self.obtaindyn(self.dynfile);
        # starting the first guess: think the polarization increase partly come from the electrical contribution
        efield=np.matmul(-1*np.linalg.inv(self.edie*self.vol),deltaP*self.vol)/0.0000154961;
        efieldzero=efield*1.0;#only count as 10% from the electrical contribution, the first step.
        print('electrical zero field=',efieldzero)
        au=0.52917721067121;
        self.writenewscf(efield,self.atomp,"ite1");
        self.writenewscfnoe(self.atomp,"itenoe1");
        accup=np.zeros([self.natoms,3])
        accue=np.zeros(3);
        for i in range(1,times):
            self.obtainforce(self.path+'/'+'ite.out'+str(i));
            Pnow=self.obtaindipolescf(self.path+'/'+'ite.out'+str(i));
            self.obtainph(self.path+'/'+'ph.out'+str(i));
            self.obtaindyn(self.path+'/'+'dyn.out'+str(i));
            deltaP=Pnow-endP;
            deltax=self.solve(self.force,deltaP);
            accup=accup+au*np.reshape(deltax[0:self.natoms*3],[self.natoms,3]);
            posit=accup+self.atomp;
            accue=accue+deltax[self.natoms*3:(self.natoms+1)*3];
            efield=accue+efieldzero;
            print('deltaE=',deltax[self.natoms*3:(self.natoms+1)*3],'accue=',efield)
            self.writenewscf(efield,posit,'ite'+str(i+1));
            self.writenewscfnoe(posit, 'itenoe'+str(i+1));
    def solve(self,force,deltaP):
        #solve the linear equation Ax=y;
        self.A=np.zeros([3*self.natoms+3,3*self.natoms+3]);
        self.x=np.zeros(3*self.natoms+3);
        self.y=np.zeros(3*self.natoms+3);
        for i in range(self.natoms):
            for j in range(self.natoms):
                self.A[3*i:3*(i+1),3*j:3*(j+1)]=np.copy(self.dynmatrix[i,j,0:3,0:3]);
        for i in range(self.natoms):
            self.A[3*i:3*(i+1),3*self.natoms:3*(self.natoms+1)]=np.copy(-1*self.atomcharge[i,0:3,0:3]*0.0003884098);
        for i in range(self.natoms):
            self.A[3*(self.natoms):3*(self.natoms+1),3*i:3*(i+1)]=np.copy(-1*self.atomcharge[i,0:3,0:3]);
        self.A[3*(self.natoms):3*(self.natoms+1),3*self.natoms:3*(self.natoms+1)]=np.copy(-1*self.edie[0:3,0:3]*0.0000154961*self.vol);
        for i in range(self.natoms):
            for j in range(3):
                self.y[i*3+j]=force[i][j];
        self.y[3*self.natoms+0]=self.vol*deltaP[0];
        self.y[3*self.natoms+1]=self.vol*deltaP[1];
        self.y[3*self.natoms+2]=self.vol*deltaP[2];
        self.x=np.matmul(np.linalg.inv(self.A),self.y);
        return np.matmul(np.linalg.inv(self.A),self.y);
pw=pwout("./PWOUT","ph.out1",'bto.dyn1','btozero.out','bto.in');
pw.obtain(5);
pw.iteratetwo(int(sys.argv[1]),    float(sys.argv[2]));
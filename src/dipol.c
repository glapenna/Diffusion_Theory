#include "head.h"

void
dipol(idstruct *serptr)
/*************************************************************************
   It evaluates electrostatic dipolar properties
*************************************************************************/
{
	_REAL   *r = serptr->idunion.ser.set->crd;
        parmstruct *prm = serptr->idunion.ser.prm;
	char atomname[ATOMTOKLEN];

	int i,ia,im,jm,count,natom,nmol;
        _REAL Mx,My,Mz,mx1,my1,mz1,mu,Msq,mu_mod,
              xc,yc,zc,q;
/*
        _REAL m1[500][3],S;
*/

        Mx = 0.; My = 0.; Mz = 0.;
        mu = 0.;

/*
	r+=3*prm->Nsolu_tot;
	prm->Charge+=prm->Nsolu_tot;
*/
	natom=prm->Natom;
	nmol=prm->Nmol;

        for(im=0; im<nmol; im++) {
         mx1 = 0.; my1 = 0.; mz1 = 0.;
         xc = 0.; yc = 0.; zc = 0.;
	 count=0;
         for(ia=0; ia<natom; ia++,r+=3) {
	  strncpy(atomname,&prm->AtomNames[ia*ATOMTOKLEN],ATOMTOKLEN);
	  if( !strncmp(atomname,"CL-",3))
          {
          /* printf("atomname(%d)= %s\n",ia,atomname); */
          /* printf("r0= %e\n",r[0]); */
	  q=prm->Charge[ia];
          mx1 += q*r[0];
          my1 += q*r[1];
          mz1 += q*r[2];
          xc += r[0];
          yc += r[1];
          zc += r[2];
          count++;
          }
         }
	 xc/=(_REAL)count;
	 yc/=(_REAL)count;
	 zc/=(_REAL)count;
	 mx1-=xc;
	 my1-=yc;
	 mz1-=zc;
         Mx  += mx1;
         My  += my1;
         Mz  += mz1;
        /*  if(im==1) printf("m1: %e %e %e\n",mx1,my1,mz1); */
/*
         printf("ATOM  %5d  S   WAT %5d    %8.3f%8.3f%8.3f\n",
          1,im+1,0.,0.,0.);
         printf("ATOM  %5d  Xe  WAT %5d    %8.3f%8.3f%8.3f\n",
          2,im+1,mx1,my1,mz1);
         printf("TER\n");
*/
         mu_mod = l_sqrt(mx1*mx1+my1*my1+mz1*mz1);
/*
         m1[im][0]=mx1/mu_mod;
         m1[im][1]=my1/mu_mod;
         m1[im][2]=mz1/mu_mod;
*/

         mu += mu_mod;
        }
        mu /= (_REAL)prm->Nmol;
        r-=3*prm->Nsit;

        Msq = Mx*Mx+My*My+Mz*Mz;

        printf("dipol: %e %e %e %e %e\n",mu,Msq,Mx,My,Mz);
        fflush(stdout);

/*
        S = 0.;
        i = 1;
        for(im=0; im<prm->Nmol-1; im++) {
         for(jm=im+1; jm<prm->Nmol; jm++,i++) {
          S+= m1[im][0]*m1[jm][0]
             +m1[im][1]*m1[jm][1]
             +m1[im][2]*m1[jm][2];
         }
        }
        S/=(_REAL)i;
        printf("P1: %e\n",S);
        fflush(stdout);
*/
}

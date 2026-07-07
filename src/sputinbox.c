#include "head.h"

void
sputinbox(idstruct *serptr)
/************************************************************************/
/* rebuild the molecule */
/************************************************************************/
{
	_REAL   *r = serptr->idunion.ser.set->crd;
	_REAL   *box = serptr->idunion.ser.set->boxc;
	_REAL   **c = serptr->idunion.ser.set->cell_to_cart;
	_REAL   **c_1 = serptr->idunion.ser.set->cart_to_cell;
	parmstruct *prm=serptr->idunion.ser.prm;;

	int i,j,ia,im,nsit;
        _REAL x0,y0,z0,r1[3],r2[3],dr[3],drp[3],lim;

        lim=0.70;
        nsit = prm->Nsolu_tot;

     printf("Sto ricostruendo con lim= %f\n",lim);
     printf("natom= %d\n",prm->Natom_solu);

        for(im=0; im<prm->Nsolu; im++) {
         for(ia=0; ia<prm->Natom_solu-1; ia++) {
/*
     printf("r(%d)= %e %e %e %e %e %e ",ia,r[0],r[1],r[2],r[3],r[4],r[5]);
     printf("%e %e %e\n",r[6],r[7],r[8]);
*/

          x0 = *r++; y0 = *r++; z0 = *r++;
          r1[0] = *r++; r1[1] = *r++; r1[2] = *r++;
          dr[0]=r1[0]-x0;
          dr[1]=r1[1]-y0;
          dr[2]=r1[2]-z0;
          r2[0]=c_1[1][1]*r1[0]+c_1[1][2]*r1[1]+c_1[1][3]*r1[2];
          r2[1]=c_1[2][1]*r1[0]+c_1[2][2]*r1[1]+c_1[2][3]*r1[2];
          r2[2]=c_1[3][1]*r1[0]+c_1[3][2]*r1[1]+c_1[3][3]*r1[2];
          drp[0]=c_1[1][1]*dr[0]+c_1[1][2]*dr[1]+c_1[1][3]*dr[2];
          drp[1]=c_1[2][1]*dr[0]+c_1[2][2]*dr[1]+c_1[2][3]*dr[2];
          drp[2]=c_1[3][1]*dr[0]+c_1[3][2]*dr[1]+c_1[3][3]*dr[2];
          if(drp[0]>lim) {
           r2[0]-=1.;
          } else if(drp[0]<-lim) {
           r2[0]+=1.;
          }
          if(drp[1]>lim) {
           r2[1]-=1.;
          } else if(drp[1]<-lim) {
           r2[1]+=1.;
          }
          if(drp[2]>lim) {
           r2[2]-=1.;
          } else if(drp[2]<-lim) {
           r2[2]+=1.;
          }
          r1[0]=c[1][1]*r2[0]+c[1][2]*r2[1]+c[1][3]*r2[2];
          r1[1]=c[2][1]*r2[0]+c[2][2]*r2[1]+c[2][3]*r2[2];
          r1[2]=c[3][1]*r2[0]+c[3][2]*r2[1]+c[3][3]*r2[2];
          r-=3;
          *r++ = r1[0];
          *r++ = r1[1];
          *r++ = r1[2];
          r-=3;
         }
         r+=3;
        }
        r-=3*nsit;
        serptr->idunion.ser.set->crd = r;
}

#include "head.h"

void euler();

/****************************************************************************/

_REAL innerpro_ext(_REAL *p)
{
      _REAL ux[3],uy[3],uz[3],ut[3],ang[3],mod,q[4];
      int   k;
      _REAL alpha,beta,gamma,beh,alefp,alefm,behchh,
            scp;
      fcomplex d0,d0t,d1,d1t,D20,D20t,D10,D10t,D00,tmp1,tmp2,tmp3,tmp4;

/* define q nh bond(NH) */
/* find the frame with z along NH */
         uz[0]= p[0];
         uz[1]= p[1];
         uz[2]= p[2];
         mod = l_sqrt(uz[0]*uz[0]+uz[1]*uz[1]+uz[2]*uz[2]);
	 for (k=0;k<3;k++) uz[k] /= mod;
         ut[0]= 1.;
         ut[1]= 0.;
         ut[2]= 0.;

         euler(uz,ut,ux,uy,ang);
       
         beh = 0.5*ang[1];
         alefp = 0.5*(ang[0]+ang[2]);
         alefm = 0.5*(ang[0]-ang[2]);
         q[0] =  cos(beh)*cos(alefp);
         q[1] =  sin(beh)*cos(alefm);
         q[2] =  sin(beh)*sin(alefm);
         q[3] =  cos(beh)*sin(alefp);
      
/*
       for (k=0,mod=0.;k<4;mod+=q[k]*q[k],k++) printf(" %e", q[k]);
       printf("\n");
       printf(" qsq = %e\n",mod);
*/

         d0.r =  q[0];
         d0.i = -q[3];
         d0t.r =  q[0];
         d0t.i =  q[3];
         d1.r =  q[2];
         d1.i = -q[1];
         d1t.r =  q[2];
         d1t.i =  q[1];
         D20.r = 0.;
         D20.i = 0.;
         D10.r = 0.;
         D10.i = 0.;
         D00.r = 0.;
         D00.i = 0.;

         tmp1 = Cmul(d1,d1);
         tmp2 = Cmul(d0,d0);
         tmp3 = Cmul(tmp1,tmp2);
         tmp4 = RCmul(2.449515,tmp3);
         D20 = Cadd(D20,tmp4);

         mod = l_sqrt(8./3.)*sin(ang[1])*sin(ang[1])*cos(-2.*ang[0]);
         printf("D20,mod: %e %e\n",D20.r,mod);
         mod = l_sqrt(8./3.)*sin(ang[1])*sin(ang[1])*sin(-2.*ang[0]);
         printf("D20,mod: %e %e\n",D20.i,mod);

         tmp1 = Cmul(d1t,d1t);
         tmp2 = Cmul(d1,d0);
         tmp3 = Cmul(tmp1,tmp2);
         tmp4 = RCmul(2.449515,tmp3); 
         D10 = Cadd(D10,tmp4);
         tmp1 = Cmul(d0t,d1t);
         tmp2 = Cmul(d0,d0);
         tmp3 = Cmul(tmp1,tmp2);
         tmp4 = RCmul(-2.449515,tmp3); 
         D10 = Cadd(D10,tmp4);

         mod = l_sqrt(3./2.)*cos(ang[1])*sin(ang[1])*cos(ang[0]);
         printf("D10,mod: %e %e\n",D10.r,mod);

         tmp1 = Cmul(d1t,d1t);
         tmp2 = Cmul(d1,d1);
         tmp3 = Cmul(tmp1,tmp2);
         tmp4 =  RCmul(1.000000,tmp3);
         D00 = Cadd(D00,tmp4);
         tmp1 = Cmul(d0t,d1t);
         tmp2 = Cmul(d1,d0);
         tmp3 = Cmul(tmp1,tmp2);
         tmp4 =  RCmul(-4.000000,tmp3);
         D00 = Cadd(D00,tmp4);
         tmp1 = Cmul(d0t,d0t);
         tmp2 = Cmul(d0,d0);
         tmp3 = Cmul(tmp1,tmp2);
         tmp4 = RCmul(1.000000,tmp3);
         D00 = Cadd(D00,tmp4);

         mod = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
         mod = 0.5*(3.*mod*mod-1.);
         printf("D00,mod: %e %e\n",D00.r,mod);

         scp = 0.;
         D20t = Conjg(D20);
         tmp1 = Cmul(D20t,D20);
         tmp2 = RCmul(2.,tmp1);
         scp += tmp2.r;
         D10t = Conjg(D10);
         tmp1 = Cmul(D10t,D10);
         tmp2 = RCmul(-2.,tmp1);
         scp += tmp2.r;
         tmp1 = Cmul(D00,D00);
         scp += tmp1.r;
         scp *= l_sqrt(5.);
         return(scp);
}


void
genquatch2(_REAL *p0,_REAL *p1,_REAL *p2,_REAL *dir,
           _REAL qfin1[],_REAL qfin2[])
/*
genquatch2(_REAL *p0,_REAL *p1,_REAL *p2,_REAL *p3)
*/
{
      _REAL ux[3],uy[3],uz[3],ut[3],ang[3],mod,q[4],qdirt[4],
            q1temp[4],q2temp[4],qxtemp[4],qytemp[4],qtemp[4],
            q1[4],q2[4],qx[4],qy[4],qz[4],uz0[4],uzt[4],q1t[4],rch;
      int   k;
      _REAL alpha,beta,gamma,beh,alefp,alefm,behchh;

/* define q ch bond(ch2) */
      rch = 1.09;
      behchh = 0.5*(_PI*109.47/180.0);

/* z on x axis */
      alpha=0.5*_PI;
      beta=0.5*_PI;
      gamma=0.;
      beh = 0.5*beta;
      alefp = 0.5*(alpha+gamma);
      alefm = 0.5*(alpha-gamma);
      qxtemp[0] =  cos(beh)*cos(alefp);
      qxtemp[1] =  sin(beh)*cos(alefm);
      qxtemp[2] =  sin(beh)*sin(alefm);
      qxtemp[3] =  cos(beh)*sin(alefp);

/* z on y axis */
      alpha=_PI;
      beta=0.5*_PI;
      gamma=0.;
      beh = 0.5*beta;
      alefp = 0.5*(alpha+gamma);
      alefm = 0.5*(alpha-gamma);
      qytemp[0] =  cos(beh)*cos(alefp);
      qytemp[1] =  sin(beh)*cos(alefm);
      qytemp[2] =  sin(beh)*sin(alefm);
      qytemp[3] =  cos(beh)*sin(alefp);

/* H1 */
      alpha=_PI+behchh;
      beta=0.5*_PI;
      gamma=0.;
      beh = 0.5*beta;
      alefp = 0.5*(alpha+gamma);
      alefm = 0.5*(alpha-gamma);
      q1temp[0] =  cos(beh)*cos(alefp);
      q1temp[1] =  sin(beh)*cos(alefm);
      q1temp[2] =  sin(beh)*sin(alefm);
      q1temp[3] =  cos(beh)*sin(alefp);

/* H2 */
      alpha=_PI-behchh;
      beta=0.5*_PI;
      gamma=0.;
      beh = 0.5*beta;
      alefp = 0.5*(alpha+gamma);
      alefm = 0.5*(alpha-gamma);
      q2temp[0] =  cos(beh)*cos(alefp);
      q2temp[1] =  sin(beh)*cos(alefm);
      q2temp[2] =  sin(beh)*sin(alefm);
      q2temp[3] =  cos(beh)*sin(alefp);

/* find the frame rigid with the 3-2 bond */
/*
*/
         uz[0]= p2[0] - p1[0];
         uz[1]= p2[1] - p1[1];
         uz[2]= p2[2] - p1[2];
/* ... or the frame rigid with the 3-1 bond */
/*
         uz[0]= p2[0] - p0[0];
         uz[1]= p2[1] - p0[1];
         uz[2]= p2[2] - p0[2];
*/
         mod = l_sqrt(uz[0]*uz[0]+uz[1]*uz[1]+uz[2]*uz[2]);
	 for (k=0;k<3;k++) uz[k] = uz[k]/mod;
         ut[0]= p0[0] - p1[0];
         ut[1]= p0[1] - p1[1];
         ut[2]= p0[2] - p1[2];
         mod = l_sqrt(ut[0]*ut[0]+ut[1]*ut[1]+ut[2]*ut[2]);
	 for (k=0;k<3;k++) ut[k] = ut[k]/mod;

         euler(uz,ut,ux,uy,ang);
       
         beh = 0.5*ang[1];
         alefp = 0.5*(ang[0]+ang[2]);
         alefm = 0.5*(ang[0]-ang[2]);
         q[0] =  cos(beh)*cos(alefp);
         q[1] =  sin(beh)*cos(alefm);
         q[2] =  sin(beh)*sin(alefm);
         q[3] =  cos(beh)*sin(alefp);

/* find q lab(dir) */
/*
*/
         uz[0]= dir[0];
         uz[1]= dir[1];
         uz[2]= dir[2];

         ut[0]=  0.;
         ut[1]=  0.;
         ut[2]= -1.;

         euler(uz,ut,ux,uy,ang);

         beh = 0.5*ang[1];
         alefp = 0.5*(ang[0]+ang[2]);
         alefm = 0.5*(ang[0]-ang[2]);
         qdirt[0] =  cos(beh)*cos(alefp);
         qdirt[1] = -sin(beh)*cos(alefm);
         qdirt[2] = -sin(beh)*sin(alefm);
         qdirt[3] = -cos(beh)*sin(alefp);

/* take the product qCH*qCH2*qdirt */
/*
*/
	 proquat(&q1temp[0],&q[0],&qtemp[0]);
	 proquat(&qtemp[0],&qdirt[0],&q1[0]);
	 proquat(&q2temp[0],&q[0],&qtemp[0]);
	 proquat(&qtemp[0],&qdirt[0],&q2[0]);
	 proquat(&qytemp[0],&q[0],&qtemp[0]);
	 proquat(&qtemp[0],&qdirt[0],&qy[0]);
	 proquat(&qxtemp[0],&q[0],&qtemp[0]);
	 proquat(&qtemp[0],&qdirt[0],&qx[0]);
	 proquat(&q[0],&qdirt[0],&qz[0]);
       for (k=0;k<4;k++) qfin1[k]=qz[k];
       for (k=0;k<4;k++) qfin2[k]=qz[k];

/*
       for (k=0;k<4;k++) printf(" %e", qdirt[k]);
       printf("\n");
*/

/* take the product qCH*qCH2 */
         proquat(q1temp,q,qtemp);
         uz0[0]= 0.;
         uz0[1]= 0.;
         uz0[2]= 0.;
         uz0[3]= 1.;
         q1[0]= qtemp[0];
         q1[1]= qtemp[1];
         q1[2]= qtemp[2];
         q1[3]= qtemp[3];
         q1t[0]= q1[0];
         q1t[1]=-q1[1];
         q1t[2]=-q1[2];
         q1t[3]=-q1[3];
         proquat(q1t,uz0,qtemp);
         proquat(qtemp,q1,uzt);
/*
 printf("ATOM  %5d  O4  MOL %5d    %8.3f%8.3f%8.3f\n",1,1,p0[0],p0[1],p0[2]);
 printf("ATOM  %5d  C1' MOL %5d    %8.3f%8.3f%8.3f\n",2,1,p1[0],p1[1],p1[2]);
 printf("ATOM  %5d  N   MOL %5d    %8.3f%8.3f%8.3f\n",3,1,p2[0],p2[1],p2[2]);
*/
/*
         p3[0]=p1[0]+rch*uzt[1];
         p3[1]=p1[1]+rch*uzt[2];
         p3[2]=p1[2]+rch*uzt[3];
*/
/*
         p3[0]=uzt[1];
         p3[1]=uzt[2];
         p3[2]=uzt[3];
*/


}

void genquatch3(_REAL *p0,_REAL *p1,_REAL *p2,
                _REAL *dir,_REAL qfin[])
{
      _REAL ux[3],uy[3],uz[3],ut[3],ang[3],mod,q[4],qdirt[4],
            q1temp[4],q2temp[4],q3temp[4],q0temp[4],qtemp[4],
            q1[4],q2[4],q3[4],qy[4],qz[4];
      int   k;
      _REAL alpha,beta,gamma,beh,alefp,alefm,behchh,betetr;

/* define q ch bond(ch2) */
      behchh = 0.5*(_PI*109.47/180.0);
      betetr = acos(1./3.);

/* z on y axis */
      alpha=0.;
      beta=-HALFPI;
      gamma=0.;
      beh = 0.5*beta;
      alefp = 0.5*(alpha+gamma);
      alefm = 0.5*(alpha-gamma);
      q0temp[0] =  cos(beh)*cos(alefp);
      q0temp[1] =  sin(beh)*cos(alefm);
      q0temp[2] =  sin(beh)*sin(alefm);
      q0temp[3] =  cos(beh)*sin(alefp);

/* H1 */
      alpha=0.;
      beta=-betetr;
      gamma=0.;
      beh = 0.5*beta;
      alefp = 0.5*(alpha+gamma);
      alefm = 0.5*(alpha-gamma);
      q1temp[0] =  cos(beh)*cos(alefp);
      q1temp[1] =  sin(beh)*cos(alefm);
      q1temp[2] =  sin(beh)*sin(alefm);
      q1temp[3] =  cos(beh)*sin(alefp);

/* H2 */
      alpha=-2.*_PI/3.;
      beta=-betetr;
      gamma=0.;
      beh = 0.5*beta;
      alefp = 0.5*(alpha+gamma);
      alefm = 0.5*(alpha-gamma);
      q2temp[0] =  cos(beh)*cos(alefp);
      q2temp[1] =  sin(beh)*cos(alefm);
      q2temp[2] =  sin(beh)*sin(alefm);
      q2temp[3] =  cos(beh)*sin(alefp);

/* H3 */
      alpha= 2.*_PI/3.;
      beta=-betetr;
      gamma=0.;
      beh = 0.5*beta;
      alefp = 0.5*(alpha+gamma);
      alefm = 0.5*(alpha-gamma);
      q3temp[0] =  cos(beh)*cos(alefp);
      q3temp[1] =  sin(beh)*cos(alefm);
      q3temp[2] =  sin(beh)*sin(alefm);
      q3temp[3] =  cos(beh)*sin(alefp);

/* find the frame rigid with a ch3 */
         uz[0]= p2[0] - p1[0];
         uz[1]= p2[1] - p1[1];
         uz[2]= p2[2] - p1[2];
         mod = size(uz);
	 for (k=0;k<3;k++) uz[k] = uz[k]/mod;
         ut[0]= p0[0] - p1[0];
         ut[1]= p0[1] - p1[1];
         ut[2]= p0[2] - p1[2];
         mod = size(ut);
	 for (k=0;k<3;k++) ut[k] = ut[k]/mod;

/* find q ch3(lab) */
         euler(uz,ut,ux,uy,ang);
       
         beh = 0.5*ang[1];
         alefp = 0.5*(ang[0]+ang[2]);
         alefm = 0.5*(ang[0]-ang[2]);
         q[0] =  cos(beh)*cos(alefp);
         q[1] =  sin(beh)*cos(alefm);
         q[2] =  sin(beh)*sin(alefm);
         q[3] =  cos(beh)*sin(alefp);
      
/* find q lab(dir) */
         uz[0]= dir[0];
         uz[1]= dir[1];
         uz[2]= dir[2];

         ut[0]=  0.;
         ut[1]=  0.;
         ut[2]= -1.;

         euler(uz,ut,ux,uy,ang);

         beh = 0.5*ang[1];
         alefp = 0.5*(ang[0]+ang[2]);
         alefm = 0.5*(ang[0]-ang[2]);
         qdirt[0] =  cos(beh)*cos(alefp);
         qdirt[1] = -sin(beh)*cos(alefm);
         qdirt[2] = -sin(beh)*sin(alefm);
         qdirt[3] = -cos(beh)*sin(alefp);

/* take the product qCH*qCH2*qdirt */
	 proquat(&q1temp[0],&q[0],&qtemp[0]);
	 proquat(&qtemp[0],&qdirt[0],&q1[0]);
	 proquat(&q2temp[0],&q[0],&qtemp[0]);
	 proquat(&qtemp[0],&qdirt[0],&q2[0]);
	 proquat(&q3temp[0],&q[0],&qtemp[0]);
	 proquat(&qtemp[0],&qdirt[0],&q3[0]);
	 proquat(&q0temp[0],&q[0],&qtemp[0]);
	 proquat(&qtemp[0],&qdirt[0],&qy[0]);
	 proquat(&q[0],&qdirt[0],&qz[0]);
       for (k=0;k<4;k++) qfin[k]=qz[k];

/*
       for (k=0;k<4;k++) printf(" %e", qfin[k]);
       printf("\n");
*/

}

void
genquatch1(_REAL *p0,_REAL *p1,_REAL *p2,_REAL *p3)
/* p3 will be filled with H2 position */
{
      _REAL ux[3],uy[3],uz[3],ut[3],ang[3],mod,q[4],
            q1temp[4],q2temp[4],q0temp[4],qtemp[4],
            q1[4],q1t[4],uz0[4],uzt[4];
      int   k;
      _REAL alpha,beta,gamma,beh,alefp,alefm,behchh,rch;

/* define q ch bond(ch2) */
      rch = 1.09;
      behchh = 0.5*(_PI*0./180.0);

/* z on y axis */
      alpha=0.;
      beta=0.5*_PI;
      gamma=0.;
      beh = 0.5*beta;
      alefp = 0.5*(alpha+gamma);
      alefm = 0.5*(alpha-gamma);
      q0temp[0] =  cos(beh)*cos(alefp);
      q0temp[1] =  sin(beh)*cos(alefm);
      q0temp[2] =  sin(beh)*sin(alefm);
      q0temp[3] =  cos(beh)*sin(alefp);

/* H1 */
      alpha=_PI+behchh;
      beta=0.5*_PI;
      gamma=0.;
      beh = 0.5*beta;
      alefp = 0.5*(alpha+gamma);
      alefm = 0.5*(alpha-gamma);
      q1temp[0] =  cos(beh)*cos(alefp);
      q1temp[1] =  sin(beh)*cos(alefm);
      q1temp[2] =  sin(beh)*sin(alefm);
      q1temp[3] =  cos(beh)*sin(alefp);

/* H2 */
      alpha=_PI-behchh;
      beta=0.5*_PI;
      gamma=0.;
      beh = 0.5*beta;
      alefp = 0.5*(alpha+gamma);
      alefm = 0.5*(alpha-gamma);
      q2temp[0] =  cos(beh)*cos(alefp);
      q2temp[1] =  sin(beh)*cos(alefm);
      q2temp[2] =  sin(beh)*sin(alefm);
      q2temp[3] =  cos(beh)*sin(alefp);

/* find the frame rigid with a ch2 */
         uz[0]= p2[0] - p0[0];
         uz[1]= p2[1] - p0[1];
         uz[2]= p2[2] - p0[2];
         mod = l_sqrt(uz[0]*uz[0]+uz[1]*uz[1]+uz[2]*uz[2]);
	 for (k=0;k<3;k++) uz[k] = uz[k]/mod;
         ut[0]= p2[0] - p1[0];
         ut[1]= p2[1] - p1[1];
         ut[2]= p2[2] - p1[2];
         mod = l_sqrt(ut[0]*ut[0]+ut[1]*ut[1]+ut[2]*ut[2]);
	 for (k=0;k<3;k++) ut[k] = ut[k]/mod;

/* find q ch2(lab) */
         euler(uz,ut,ux,uy,ang);
       
         beh = 0.5*ang[1];
         alefp = 0.5*(ang[0]+ang[2]);
         alefm = 0.5*(ang[0]-ang[2]);
         q[0] =  cos(beh)*cos(alefp);
         q[1] =  sin(beh)*cos(alefm);
         q[2] =  sin(beh)*sin(alefm);
         q[3] =  cos(beh)*sin(alefp);

/* take the product qCH*qCH2 */
	 proquat(q1temp,q,qtemp);
         uz0[0]= 0.;
         uz0[1]= 0.;
         uz0[2]= 0.;
         uz0[3]= 1.;
         q1[0]= qtemp[0];
         q1[1]= qtemp[1];
         q1[2]= qtemp[2];
         q1[3]= qtemp[3];
         q1t[0]= q1[0];
         q1t[1]=-q1[1];
         q1t[2]=-q1[2];
         q1t[3]=-q1[3];
	 proquat(q1t,uz0,qtemp);
	 proquat(qtemp,q1,uzt);
/*
 printf("ATOM  %5d  O4  MOL %5d    %8.3f%8.3f%8.3f\n",1,1,p0[0],p0[1],p0[2]);
 printf("ATOM  %5d  C1' MOL %5d    %8.3f%8.3f%8.3f\n",2,1,p1[0],p1[1],p1[2]);
 printf("ATOM  %5d  N   MOL %5d    %8.3f%8.3f%8.3f\n",3,1,p2[0],p2[1],p2[2]);
*/
/*
*/
         p3[0]=p1[0]+rch*uzt[1];
         p3[1]=p1[1]+rch*uzt[2];
         p3[2]=p1[2]+rch*uzt[3];
/*
         p3[0]=uzt[1];
         p3[1]=uzt[2];
         p3[2]=uzt[3];
*/

/*
 printf("ATOM  %5d  H1  MOL %5d    %8.3f%8.3f%8.3f\n",4,1,p3[0],p3[1],p3[2]);
*/

}

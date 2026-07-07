#include "head.h"

#define TOL     .000000001


/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
void euler(_REAL u32[3], _REAL u12[3],
           _REAL vx[3], _REAL vy[3], _REAL ang[3])
{
      _REAL v[3],vt[3],rz[3][3],rx1[3][3],rz2[3][3],mod,sum;
      _REAL al,be,ga,scp;
      _REAL phase();
      int k,i;

      _REAL tol = 1.e-15;
      _REAL zero = 0.0;
      _REAL pi = _PI;
      
      scp = u32[0]*u12[0]+u32[1]*u12[1]+u32[2]*u12[2];
      if (Nabs(scp)>=1.-tol && Nabs(scp)<=1.+tol) {
       al = 0.;
       ang[0] = al;
       be = acos(u32[2]);
       ang[1] = be;
       ga = 0.;
       ang[2] = ga;
       vx[0]=1.;
       vx[1]=0.;
       vx[2]=0.;
       vy[0]=0.;
       vy[1]=1.;
       vy[2]=0.;
       return;
      }

      be = acos(u32[2]);
      ang[1] = be;

      vx[0] = u32[1]*u12[2] - u32[2]*u12[1];
      vx[1] = u32[2]*u12[0] - u32[0]*u12[2];
      vx[2] = u32[0]*u12[1] - u32[1]*u12[0];
      mod = sqrt(vx[0]*vx[0]+vx[1]*vx[1]+vx[2]*vx[2]);
      for(k=0;k<3;k++) vx[k]=vx[k]/mod;
      vy[0] = u32[1]*vx[2] - u32[2]*vx[1];
      vy[1] = u32[2]*vx[0] - u32[0]*vx[2];
      vy[2] = u32[0]*vx[1] - u32[1]*vx[0];
      if ((u32[1]>=zero-tol && u32[1]<=zero+tol) &&
     +    (u32[0]>=zero-tol && u32[0]<=zero+tol)) {
       al = 0.0;
       ga=phase(vx[0],vx[1]);
       ang[0] = al;
       ang[2] = ga;
       return;
      }
      v[0] = -u32[1];
      v[1] =  u32[0];
      v[2] =  0.0;
      mod = sqrt(v[0]*v[0]+v[1]*v[1]);
      for(k=0;k<3;k++) v[k]=v[k]/mod;
      al=phase(v[0],v[1]);

      rz[0][0] = cos(al);
      rz[1][0] = -sin(al);
      rz[2][0] = 0.0;
      rz[0][1] = -rz[1][0];
      rz[1][1] =  rz[0][0];
      rz[2][1] =  rz[2][0];
      rz[0][2] =  rz[2][0];
      rz[1][2] =  rz[2][1];
      rz[2][2] = 1.0;

      rx1[0][0] = 1.0;
      rx1[1][0] = 0.0;
      rx1[2][0] = 0.0;
      rx1[0][1] = rx1[1][0];
      rx1[1][1] = cos(be);
      rx1[2][1] = -sin(be);
      rx1[0][2] = rx1[2][0];
      rx1[1][2] = -rx1[2][1];
      rx1[2][2] = rx1[1][1];

      for (i=0;i<3;i++) {
        sum=0.;
        for (k=0;k<3;k++) sum+=rz[i][k]*vx[k];
        vt[i]=sum;
      }
      for (i=0;i<3;i++) {
        sum=0.;
        for (k=0;k<3;k++) sum+=rx1[i][k]*vt[k];
        v[i]=sum;
      }

      ga=phase(v[0],v[1]);
      ang[0] = al;
      ang[2] = ga;

/*
      printf("al,be,ga in euler: %e %e %e\n",al,be,ga);
*/

      return;
}
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
# define TRUE 1
# define FALSE 0
struct quadrantstr
   {
   unsigned xp:1;
   unsigned yp:1;
   unsigned xn:1;
   unsigned yn:1;
   unsigned xz:1;
   unsigned yz:1;
   unsigned x1:1;
   unsigned y1:1;
   unsigned xm1:1;
   unsigned ym1:1;
   };

_REAL phase(_REAL x,_REAL y)
{
/*
      unsigned xp:1,yp:1,xn:1,yn:1,xz:1,yz:1,x1:1,y1:1,xm1:1,ym1:1;
*/
      struct quadrantstr quadrant;
      _REAL xy,fi=0.;
      _REAL pi = _PI;
      _REAL tol = 1.e-15;
      _REAL zero = 0.0;
      _REAL one = 1.0;
      _REAL onem = -1.0;
      
      quadrant.xp = FALSE;
      quadrant.yp = FALSE;
      quadrant.xn = FALSE;
      quadrant.yn = FALSE; 
      quadrant.xz = FALSE; 
      quadrant.yz = FALSE; 
      quadrant.x1 = FALSE; 
      quadrant.y1 = FALSE; 
      quadrant.xm1 = FALSE; 
      quadrant.ym1 = FALSE; 

      fi = 0.0;
      xy = sqrt(x*x + y*y);
      x = x/xy;
      y = y/xy;
      if (x>=zero-tol && x<=zero+tol) quadrant.xz=TRUE;
      if (x>=one-tol && x<=one+tol) quadrant.x1=TRUE;
      if (x>=onem-tol && x<=onem+tol) quadrant.xm1=TRUE;
      if (y>=zero-tol && y<=zero+tol) quadrant.yz=TRUE;
      if (y>=one-tol && y<=one+tol) quadrant.y1=TRUE;
      if (y>=onem-tol && y<=onem+tol) quadrant.ym1=TRUE;
      if (x>zero+tol) quadrant.xp=TRUE;
      if (x<zero-tol) quadrant.xn=TRUE;
      if (y>zero+tol) quadrant.yp=TRUE;
      if (y<zero-tol) quadrant.yn=TRUE;

      if (quadrant.xz && quadrant.yz) {
       fi = zero;
       return(fi);
      }
      if (quadrant.x1 && quadrant.yz) {
       fi = zero;
       return(fi);
      }
      if (quadrant.y1 && quadrant.xz) {
       fi = pi/2.0;
       return(fi);
      }
      if (quadrant.xm1 && quadrant.yz) {
       fi = pi;
       return(fi);
      }
      if (quadrant.ym1 && quadrant.xz) {
       fi = pi/2.0+pi;
       return(fi);
      }
      if (quadrant.xp && quadrant.yp) {
       fi = asin(y);
       return(fi);
      }
      if (quadrant.xn && quadrant.yp) {
       fi = acos(x);
       return(fi);
      }
      if (quadrant.xp && quadrant.yn) {
       fi = asin(y);
       fi = fi + 2.*pi;
       return(fi);
      }
      if (quadrant.xn && quadrant.yn) {
       y = -1.*y;
       x = -1.*x;
       fi = acos(x) + pi;
      }
      return(fi);
}

/*****************************************************************************
                                                                proquat()
******************************************************************************/

void proquat(_REAL a[4], _REAL b[4], _REAL c[4])
{
      c[0] = a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];
      c[1] = a[0]*b[1]+a[1]*b[0]-a[2]*b[3]+a[3]*b[2];
      c[2] = a[0]*b[2]+a[2]*b[0]-a[3]*b[1]+a[1]*b[3];
      c[3] = a[0]*b[3]+a[3]*b[0]-a[1]*b[2]+a[2]*b[1];
}

void r2q(_REAL *r, _REAL *q)
{
 _REAL x,y,z,xy,mod,cbh,sbh,cah,sah;

 x = r[0]; y = r[1]; z = r[2];
 mod = l_sqrt(x*x+y*y+z*z);
 if( Nabs(mod)<=ZERO ) {
  q[0]=1.;
  q[1]=0.;
  q[2]=0.;
  q[3]=0.;
  return;
 }
 x /= mod; y /= mod; z /= mod;
 xy = l_sqrt(x*x+y*y);
 if( Nabs(xy)<=ZERO ) {
  q[0]=1.;
  q[1]=0.;
  q[2]=0.;
  q[3]=0.;
  return;
 }
 x /= xy; y /= xy;
 cbh = l_sqrt(0.5*(1.+z));
 sbh = l_sqrt(0.5*(1.-z));
 cah = l_sqrt(0.5*(1.+x));
 sah = SIGN(1,y)*l_sqrt(0.5*(1.-x));
 q[0] = cbh*cah;
 q[2] = sbh*cah;
 q[1] =-sbh*sah;
 q[3] = cbh*sah;
}


void procomp(a,b,c)
_REAL a[2],b[2],c[2];
{
   c[0] = a[0]*b[0]-a[1]*b[1];
   c[1] = a[0]*b[1]+a[1]*b[0];
}


/************************************************************************
							DIAGESORT()
*************************************************************************/
/*
 *  diagEsort() - diagonalize 3x3 matrix & sort eigenvectors by eigenvalues
 */

int
diagEsort(_REAL *mat, _REAL *Emat, _REAL *Evec[], 
          _REAL *Eigenvalue, int isort)
{
	int	njrot;
	int 	i, j, k, i3;
	_REAL	eigenvector[9], *eA, v;

	if (!jacobi3(mat, Eigenvalue, eigenvector, &njrot)) {
		printf("convergence failed\n");
		return(0);
	}
/*
	printf("\t[convergence in %d cycles]\n",njrot);
	printf("diagonalized matrix:\n");
        for(i=0;i<3;i++) {
         for(j=0;j<3;j++) {
          printf("%e ",*mat++);
         }
         printf("\n");
        }
        mat-=9;
*/
/*
	parr(mat);
	printf("eigenvectors:\n");
	parr(eigenvector);
	printf("eigenvalues:\n %e  %e  %e\n",
		Eigenvalue[0], Eigenvalue[1], Eigenvalue[2]);
*/

	for (i=i3=0; i<3; i++, i3+=3)
		for (j=0; j<3; j++)
			Emat[i3+j] = eigenvector[j*3+i];
/*
	printf("eigenvectors flopped:\n");
	parr(Emat);
*/
	/*
	 *  set up vector pointers to Ematrix
	 */

	for (i=0; i<3; i++)
		Evec[i] = (_REAL *) &Emat[i*3];

        if (isort) {
	/*
	 *  sort eigenvectors by eigenvalue
	 */

	for (i=0; i<2; i++) {
		v = Eigenvalue[k=i];
		for (j=i+1; j<3; j++)
			if (Eigenvalue[j] > v)  
				v = Eigenvalue[k=j];
		if (k != i) {

			/*
			 *  swap eigenvalues & eigenvector pointers
			 */

			Eigenvalue[k] = Eigenvalue[i];
			Eigenvalue[i] = v;
			eA = Evec[i];
			Evec[i] = Evec[k];
			Evec[k] = eA;
		}
	}
/*
	printf("e sorted\n%e   %e   %e\n%e   %e   %e\n%e   %e   %e\n",
		Evec[0][0], Evec[0][1], Evec[0][2],
		Evec[1][0], Evec[1][1], Evec[1][2],
		Evec[2][0], Evec[2][1], Evec[2][2]);
*/
        }
	return(1);
}


/*
 *  jacobi3() - get jacobian of 3x3 matrix
 */

#define ROTATE(ARR,MAJ1,MIN1,MAJ2,MIN2) { \
		g = ARR[MAJ1 + MIN1]; \
		h = ARR[MAJ2 + MIN2]; \
		ARR[MAJ1 + MIN1] = g - s*(h+g*tau); \
		ARR[MAJ2 + MIN2] = h + s*(g-h*tau); }

int
jacobi3(_REAL *a,_REAL *d,_REAL *v,int *nrot)
/* n must be 3.  see b[3] and z[3] below */
{
	/* routine adapted from Numerical Recipes in C */
	int 	i, j, ip, iq, p3, j3;
	_REAL 	tresh, theta, tau, t, sm, s, h, g, c, b[3], z[3];

	/*
	 *  initialize
	 */

	for (ip=p3=0; ip<3; ip++,p3+=3) {	
		/* initialize the identity matrix */
		for (iq=0; iq<3; iq++) 
			v[p3 + iq] = 0.0;
		v[p3 + ip] = 1.0;
		/* initialize b and d to diagonal of a */
		b[ip] = d[ip] = a[p3 + ip];
		z[ip] = 0.0;
	}
	*nrot = 0;
	for (i=0; i<50; i++) {		/* 50 tries */

		/* 
		 *  sum off-diagonal elements 
		 */

		sm = 0.0;
		for (ip=p3=0; ip<2; ip++,p3+=3) {	
			for (iq=ip+1; iq<3; iq++)
				sm += Nabs(a[p3 + iq]);
		}
/*
	printf("sm %e\n", sm);
*/

		/* 
		 *  normal return, which relies on quadratic 
		 * 	convergence to machine underflow
		 */

		if (sm == 0.0) {
			return(1);
		}
		if (i < 3) 
			tresh = sm * 0.2 / 9.0;	  /* on 1st three sweeps... */
		else       
			tresh = 0.0;		  /* thereafter... */
		for (ip=p3=0; ip<2; ip++,p3+=3) {
			for (iq=ip+1; iq<3; iq++) {
				g = 100.0 * Nabs(a[p3 + iq]);

				/* 
				 *  after four sweeps, skip the rotation if
		   		 * 	the off-diagonal element is small 
				 */

				if ( i > 3  &&  Nabs(d[ip])+g == Nabs(d[ip])
				   && Nabs(d[iq])+g == Nabs(d[iq])) {
					a[p3 + iq] = 0.0;
				} else if (Nabs(a[p3 + iq]) > tresh) {
					h = d[iq]-d[ip];
					if (Nabs(h)+g==Nabs(h))
						t = a[p3 + iq] / h;
					else {
						theta = 0.5 * h / a[p3 + iq];
						t = 1.0 / (Nabs(theta)+
						  (_REAL)sqrt(1.0+theta*theta));
						if (theta < 0.0) 
							t = -t;
					}
					c = 1.0 / (_REAL)sqrt(1+t*t);
					s = t * c;
					tau = s / (1.0+c);
					h = t * a[p3 + iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[p3 + iq] = 0.0;
					for (j=j3=0; j<=ip-1; j++,j3+=3) 
						ROTATE(a,j3,ip,j3,iq)
					for (j=ip+1; j<=iq-1; j++) 
						ROTATE(a,p3,j,j*3,iq)
					for (j=iq+1; j<3; j++) 
						ROTATE(a,p3,j,iq*3,j)
/*
printf("t %e  c %e  s %e  tau %e\n", t, c, s, tau);
*/
					for (j3=0; j3<9; j3+=3) 
						ROTATE(v,j3,ip,j3,iq)
/*
parr(v);
*/
					++(*nrot);
				}
			}
		}
		for (ip=0; ip<3; ip++) {
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	printf("Too many iterations in routine JACOBI\n");
	return(0);
}	

#undef ROTATE

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void jacobi(_REAL **a,int n,_REAL d[],_REAL **v)
{
	int j,iq,ip,i,nrot=0;
	_REAL tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	_REAL *vector();
	void free_vector();

	b=vector(1,n);
	z=vector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += Nabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*Nabs(a[ip][iq]);
				if (i > 4 && Nabs(d[ip])+g == Nabs(d[ip])
					&& Nabs(d[iq])+g == Nabs(d[iq]))
					a[ip][iq]=0.0;
				else if (Nabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if (Nabs(h)+g == Nabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(Nabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	printf("Too many iterations in routine JACOBI");
        exit(1);
}

#undef ROTATE

void eigsrt(_REAL d[],_REAL **v,int n)
{
	int k,j,i;
	_REAL p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}

void qeuler(_REAL q[4],
   _REAL *cosbh,_REAL *cosa,_REAL *sina,_REAL *cosg,_REAL *sing)
/* extract cos and sin of euler (Rose) angles from a true quaternion */
{
   _REAL sinbh,q0,q1,q2,q3;

   q0 = q[0];
   q1 = q[1];
   q2 = q[2];
   q3 = q[3];
   *cosbh = l_sqrt(q0*q0+q3*q3);
   sinbh = l_sqrt(q1*q1+q2*q2);

/*
   if (*cosbh ==1. || *cosbh == 0.) {
*/
   if ((*cosbh >=1.-TOL && *cosbh <= 1.+TOL)) { 
    *cosbh = 1.;
    *cosa = 1.;
    *sina = 0.;
    *cosg = 1.;
    *sing = 0.;
    return;
   }
   if ((*cosbh >=0.-TOL && *cosbh <= 0.+TOL)) {
    *cosbh = 0.;
    *cosa = 1.;
    *sina = 0.;
    *cosg = 1.;
    *sing = 0.;
    return;
   }
   q0 /= *cosbh;
   q1 /= sinbh;
   q2 /= sinbh;
   q3 /= *cosbh;
   *cosa = q0*q2+q1*q3;
   *cosg = q0*q2-q1*q3;
   *sina = q2*q3-q0*q1;
   *sing = q2*q3+q0*q1;
   return;
}

_REAL determinant(_REAL **a)
/* compute the determinant of 3x3 matrix NR matrix */
// non so se funziona (20 gen 2003)
{
  _REAL det;
  det= a[1][1]*a[2][2]*a[3][3]
      +a[1][2]*a[2][3]*a[3][1]
      +a[1][3]*a[2][1]*a[3][2]
      -a[1][1]*a[3][2]*a[2][3]
      -a[2][1]*a[1][2]*a[3][3]
      -a[3][1]*a[2][2]*a[1][3];
   return(det);
} 

void sort(int n,_REAL ra[])
{
	int l,j,ir,i;
	_REAL rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1)
			rra=ra[--l];
		else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) ++j;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		ra[i]=rra;
	}
}

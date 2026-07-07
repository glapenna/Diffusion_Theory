#include "head.h"

#define MAXTRY 3

void eigvec_nr(_REAL **a, int n, _REAL *wr, _REAL *z)
{
 int i,j,k,ia,it,ib;
 _REAL **aa,**u,**v,*w,*wi,*tau,*b,*y,lambda,taut,taun,tol,wmax,thresh,
       *zt;
 void srt();
 double drand48();
 void balancc(_REAL **a,int n);

 _REAL b0,sum,dev,err,bmod;

 aa=matrix(1,n,1,n);
 for(i=1;i<=n;i++) {
  for(j=1;j<=n;j++) {
   aa[i][j]=a[i][j];
  }
 }

 wi=vector(1,n);
 balancc(a,n);
 elmhes(a,n);
 hqr(a,n,wr,wi);
 srt(wr,n);

 tau=vector(1,n);
 for(i=1,k=0;i<=n;i++) {
/*  printf("eigvec: lambda(%d)= %e\n",i,wr[i]); */
  if(wr[i]!=0.) {
   k++;
   tau[k]=wr[i];
  }
 }
 /* return; */

/* u=matrix(1,n,1,n); */
 v=matrix(1,n,1,n);
 w=vector(1,n);
 b=vector(1,n);
 y=vector(1,n);
 for(ia=1;ia<=n;ia++)
 {
 lambda=tau[ia];
 err=0.001;
 dev=(1.-2.*drand48())*err;
/* taut = lambda*(1.+dev); */
 taut = lambda+dev;
/* printf("ia= %d, lambda= %e, dev= %e\n",ia,lambda,dev); */
 err=0.001;
 b0=1.00;
 for(i=1,sum=0.;i<=n;i++) {
  dev=(1.-2.*drand48())*err;
  b[i] = b0*(1.+dev);
  sum+=b[i]*b[i];
 }
 sum = l_sqrt(sum);
/* printf("b0mod = %e\n",sum); */
 for(i=1;i<=n;i++) {
  b[i] /= sum;
 }

 it=0;
 while((Nabs(taut-lambda)>ZERO || it==0) && it<MAXTRY) {

 for(i=1;i<=n;i++) {
  for(j=1;j<=n;j++) {
   a[i][j]=aa[i][j];
  }
  a[i][i]=aa[i][i]-taut;
 }
/*
 for(i=1;i<=n;i++) {
  for(j=1;j<=n;j++) {
   printf("%d %d %e\n",i,j,a[i][j]);
  }
 }
*/
 svdcmp(a,n,n,w,v);

 tol=TOLSVD;
 wmax=0.;
 for(j=1;j<=n;j++) {
  if(w[j]>wmax) wmax=w[j];
 }
/* printf("wmax: %e\n",wmax); */
 thresh=tol*wmax;
 for(j=1;j<=n;j++) {
  if(w[j]<thresh) {
   printf("eigvec: it = %d :",it);
   printf("A-tI is singular because w(%d) = %e\n",j,w[j]);
   printf("tau=%e, lambda=%e, tau-lambda=%e\n",taut,lambda,taut-lambda);
/*
   exit(1);
*/
   w[j]=0.;
  }
 }
 for(j=1;j<=n;j++) {
  if(w[j]!=0.) {
   w[j]=1./w[j];
  }
 }

 for(ib=1;ib<=MAXTRY;ib++) {

 for(i=1,bmod=0.,thresh=0.;i<=n;i++) {
  for(j=1,wmax=0.;j<=n;j++) {
   for(k=1,sum=0.;k<=n;k++) {
    sum+=v[i][k]*a[j][k]*w[k];
   }
   wmax+=sum*b[j];
  }
  y[i]=wmax;
  thresh+=y[i]*b[i];
/*  printf("y(%d) = %e\n",i,y[i]); */
  bmod+=wmax*wmax;
 }

 bmod=1./l_sqrt(bmod);
/* printf("bmod(%d) = %e\n",it,bmod); */
 for(i=1,dev=0.,sum=0.;i<=n;i++) {
  y[i]*=bmod;
  dev+=(y[i]-b[i])*(y[i]-b[i]);
  sum+=y[i]*y[i];
  b[i]=y[i];
 }
 sum=l_sqrt(sum);
 dev=l_sqrt(dev);
/*
 printf("vmod(%d) = %e\n",ib,sum);
 printf("dev(%d) = %e\n",ib,dev);
*/

 }

 taun = taut+1./thresh;
 it++;
/* printf("it = %d: old tau = %e  new tau = %e\n",it,taut,taun); */
 taut = taun;

 }
/*
 printf("ia = %d: A*v, lambda*v:\n",ia);
 for(i=1;i<=n;i++) {
  for(j=1,sum=0.;j<=n;j++) {
   sum+=aa[i][j]*b[j];
  }
  printf("%d %e %e\n",i,sum,lambda*b[i]);
 }
*/
 zt=z+ia-1;
 for(i=1;i<=n;i++,zt+=n) {
  *zt = b[i];
 }
  
 }
 zt=z;
/* free_matrix(u,1,n,1,n); */
/*
 free_matrix(aa,1,n,1,n);
 free_matrix(v,1,n,1,n);
 free_vector(tau,1,n);
 free_vector(w,1,n);
 free_vector(wi,1,n);
 free_vector(b,1,n);
 free_vector(y,1,n);
*/

}

#define RADIX 2.0

void balancc(_REAL **a,int n)
{
	int last,j,i;
	_REAL s,r,g,f,c,sqrdx;

	sqrdx=RADIX*RADIX;
	last=0;
	while (last == 0) {
		last=1;
		for (i=1;i<=n;i++) {
			r=c=0.0;
			for (j=1;j<=n;j++)
				if (j != i) {
					c += Nabs(a[j][i]);
					r += Nabs(a[i][j]);
				}
			if (c && r) {
				g=r/RADIX;
				f=1.0;
				s=c+r;
				while (c<g) {
					f *= RADIX;
					c *= sqrdx;
				}
				g=r*RADIX;
				while (c>g) {
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c+r)/f < 0.95*s) {
					last=0;
					g=1.0/f;
					for (j=1;j<=n;j++) a[i][j] *= g;
					for (j=1;j<=n;j++) a[j][i] *= f;
				}
			}
		}
	}
}

#undef RADIX

#define SWAP(g,h) {y=(g);(g)=(h);(h)=y;}

void elmhes(_REAL **a,int n)
{
	int m,j,i;
	_REAL y,x;

	for (m=2;m<n;m++) {
		x=0.0;
		i=m;
		for (j=m;j<=n;j++) {
			if (Nabs(a[j][m-1]) > Nabs(x)) {
				x=a[j][m-1];
				i=j;
			}
		}
		if (i != m) {
			for (j=m-1;j<=n;j++) SWAP(a[i][j],a[m][j])
			for (j=1;j<=n;j++) SWAP(a[j][i],a[j][m])
		}
		if (x) {
			for (i=m+1;i<=n;i++) {
				if (y=a[i][m-1]) {
					y /= x;
					a[i][m-1]=y;
					for (j=m;j<=n;j++)
						a[i][j] -= y*a[m][j];
					for (j=1;j<=n;j++)
						a[j][m] += y*a[j][i];
				}
			}
		}
	}
}

void hqr(_REAL **a,int n,_REAL wr[],_REAL wi[])
{
	int nn,m,l,k,j,its,i,mmin;
	_REAL z,y,x,w,v,u,t,s,r,q,p,anorm;
	void nrerror();

	anorm=Nabs(a[1][1]);
	for (i=2;i<=n;i++)
		for (j=(i-1);j<=n;j++)
			anorm += Nabs(a[i][j]);
	nn=n;
	t=0.0;
	while (nn >= 1) {
		its=0;
		do {
			for (l=nn;l>=2;l--) {
				s=Nabs(a[l-1][l-1])+Nabs(a[l][l]);
				if (s == 0.0) s=anorm;
				if (Nabs(a[l][l-1]) + s == s) break;
			}
			x=a[nn][nn];
			if (l == nn) {
				wr[nn]=x+t;
				wi[nn--]=0.0;
			} else {
				y=a[nn-1][nn-1];
				w=a[nn][nn-1]*a[nn-1][nn];
				if (l == (nn-1)) {
					p=0.5*(y-x);
					q=p*p+w;
					z=sqrt(Nabs(q));
					x += t;
					if (q >= 0.0) {
						z=p+SIGN(z,p);
						wr[nn-1]=wr[nn]=x+z;
						if (z) wr[nn]=x-w/z;
						wi[nn-1]=wi[nn]=0.0;
					} else {
						wr[nn-1]=wr[nn]=x+p;
						wi[nn-1]= -(wi[nn]=z);
					}
					nn -= 2;
				} else {
					if (its == 30) nrerror("Too many iterations in HQR");
					if (its == 10 || its == 20) {
						t += x;
						for (i=1;i<=nn;i++) a[i][i] -= x;
						s=Nabs(a[nn][nn-1])+Nabs(a[nn-1][nn-2]);
						y=x=0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					for (m=(nn-2);m>=l;m--) {
						z=a[m][m];
						r=x-z;
						s=y-z;
						p=(r*s-w)/a[m+1][m]+a[m][m+1];
						q=a[m+1][m+1]-z-r-s;
						r=a[m+2][m+1];
						s=Nabs(p)+Nabs(q)+Nabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u=Nabs(a[m][m-1])*(Nabs(q)+Nabs(r));
						v=Nabs(p)*(Nabs(a[m-1][m-1])+Nabs(z)+Nabs(a[m+1][m+1]));
						if (u+v == v) break;
					}
					for (i=m+2;i<=nn;i++) {
						a[i][i-2]=0.0;
						if  (i != (m+2)) a[i][i-3]=0.0;
					}
					for (k=m;k<=nn-1;k++) {
						if (k != m) {
							p=a[k][k-1];
							q=a[k+1][k-1];
							r=0.0;
							if (k != (nn-1)) r=a[k+2][k-1];
							if (x=Nabs(p)+Nabs(q)+Nabs(r)) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if (s=SIGN(sqrt(p*p+q*q+r*r),p)) {
							if (k == m) {
								if (l != m)
								a[k][k-1] = -a[k][k-1];
							} else
								a[k][k-1] = -s*x;
							p += s;
							x=p/s;
							y=q/s;
							z=r/s;
							q /= p;
							r /= p;
							for (j=k;j<=nn;j++) {
								p=a[k][j]+q*a[k+1][j];
								if (k != (nn-1)) {
									p += r*a[k+2][j];
									a[k+2][j] -= p*z;
								}
								a[k+1][j] -= p*y;
								a[k][j] -= p*x;
							}
							mmin = nn<k+3 ? nn : k+3;
							for (i=l;i<=mmin;i++) {
								p=x*a[i][k]+y*a[i][k+1];
								if (k != (nn-1)) {
									p += z*a[i][k+2];
									a[i][k+2] -= p*r;
								}
								a[i][k+1] -= p*q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		} while (l < nn-1);
	}
}


void srt(_REAL d[],int n)
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
		}
	}
}

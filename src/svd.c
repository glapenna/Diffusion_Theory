#include "head.h"

int invert(_REAL **a, int m, int n, _REAL **b)
{
 _REAL **at,**v, *w, tol, wmax, thresh, sum;
 int i,j,k,rank;

/*
 for(i=1;i<=3;i++) {
  for(j=1;j<=3;j++) {
   printf("a(%d,%d)= %f\n",i,j,a[i][j]);
  }
 }
*/

/*
*/
 at=matrix(1,n,1,m);
 for(i=1;i<=n;i++) {
  for(j=1;j<=n;j++) {
   at[i][j]=a[i][j];
  }
 }

 v=matrix(1,n,1,m);
 w=vector(1,m);
 printf("entro in svdcmp...\n");

 svdcmp(at,m,n,w,v);
/* svdcmp(a,m,n,w,v); */

 tol=TOLSVD;
 wmax=0.;
 for(j=1;j<=m;j++) {
  if(w[j]>wmax) wmax=w[j];
 }
 printf("wmax: %e\n",wmax);
 thresh=tol*wmax;
 rank=m;
 for(j=1;j<=m;j++) {
  if(w[j]<thresh) {
   rank--;
   printf("invert: a is singular because w(%d) = %e\n",j,w[j]);
   w[j]=0.;
  }
 }

 for(j=1;j<=m;j++) {
  if(w[j]!=0.) {
   w[j]=1./w[j];
  }
 }
 for(i=1;i<=n;i++) {
  for(j=1;j<=m;j++) {
   for(k=1,sum=0.;k<=m;k++) {
    sum+=v[i][k]*at[j][k]*w[k];
/*    sum+=v[i][k]*a[j][k]*w[k]; */
   }
   b[i][j]=sum;
  }
 }
 free_vector(w,1,m);
 free_matrix(v,1,n,1,m);
/*
*/
 free_matrix(at,1,n,1,m);
 return(rank);

}

#define MAXSVD 30

static _REAL at,bt,ct;
#define PYTHAG(a,b) ((at=Nabs(a)) > (bt=Nabs(b)) ? \
(ct=bt/at,at*l_sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*l_sqrt(1.0+ct*ct)): 0.0))

static _REAL maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))

void svdcmp(_REAL **a,int m,int n,_REAL *w,_REAL **v)
{
	int flag,i,its,j,jj,k,l,nm;
	_REAL c,f,h,s,x,y,z;
	_REAL anorm=0.0,g=0.0,scale=0.0;
	_REAL *rv1;

/*
        printf("svdcmp: m,n: %d\n",m,n);
        printf("svdcmp: a:\n");
        for(i=1;i<=n;i++) {
         for(j=1;j<=n;j++) {
          printf("i,j,a: %d %d %e\n",i,j,a[i][j]);
         }
        }
*/

	if (m < n) nrerror("SVDCMP: You must augment A with extra zero rows");
	rv1=vector(1,n);
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += Nabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(l_sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
					}
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += Nabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(l_sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
						for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(Nabs(w[i])+Nabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else {
			for (j=i;j<=m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=MAXSVD;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if (Nabs(rv1[l])+anorm == anorm) {
					flag=0;
					break;
				}
				if (Nabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					if (Nabs(f)+anorm != anorm) {
						g=w[i];
						h=PYTHAG(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-f*h);
						for (j=1;j<=m;j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
						}
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == MAXSVD) nrerror("No convergence in MAXSVD SVDCMP iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=PYTHAG(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=PYTHAG(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=PYTHAG(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);

/*
         for(i=1;i<=m;i++) {
           printf("svdcmp: w: %d %e\n",i,w[i]);
         }

        { _REAL sum;
         for(i=1;i<=m;i++) {
          for(j=1;j<=n;j++) {
           for(k=1,sum=0.;k<=n;k++) {
            sum+=a[i][k]*w[k]*v[j][k];
           }
           printf("svdcmp: u,v,a: %d %d %e %e %e\n",i,j,a[i][j],v[i][j],sum);
          }
         }
        }
*/

}

#undef MAX
#undef PYTHAG


void svbksb(_REAL **u,_REAL w[],_REAL **v,int m,int n,_REAL b[],_REAL x[])
{
	int jj,j,i;
	_REAL s,*tmp;

	tmp=vector(1,n);
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	free_vector(tmp,1,n);
}

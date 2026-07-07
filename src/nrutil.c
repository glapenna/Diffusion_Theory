#include "head.h"

void nrerror(char error_text[])
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

int *ivector(int nl,int nh)
{
	int *v;

	v=(int *)calloc((unsigned) (nh-nl+1),sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

_REAL *vector(int nl,int nh)
{
	_REAL *v;

	v=(_REAL *)calloc((unsigned) (nh-nl+1),sizeof(_REAL));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

_REAL **matrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	_REAL **m;

	m=(_REAL **) calloc((unsigned) (nrh-nrl+1),sizeof(_REAL*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(_REAL *) calloc((unsigned) (nch-ncl+1),sizeof(_REAL));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

fcomplex **Cmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	fcomplex **m;

	m=(fcomplex **) calloc((unsigned) (nrh-nrl+1),sizeof(fcomplex*));
	if (!m) nrerror("allocation failure 1 in Cmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(fcomplex *) calloc((unsigned) (nch-ncl+1),sizeof(fcomplex));
		if (!m[i]) nrerror("allocation failure 2 in Cmatrix()");
		m[i] -= ncl;
	}
	return m;
}


rdiadstruct **rdiad_matrix(int nrl,int nrh,int ncl,int nch)
{
 int i;
 rdiadstruct **m;

 m=(rdiadstruct **) calloc((unsigned) (nrh-nrl+1),sizeof(rdiadstruct *));
 if (!m) nrerror("allocation failure 1 in rdiad_matrix()");
 m -= nrl;

 for(i=nrl;i<=nrh;i++) {
  m[i]=(rdiadstruct *) calloc((unsigned) (nch-ncl+1),sizeof(rdiadstruct));
  if (!m[i]) nrerror("allocation failure 2 in rdiad_matrix()");
  m[i] -= ncl;
 }
 return m;
}

rotstruct **rot_matrix(int nrl,int nrh,int ncl,int nch)
{
 int i;
 rotstruct **m;

 m=(rotstruct **) calloc((unsigned) (nrh-nrl+1),sizeof(rotstruct *));
 if (!m) nrerror("allocation failure 1 in rot_matrix()");
 m -= nrl;

 for(i=nrl;i<=nrh;i++) {
  m[i]=(rotstruct *) calloc((unsigned) (nch-ncl+1),sizeof(rotstruct));
  if (!m[i]) nrerror("allocation failure 2 in rot_matrix()");
  m[i] -= ncl;
 }
 return m;
}

rot1struct **rot1_matrix(int nrl,int nrh,int ncl,int nch)
{
 int i;
 rot1struct **m;

 m=(rot1struct **) calloc((unsigned) (nrh-nrl+1),sizeof(rot1struct *));
 if (!m) nrerror("allocation failure 1 in rot1_matrix()");
 m -= nrl;

 for(i=nrl;i<=nrh;i++) {
  m[i]=(rot1struct *) calloc((unsigned) (nch-ncl+1),sizeof(rot1struct));
  if (!m[i]) nrerror("allocation failure 2 in rot1_matrix()");
  m[i] -= ncl;
 }
 return m;
}

/* tmatrix non funziona */
_REAL **tmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i,shift;
	_REAL **m;

	m=(_REAL **) calloc((unsigned) (nrh-nrl+1),sizeof(_REAL*));
	if (!m) nrerror("allocation failure 1 in tmatrix()");
	m -= nrl;

        if((nrh-nrl)>(nch-ncl)) {
         nrerror("more rows than cols in tmatrix()");
        }

	for(i=nrl,shift=0;i<=nrh;i++,shift++) {
		m[i]=(_REAL *)
                  calloc((unsigned) (nch-ncl-shift+1),sizeof(_REAL));
		if (!m[i]) nrerror("allocation failure 2 in tmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(int nrl,int nrh,int ncl,int nch)
{
	int i,**m;

	m=(int **)calloc((unsigned) (nrh-nrl+1),sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)calloc((unsigned) (nch-ncl+1),sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}



_REAL **submatrix(_REAL **a,int oldrl,int oldrh,
                  int oldcl,int oldch,int newrl,int newcl)
{
	int i,j;
	_REAL **m;

	m=(_REAL **) calloc((unsigned) (oldrh-oldrl+1),sizeof(_REAL*));
	if (!m) nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	return m;
}



void free_vector(_REAL *v,int nl,int nh)
{
	free((char*) (v+nl));
}

void free_ivector(int *v,int nl,int nh)
{
	free((char*) (v+nl));
}

void free_matrix(_REAL **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


void free_submatrix(_REAL **b,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (b+nrl));
}


_REAL **convert_matrix(_REAL *a,int nrl,int nrh,int ncl,int nch)
{
	int i,j,nrow,ncol;
	_REAL **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	m = (_REAL **) calloc((unsigned) (nrow),sizeof(_REAL*));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	return m;
}



void free_convert_matrix(_REAL **b,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (b+nrl));
}

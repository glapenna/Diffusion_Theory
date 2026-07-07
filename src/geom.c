#include "head.h"

/*****************************************************************************
								DISTANCESQ()
******************************************************************************/

_REAL
distancesq(_REAL p1[3], _REAL p2[3])
{
        _REAL	x, y, z;

/*
	printf("distsq / p1: %f %f %f  p2: %f %f %f\n", 
	p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
*/
        x = p1[0] - p2[0];
        x = x * x;
 
        y = p1[1] - p2[1];
        y = y * y;
 
        z = p1[2] - p2[2];
        z = z * z;

        return (x + y + z);
}

/*****************************************************************************
								DISTANCE()
******************************************************************************/

_REAL
distance(_REAL p1[3], _REAL p2[3])
{
        _REAL  x, y, z;

/*
	printf("dist / p1: %f %f %f  p2: %f %f %f\n", 
	p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
*/
        x = (_REAL) (p1[0] - p2[0]);
        x = x * x;
 
        y = (_REAL) (p1[1] - p2[1]);
        y = y * y;
 
        z = (_REAL) ( p1[2] - p2[2]);
        z = z * z;
 
        return ((_REAL)l_sqrt (x + y + z));
}

/*****************************************************************************
							    BOXDISTANCESQ()
******************************************************************************/

_REAL
boxdistancesq( _REAL p1[3], _REAL p2[3], _REAL box[3])
{
        _REAL	x, y, z;

/*
	printf("boxdistsq / p1: %f %f %f  p2: %f %f %f  box: %f %f %f\n", 
	p1[0], p1[1], p1[2], p2[0], p2[1], box[0], box[1], box[2]);
*/
        x = p1[0] - p2[0];
	x = x - (Nabs(x) + box[0] / 2.0) * ((x < 0) ? -1 : 1);
        x = x * x;
 
        y = (p1[1] - p2[1]);
	y = y - (Nabs(y) + box[1] / 2.0) * ((y < 0) ? -1 : 1);
        y = y * y;
 
        z = ( p1[2] - p2[2]);
	z = z - (Nabs(z) + box[2] / 2.0) * ((z < 0) ? -1 : 1);
        z = z * z;
 
        return (x + y + z);
}


/*****************************************************************************
								BOXDISTANCE()
******************************************************************************/

_REAL
distancebox( _REAL p1[3], _REAL p2[3], _REAL box[3],_REAL boxm1[3])
{
        _REAL  a[3], sum;
        int k;

/*
	printf("boxdist / p1: %f %f %f  p2: %f %f %f  box: %f %f %f\n", 
	p1[0], p1[1], p1[2], p2[0], p2[1], box[0], box[1], box[2]);
*/
        for(k=0;k<3;k++) a[k] = p2[k]-p1[k];

        for (k=0,sum=0.;k<3;k++) {
         a[k] -= box[k]*(_REAL)(int)(a[k]*boxm1[k]+((a[k]<0.)?-1:1)*0.5);
         sum += a[k]*a[k];
        }

        return ((_REAL)l_sqrt (sum));
}


/*****************************************************************************
								VECTOR()
******************************************************************************/

dr(_REAL a[3], _REAL b[3], _REAL c[3])
{
	a[0] = b[0] - c[0];
	a[1] = b[1] - c[1];
	a[2] = b[2] - c[2];
}

/*****************************************************************************
								SIZE()
******************************************************************************/
/* TODO make inline */
_REAL
size(_REAL a[3])
{
	return(l_sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
}

/*****************************************************************************
								NORMALIZE()
******************************************************************************/
/* TODO make inline */
void
normalize(_REAL a[3])
{
	_REAL	b;

	b = size(a);
/*
	printf("norm %e  %e  %e  size %e\n", a[0],a[1],a[2],b);
*/
	a[0] /= b;
	a[1] /= b;
	a[2] /= b;
}

/*****************************************************************************
								DOT()
******************************************************************************/

/* TODO make inline */
_REAL
dot(_REAL a[3], _REAL b[3])
{
/*
	printf ("dot %e = %e * %e  +  %e * %e  +  %e * %e\n",
		a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
		a[0],b[0],a[1],b[1],a[2],b[2]);
*/
	return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

/*****************************************************************************
								CROSS()
******************************************************************************/

/*
 *  cross() - a = b x c
 */

cross(_REAL a[3], _REAL b[3], _REAL c[3])
{
	a[0] = b[1] * c[2] - b[2] * c[1];
	a[1] = b[2] * c[0] - b[0] * c[2];
	a[2] = b[0] * c[1] - b[1] * c[0];
}


/*****************************************************************************
								ANGLE()
******************************************************************************/

_REAL
angle(_REAL *p0, _REAL *p1, _REAL *p2)
{
	if (p2 ) {	/* TODO is this ok for vectors? what about p0=p1 etc */

		/* 
		 *  3 points 
		 */

		_REAL	a[3], b[3];

		a[0] = p0[0] - p1[0];
		a[1] = p0[1] - p1[1];
		a[2] = p0[2] - p1[2];

		b[0] = p2[0] - p1[0];
		b[1] = p2[1] - p1[1];
		b[2] = p2[2] - p1[2];

		return(l_acos(dot(a, b)/(size(a)*size(b))));
	}
	/*
	 *  else 2 axes (including derived from plane)
	 */

	return(l_acos(dot(p0, p1)/(size(p0)*size(p1))));
}

/*****************************************************************************/

_REAL
anglebox(_REAL *p0, _REAL *p1, _REAL *p2, _REAL *box,_REAL *boxm1)
{
        int k;

	if (p2 ) {	/* TODO is this ok for vectors? what about p0=p1 etc */

		/* 
		 *  3 points 
		 */

		_REAL	a[3], b[3];

		a[0] = p0[0] - p1[0];
		a[1] = p0[1] - p1[1];
		a[2] = p0[2] - p1[2];

		b[0] = p2[0] - p1[0];
		b[1] = p2[1] - p1[1];
		b[2] = p2[2] - p1[2];

#ifdef PBC
    { _REAL x2,y2,z2;
     x2=a[0]*boxm1[0];
     y2=a[1]*boxm1[1];
     z2=a[2]*boxm1[2];
     a[0]-=box[0]*rint(x2);
     a[1]-=box[1]*rint(y2);
     a[2]-=box[2]*rint(z2);
     x2=b[0]*boxm1[0];
     y2=b[1]*boxm1[1];
     z2=b[2]*boxm1[2];
     b[0]-=box[0]*rint(x2);
     b[1]-=box[1]*rint(y2);
     b[2]-=box[2]*rint(z2);
    }
#endif
		return((_REAL)l_acos((_REAL)(dot(a, b)/(size(a)*size(b)))));
	}
	/*
	 *  else 2 axes (including derived from plane)
	 */

	return((_REAL)l_acos((_REAL)(dot(p0, p1)/(size(p0)*size(p1)))));
}

/*****************************************************************************
								TORSION()
******************************************************************************/

_REAL
torsion(_REAL *p1, _REAL *p2, _REAL *p3, _REAL *p4)
{
	_REAL	a[3], b[3], c[3], A[3], B[3], C[3], cp, sign, phase;

	/* 
	 *  4 points:  1	2	3	4
	 *		  -a->     <-b-    -c->  
	 *			A	B
	 */

	a[0] = p1[0] - p2[0];
	a[1] = p1[1] - p2[1];
	a[2] = p1[2] - p2[2];

	b[0] = p3[0] - p2[0];
	b[1] = p3[1] - p2[1];
	b[2] = p3[2] - p2[2];

	c[0] = p4[0] - p3[0];
	c[1] = p4[1] - p3[1];
	c[2] = p4[2] - p3[2];

	A[0] = a[1] * b[2] - a[2] * b[1];
	A[1] = a[2] * b[0] - a[0] * b[2];
	A[2] = a[0] * b[1] - a[1] * b[0];

	B[0] = c[1] * b[2] - c[2] * b[1];
	B[1] = c[2] * b[0] - c[0] * b[2];
	B[2] = c[0] * b[1] - c[1] * b[0];

	cp = ( A[0] * B[0]  +  A[1] * B[1]  +  A[2] * B[2] ) 
		/ (size(A) * size(B));

	C[0] = A[1] * B[2] - A[2] * B[1];
	C[1] = A[2] * B[0] - A[0] * B[2];
	C[2] = A[0] * B[1] - A[1] * B[0];

	if ( b[0] * C[0]  +  b[1] * C[1]  +  b[2] * C[2]   > 0.0) {
		sign = 1.0;
		phase = 0.0;
	} else {
		sign = -1.0;
		phase = 2.*_PI;
	}

	return(sign*(_REAL)l_acos((_REAL)cp)+phase);
}

/*****************************************************************************/

_REAL
torsionbox(_REAL *p1, _REAL *p2, _REAL *p3, _REAL *p4,
           _REAL *box,_REAL *boxm1)
{
        int k;
	_REAL	a[3], b[3], c[3], A[3], B[3], C[3], cp,sp, sign, phase,phi;

	/* 
	 *  4 points:  1	2	3	4
	 *		  -a->     -b->    -c->  
	 *			A	B
	 */

	a[0] = p2[0] - p1[0];
	a[1] = p2[1] - p1[1];
	a[2] = p2[2] - p1[2];

	b[0] = p3[0] - p2[0];
	b[1] = p3[1] - p2[1];
	b[2] = p3[2] - p2[2];

	c[0] = p4[0] - p3[0];
	c[1] = p4[1] - p3[1];
	c[2] = p4[2] - p3[2];

    for (k=0;k<3;k++) {
     a[k] -= box[k]*(_REAL)(int)(a[k]*boxm1[k]+((a[k]<0.)?-1:1)*0.5);
     b[k] -= box[k]*(_REAL)(int)(b[k]*boxm1[k]+((b[k]<0.)?-1:1)*0.5);
     c[k] -= box[k]*(_REAL)(int)(c[k]*boxm1[k]+((c[k]<0.)?-1:1)*0.5);
    }

	A[0] = a[1] * b[2] - a[2] * b[1];
	A[1] = a[2] * b[0] - a[0] * b[2];
	A[2] = a[0] * b[1] - a[1] * b[0];

	B[0] = b[1] * c[2] - b[2] * c[1];
	B[1] = b[2] * c[0] - b[0] * c[2];
	B[2] = b[0] * c[1] - b[1] * c[0];

	cp = ( A[0] * B[0]  +  A[1] * B[1]  +  A[2] * B[2] ) 
		/ (size(A) * size(B));
        if(cp>1.) cp=1.;
        else if(cp<-1.) cp=-1.;

	C[0] = b[1] * A[2] - b[2] * A[1];
	C[1] = b[2] * A[0] - b[0] * A[2];
	C[2] = b[0] * A[1] - b[1] * A[0];

	sp = ( C[0] * B[0]  +  C[1] * B[1]  +  C[2] * B[2] )
		/ (size(C) * size(B));
        if(sp>1.) sp=1.;
        else if(sp<-1.) sp=-1.;
        phi = l_acos(cp);
        phi = copysign(phi,sp);

	if ( sp > 0.0 ) {
		phase = 0.0;
	} else {
		phase = 2.*_PI;
	}
        phi+=phase;
/*
        phi=sign*(l_acos(cp))+phase;
        phi=sign*(l_acos(cp))+phase;
        if(phi>2.*_PI) phi-=2.*_PI;
        if(phi<0.) phi+=2.*_PI;
*/
	return(phi);
}

/*****************************************************************************
                                                                VECTOR()
******************************************************************************/

dvector(_REAL a[3], _REAL b[3], _REAL c[3])
{
        a[0] = b[0] - c[0];
        a[1] = b[1] - c[1];
        a[2] = b[2] - c[2];
}
/**************************************************************************
                                                    IMAGES()
***************************************************************************/

void images(_REAL *dist, idstruct *serptr)
{
   _REAL sx,sy,sz,xs,ys,zs;
   _REAL **c=serptr->idunion.ser.set->cell_to_cart;
   _REAL **c_1=serptr->idunion.ser.set->cart_to_cell;

   sx=dist[0]*c_1[1][1]+dist[1]*c_1[1][2]+dist[2]*c_1[1][3]; 
   sy=dist[0]*c_1[2][1]+dist[1]*c_1[2][2]+dist[2]*c_1[2][3]; 
   sz=dist[0]*c_1[3][1]+dist[1]*c_1[3][2]+dist[2]*c_1[3][3]; 
     
   xs=sx-(int)(sx+((dist[0]<0.)?-1:1)*0.5);
   ys=sy-(int)(sy+((dist[1]<0.)?-1:1)*0.5);
   zs=sz-(int)(sz+((dist[2]<0.)?-1:1)*0.5);

   dist[0]=xs*c[1][1]+ys*c[1][2]+zs*c[1][3]; 
   dist[1]=xs*c[2][1]+ys*c[2][2]+zs*c[2][3]; 
   dist[2]=xs*c[3][1]+ys*c[3][2]+zs*c[3][3]; 
}

void images_to(_REAL *dist, _REAL *box, _REAL *boxm1)
{
   _REAL xs,ys,zs,corr;
   const _REAL r75=1.3333333333333333;

   xs=dist[0]*boxm1[0];
   ys=dist[1]*boxm1[1];
   zs=dist[2]*boxm1[2];
   xs=xs-rint(xs);
   ys=ys-rint(ys);
   zs=zs-rint(zs);
   corr = 0.5* (int)(r75*(Nabs(xs)+Nabs(ys)+Nabs(zs)));
   xs=xs-SIGN(corr,xs);
   ys=ys-SIGN(corr,ys);
   zs=zs-SIGN(corr,zs);
   dist[0]=xs*box[0];
   dist[1]=ys*box[1];
   dist[2]=zs*box[2];
}

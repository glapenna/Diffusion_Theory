#define _REAL double

#include <math.h>
#include "complex.h"

fcomplex Cadd(fcomplex a,fcomplex b)
{	fcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

fcomplex Csub(fcomplex a,fcomplex b)
{	fcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}

fcomplex Cmul(fcomplex a,fcomplex b)
{	fcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

fcomplex Cmul4(fcomplex a,fcomplex b,fcomplex c,fcomplex d)
{	fcomplex e,f,g;
        e=Cmul(a,b);
        f=Cmul(c,d);
        g=Cmul(e,f);
	return g;
}

fcomplex Cmul6(fcomplex a,fcomplex b,fcomplex c,
               fcomplex d,fcomplex e,fcomplex f)
{	fcomplex t1,t2,t3;
        t1=Cmul4(a,b,c,d);
        t2=Cmul(e,f);
        t3=Cmul(t1,t2);
	return t3;
}

fcomplex Cmul8(fcomplex a,fcomplex b,fcomplex c,fcomplex d,
               fcomplex e,fcomplex f,fcomplex g,fcomplex h)
{	fcomplex t1,t2,t3;
        t1=Cmul4(a,b,c,d);
        t2=Cmul4(e,f,g,h);
        t3=Cmul(t1,t2);
	return t3;
}

fcomplex Complex(_REAL re,_REAL im)
{	fcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

fcomplex Conjg(fcomplex z)
{	fcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

fcomplex Ccopy(fcomplex z)
{	fcomplex c;
	c.r=z.r;
	c.i=z.i;
	return c;
}

_REAL Cmod(fcomplex z)
{	_REAL a;
	a = sqrt(z.r*z.r+z.i*z.i);
	return a;
}

fcomplex Cdiv(fcomplex a,fcomplex b)
{	fcomplex c;
	_REAL r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

_REAL Cabs(fcomplex z)
{	_REAL x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

fcomplex Csqrt(fcomplex z)
{	fcomplex c;
	_REAL x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

fcomplex RCmul(_REAL x,fcomplex a)
{	fcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}

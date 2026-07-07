#include "head.h"

void D2mn(_REAL *qld,fcomplex **D)
{
   fcomplex d0,d1,d0cc,d1cc,Dt,
     dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8,dd9,dd10;
   _REAL sqrt_six,sqrt_two,msqrt_two;
   int n1,n2;
 
   sqrt_six=l_sqrt(6.);
   sqrt_two=l_sqrt(2.);
   msqrt_two = -sqrt_two;
   n1=3;
   n2=5;

/*  d0 = Complex(qdl[0],-qdl[3]); d1 = Complex(qdl[2],-qdl[1]); */
    d0 = Complex(qld[0],-qld[3]); d1 = Complex(qld[2],-qld[1]);
    d0cc = Conjg(d0); d1cc = Conjg(d1);
    dd1 = Cmul(d0,d0);
    dd2 = Cmul(d0,d1);
    dd3 = Cmul(d0,d1cc);
    dd4 = Cmul(d0,d0cc);
    dd5 = Cmul(d1,d1);
    dd6 = Cmul(d1,d1cc);
    dd7 = Cmul(d1,d0cc);
    dd8 = Cmul(d1cc,d1cc);
    dd9 = Cmul(d1cc,d0cc);
    dd10 = Cmul(d0cc,d0cc);

/* D^{(1)} */
/*
    D[1][1] = dd1;
    D[1][2] = RCmul(msqrt_two,dd3);
    D[1][3] = dd8;
    D[2][1] = RCmul(sqrt_two,dd2);
    D[2][2] = Csub(dd4,dd6);
    D[2][3] = RCmul(msqrt_two,dd9);
    D[3][1] = dd5;
    D[3][2] = RCmul(sqrt_two,dd7);
    D[3][3] = dd10;
*/

/* D^{(2)} */
    D[1][1] = Complex(0.,0.);
    Dt = Cmul(dd1,dd1);
    D[1][1] = Cadd(D[1][1],Dt);

    D[1][2] = Complex(0.,0.);
    Dt = Cmul(dd3,dd1);
    Dt = RCmul(-2.,Dt);
    D[1][2] = Cadd(D[1][2],Dt);

    D[1][3] = Complex(0.,0.);
    Dt = Cmul(dd8,dd1);
    Dt = RCmul(sqrt_six,Dt);
    D[1][3] = Cadd(D[1][3],Dt);

    D[1][4] = Complex(0.,0.);
    Dt = Cmul(dd8,dd3);
    Dt = RCmul(-2.0,Dt);
    D[1][4] = Cadd(D[1][4],Dt);

    D[1][5] = Complex(0.,0.);
    Dt = Cmul(dd8,dd8);
    D[1][5] = Cadd(D[1][5],Dt);

    D[2][1] = Complex(0.,0.);
    Dt = Cmul(dd2,dd1);
    Dt = RCmul(2.,Dt);
    D[2][1] = Cadd(D[2][1],Dt);

    D[2][2] = Complex(0.,0.);
    Dt = Cmul(dd6,dd1);
    Dt = RCmul(-3.,Dt);
    D[2][2] = Cadd(D[2][2],Dt);
    Dt = Cmul(dd4,dd1);
    D[2][2] = Cadd(D[2][2],Dt);

    D[2][3] = Complex(0.,0.);
    Dt = Cmul(dd8,dd2);
    Dt = RCmul(sqrt_six,Dt);
    D[2][3] = Cadd(D[2][3],Dt);
    Dt = Cmul(dd9,dd1);
    Dt = RCmul(-sqrt_six,Dt);
    D[2][3] = Cadd(D[2][3],Dt);

    D[2][4] = Complex(0.,0.);
    Dt = Cmul(dd8,dd6);
    Dt = RCmul(-1.,Dt);
    D[2][4] = Cadd(D[2][4],Dt);
    Dt = Cmul(dd9,dd3);
    Dt = RCmul(3.,Dt);
    D[2][4] = Cadd(D[2][4],Dt);

    D[2][5] = Complex(0.,0.);
    Dt = Cmul(dd9,dd8);
    Dt = RCmul(-2.,Dt);
    D[2][5] = Cadd(D[2][5],Dt);

    D[3][1] = Complex(0.,0.);
    Dt = Cmul(dd5,dd1);
    Dt = RCmul(sqrt_six,Dt);
    D[3][1] = Cadd(D[3][1],Dt);

    D[3][2] = Complex(0.,0.);
    Dt = Cmul(dd6,dd2);
    Dt = RCmul(-sqrt_six,Dt);
    D[3][2] = Cadd(D[3][2],Dt);
    Dt = Cmul(dd7,dd1);
    Dt = RCmul(sqrt_six,Dt);
    D[3][2] = Cadd(D[3][2],Dt);

    D[3][3] = Complex(0.,0.);
    Dt = Cmul(dd8,dd5);
    D[3][3] = Cadd(D[3][3],Dt);
    Dt = Cmul(dd9,dd2);
    Dt = RCmul(-4.,Dt);
    D[3][3] = Cadd(D[3][3],Dt);
    Dt = Cmul(dd10,dd1);
    D[3][3] = Cadd(D[3][3],Dt);

    D[3][4] = Complex(0.,0.);
    Dt = Cmul(dd9,dd6);
    Dt = RCmul(sqrt_six,Dt);
    D[3][4] = Cadd(D[3][4],Dt);
    Dt = Cmul(dd10,dd3);
    Dt = RCmul(-sqrt_six,Dt);
    D[3][4] = Cadd(D[3][4],Dt);

    D[3][5] = Complex(0.,0.);
    Dt = Cmul(dd10,dd8);
    Dt = RCmul(sqrt_six,Dt);
    D[3][5] = Cadd(D[3][5],Dt);

    D[4][1] = Complex(0.,0.);
    Dt = Cmul(dd5,dd2);
    Dt = RCmul(2.,Dt);
    D[4][1] = Cadd(D[4][1],Dt);

    D[4][2] = Complex(0.,0.);
    Dt = Cmul(dd6,dd5);
    Dt = RCmul(-1.,Dt);
    D[4][2] = Cadd(D[4][2],Dt);
    Dt = Cmul(dd7,dd2);
    Dt = RCmul(3.,Dt);
    D[4][2] = Cadd(D[4][2],Dt);

    D[4][3] = Complex(0.,0.);
    Dt = Cmul(dd9,dd5);
    Dt = RCmul(-sqrt_six,Dt);
    D[4][3] = Cadd(D[4][3],Dt);
    Dt = Cmul(dd10,dd2);
    Dt = RCmul(sqrt_six,Dt);
    D[4][3] = Cadd(D[4][3],Dt);

    D[4][4] = Complex(0.,0.);
    Dt = Cmul(dd10,dd6);
    Dt = RCmul(-3.,Dt);
    D[4][4] = Cadd(D[4][4],Dt);
    Dt = Cmul(dd10,dd4);
    D[4][4] = Cadd(D[4][4],Dt);

    D[4][5] = Complex(0.,0.);
    Dt = Cmul(dd10,dd9);
    Dt = RCmul(-2.,Dt);
    D[4][5] = Cadd(D[4][5],Dt);

    D[5][1] = Complex(0.,0.);
    Dt = Cmul(dd5,dd5);
    D[5][1] = Cadd(D[5][1],Dt);

    D[5][2] = Complex(0.,0.);
    Dt = Cmul(dd7,dd5);
    Dt = RCmul(2.,Dt);
    D[5][2] = Cadd(D[5][2],Dt);

    D[5][3] = Complex(0.,0.);
    Dt = Cmul(dd10,dd5);
    Dt = RCmul(sqrt_six,Dt);
    D[5][3] = Cadd(D[5][3],Dt);

    D[5][4] = Complex(0.,0.);
    Dt = Cmul(dd10,dd7);
    Dt = RCmul(2.,Dt);
    D[5][4] = Cadd(D[5][4],Dt);

    D[5][5] = Complex(0.,0.);
    Dt = Cmul(dd10,dd10);
    D[5][5] = Cadd(D[5][5],Dt);

}

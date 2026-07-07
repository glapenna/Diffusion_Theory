#include "head.h"
float images(_REAL x,_REAL y,_REAL z){
   idstruct *serptr;
   _REAL sx,sy,sz,xs,ys,zs;
   _REAL **c=serptr->idunion.ser.set->cell_to_cart;
   _REAL **c_1=serptr->idunion.ser.set->cell_to_cart;
   

   sx=x*c_1[1][1]+y*c_1[1][2]+z*c_1[1][3]; 
   sy=x*c_1[2][1]+y*c_1[2][2]+z*c_1[2][3]; 
   sz=x*c_1[3][1]+y*c_1[3][2]+z*c_1[3][3]; 
     
   xs=sx-(_REAL)(int)(sx+((sx<o.)?.-1:1)*0.5);
   ys=sy-(_REAL)(int)(sy+((sy<o.)?.-1:1)*0.5);
   zs=sz-(_REAL)(int)(sz+((sz<o.)?.-1:1)*0.5);

   x=xs*c[1][1]+ys*c[2][1]+zs*c[3][1]; 
   y=xs*c[1][2]+ys*c[2][2]+zs*c[3][2]; 
   z=xs*c[1][3]+ys*c[2][3]+zs*c[3][3]; 
 }

#include <limits.h>
#include <malloc.h>
#include <math.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <stdio.h>
#include <stdlib.h>

int get_bin(FILE *,const int *natoms,float *coord,float *box);
int put_bin(FILE *,const int *natoms,float *coord,float *box);

main(int argc,char *argv[])
{
 float t,prec,*coord,box[9];
 int i,j,k,natoms,uz,reverse,xdrid,step,xtcmagic,skip,
     freq,init,err;
 char inputname[512],outputname[512];
 FILE *input,*output;
 XDR *xdrs;

/* default */
 init=1;
 prec=1000.;
 reverse=1;
 xtcmagic=1995;
 freq=1;
 skip=0;

 if(argc!=7) {
  printf("Error: wrong number of arguments (%d).\n",argc);
  printf("Usage:\n");
  printf("bin2xtc bin_file xtc_file direction precision skip frequency\n");
  exit(1);
 }
 strcpy(inputname,argv[1]);
 strcpy(outputname,argv[2]);
 reverse=atoi(argv[3]);
 prec=atof(argv[4]);
 skip=atoi(argv[5]);
 freq=atoi(argv[6]);

 if(reverse==1)
 {

 input=fopen(inputname,"r");
 if(input==NULL)
 {
  printf("Error opening input binary file.\n");
  exit(1);
 }
 fread(&natoms,sizeof(int),1,input);
/* overwrite natoms!!! */
// natoms=76;
 rewind(input);
 coord = calloc(3*natoms,sizeof(float));

/* set default box */ 
 box[0]=10000.;
 box[1]=0.;
 box[2]=0.;
 box[3]=0.;
 box[4]=10000.;
 box[5]=0.;
 box[6]=0.;
 box[7]=0.;
 box[8]=10000.;

 xdrs = (XDR *)calloc(1,sizeof(XDR));
 xdrid = xdropen(xdrs,outputname,"w");

 step=0;
 while(step<skip)
 {
  get_bin(input,&natoms,coord,box);
  step++;
 }

 t=0.;
 step=0;
 while(get_bin(input,&natoms,coord,box))
 {
  if((step+1)%freq==0)
  {
   printf("t= %f\n",t);
// overwriting box
/*
   box[0]=10.;
   box[4]=10.;
   box[8]=10.;
*/
   printf("boxz= %f\n",box[8]);
   xdr_int(xdrs,&xtcmagic);
   xdr_int(xdrs,&natoms);
   xdr_int(xdrs,&step);
   xdr_float(xdrs,&t);
   for(k=0;k<9;k++) xdr_float(xdrs,&box[k]);
   xdr3dfcoord(xdrs,coord,&natoms,&prec);
  }
  t+=1.;
  step++;
 }
 xdrclose(xdrs);

 } 
 else
 {

 xdrs = (XDR *)calloc(1,sizeof(XDR));
 xdrid = xdropen(xdrs,inputname,"r");

 err=1;
 if(init == 1)
 {
  err=xdr_int(xdrs,&xtcmagic);
  if(err==0) exit(1);
  xdr_int(xdrs,&natoms);
  printf("err= %d\n",err);
  printf("magic= %d\n",xtcmagic);
  printf("natoms= %d\n",natoms);
  coord = calloc(3*natoms,sizeof(float));
 }

 output=fopen(outputname,"w");
 if(output==NULL)
 {
  printf("Error opening output binary file.\n");
  exit(1);
 }

 j=1;
 while(err!=0)
 {
  if(j%freq==0)
  {
   if(init != 1)
   {
    err=xdr_int(xdrs,&xtcmagic);
    if(err==0) exit(1);
    xdr_int(xdrs,&natoms);
    printf("err= %d\n",err);
    printf("magic= %d\n",xtcmagic);
    printf("natoms= %d\n",natoms);
   }
   init=0;
   xdr_int(xdrs,&step);
   xdr_float(xdrs,&t);
   fflush(stdout);
   printf("t= %f\n",t);
   for(k=0;k<9;k++) xdr_float(xdrs,&box[k]);
   printf("box= %f %f %f\n",box[0],box[4],box[8]);
   xdr3dfcoord(xdrs,coord,&natoms,&prec);
   printf("rf= %f %f %f\n",
    coord[3*natoms-3],coord[3*natoms-2],coord[3*natoms-1]);
   if(err==1 && j>skip) put_bin(output,&natoms,coord,box);
  }
  j++;
 }
 xdrclose(xdrs);

 }

}

int get_bin(FILE *input,const int *natoms,float *coord,float *box)
{
 int i,natoms_t,success,iread;
 float box_t[3];

 success=0;
 fread(&natoms_t,sizeof(int),1,input);
 natoms_t=*natoms;
 if(natoms_t!=*natoms)
 {
  return(success);
 }
 /*
 */
 iread=fread(&box_t,3*sizeof(float),1,input);
 if(iread==0)
 {
  return(success);
 }
 *box++=0.1*box_t[0];
 *box++=0.;
 *box++=0.;
 *box++=0.;
 *box++=0.1*box_t[1];
 *box++=0.;
 *box++=0.;
 *box++=0.;
 *box++=0.1*box_t[2];
 box-=9;
 iread=fread(coord,3*natoms_t*sizeof(float),1,input);
/*
*/
 for(i=0;i<3*natoms_t;i++)
 {
  *coord++ *= 0.1;
 }
 coord-=3*natoms_t;
 if(iread==0)
 {
  return(success);
 }
 success=1;
 return(success);
}

int put_bin(FILE *output,const int *natoms,float *coord,float *box)
{
 int i,natoms_t,success,iwrite;
 float box_t[3];

 success=0;
 natoms_t=*natoms;
 fwrite(&natoms_t,sizeof(int),1,output);
 box_t[0]=10.*box[0];
 box_t[1]=10.*box[4];
 box_t[2]=10.*box[8];
 iwrite=fwrite(box_t,3*sizeof(float),1,output);
 if(iwrite==0)
 {
  return(success);
 }

 for(i=0;i<3*natoms_t;i++)
 {
  *coord++ *= 10.;
 }
 coord-=3*natoms_t;
 iwrite=fwrite(coord,3*natoms_t*sizeof(float),1,output);
 if(iwrite==0)
 {
  return(success);
 }
 success=1;
 return(success);
}

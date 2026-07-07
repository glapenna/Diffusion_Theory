#include "head.h"
#include <time.h>

// extern char tok[];
// char tok[];
// void *get2(int,int);

extern idstruct *startid(),*getstruct();
rotstruct **rot_matrix();
rot1struct **rot1_matrix();

/************************************************************************/

void
init_mca_rl()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(MCA_RL, id);

	dptr->idunion.mca_rl.cases = 0;
	dptr->idunion.mca_rl.dump = 0;
	dptr->idunion.mca_rl.rigid = 0;
	dptr->idunion.mca_rl.blocks = 72;
	dptr->idunion.mca_rl.nt12 = 3;
	dptr->idunion.mca_rl.nt1234 = 9;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("MCA_RL: expected file name","");
	strcpy(dptr->idunion.mca_rl.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.mca_rl.fname, 0, &exists, "MCA_RL"))
		exit(1);
	if ((dptr->idunion.mca_rl.file = 
			 fopen(dptr->idunion.mca_rl.fname, "w")) == NULL){
		perror(dptr->idunion.mca_rl.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/
void
setup_mca_rl()
{
	idstruct	*dptr,*mptrr,*mptrl,*bptrr,*bptrl;
	int i,j,k,m1,m2,kk,blocki,blockj,skip;
        _REAL bavr,bavl;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != MCA_RL)
                        inerr("setup_mca_rl: expected id or keyword\n","");
        } else {
                printf("setup_mca_rl: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((mptrr = getstruct(tok)) == NULL) 
                inerr("setup_mca_rl: no such id\n","");
        if (mptrr->id_type != MCA_P2P1)
                inerr("setup_mca_rl: mca_p2p1 1 id expected\n","");
        dptr->idunion.mca_rl.mcar = mptrr;
        gettok();
        if ((mptrl = getstruct(tok)) == NULL) 
                inerr("setup_mca_rl: no such id\n","");
        if (mptrl->id_type != MCA_P2P1)
                inerr("setup_mca_rl: mca_p2p1 2 id expected\n","");
        dptr->idunion.mca_rl.mcal = mptrl;

        gettok();
        for(;;) {
         if(tok[0]==';') {
          break;
         } else if(!strcmp(tok,"READFS")) {
          dptr->idunion.mca_rl.dump=-1;
          gettok();
         } else if(!strcmp(tok,"RIGID")) {
          dptr->idunion.mca_rl.rigid=1;
          gettok();
         } else {
          printf("setup_mca_rl: option not found.\n");
          exit(1);
         }
        }

        dptr->idunion.mca_rl.type = 0;

        dptr->idunion.mca_rl.n0=mptrr->idunion.mca_p2p1.n0+
                                mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n12_rl=mptrr->idunion.mca_p2p1.n12*
                                    mptrl->idunion.mca_p2p1.n12;

        dptr->idunion.mca_rl.n1234_rrrl =mptrr->idunion.mca_p2p1.n12+
                                         dptr->idunion.mca_rl.n12_rl;
        dptr->idunion.mca_rl.n1234_rrll =mptrr->idunion.mca_p2p1.n12+
                                         mptrl->idunion.mca_p2p1.n12;
        dptr->idunion.mca_rl.n1234_rlrr =dptr->idunion.mca_rl.n12_rl+
                                         mptrr->idunion.mca_p2p1.n12;
        dptr->idunion.mca_rl.n1234_rlrl =dptr->idunion.mca_rl.n12_rl+
                                         dptr->idunion.mca_rl.n12_rl;
        dptr->idunion.mca_rl.n1234_rlll =dptr->idunion.mca_rl.n12_rl+
                                         mptrl->idunion.mca_p2p1.n12;
        dptr->idunion.mca_rl.n1234_llrr =mptrl->idunion.mca_p2p1.n12+
                                         mptrr->idunion.mca_p2p1.n12;
        dptr->idunion.mca_rl.n1234_llrl =mptrl->idunion.mca_p2p1.n12+
                                         dptr->idunion.mca_rl.n12_rl;

        blocki = dptr->idunion.mca_rl.nt12+dptr->idunion.mca_rl.nt1234;
        blockj = dptr->idunion.mca_rl.nt12+dptr->idunion.mca_rl.nt1234;
        dptr->idunion.mca_rl.dim= (int *)get2(dptr->idunion.mca_rl.blocks,sizeof(int));
        dptr->idunion.mca_rl.block= (int **)get2(blocki,sizeof(int *));
        for(i=0;i<blocki;i++) {
         dptr->idunion.mca_rl.block[i] = (int *)get2(blockj-i,sizeof(int));
         dptr->idunion.mca_rl.block[i]-=i;
        }
        dptr->idunion.mca_rl.diag= (int **)get2(blocki,sizeof(int *));
        for(i=0;i<blocki;i++) {
         dptr->idunion.mca_rl.diag[i] = (int *)get2(blockj-i,sizeof(int));
         dptr->idunion.mca_rl.diag[i]-=i;
        }

        for(i=0,k=0,skip=0;i<blocki;i++) {
         for(j=i;j<blockj;j++,k++) {
          if( (i==0 && j==0) ||
              (i==0 && j==3) ||
              (i==2 && j==2) ||
              (i==2 && j==11) ||
              (i==3 && j==3) ||
              (i==11 && j==11) ) {
           dptr->idunion.mca_rl.block[i][j] = -1;
           dptr->idunion.mca_rl.diag[i][j] = -1;
           skip++;
          } else {
           dptr->idunion.mca_rl.block[i][j] = k-skip;
           dptr->idunion.mca_rl.diag[i][j] = 0;
           if(i==j) {
            dptr->idunion.mca_rl.diag[i][j] = 1;
           }
          }
         }
        }


/* rr- ... */
        dptr->idunion.mca_rl.dim[0]=
         mptrr->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n12_rl;
        dptr->idunion.mca_rl.dim[1]=
         mptrr->idunion.mca_p2p1.n12*mptrl->idunion.mca_p2p1.n12;
        dptr->idunion.mca_rl.dim[2]=
         mptrr->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_rrrl;
        dptr->idunion.mca_rl.dim[3]=
         mptrr->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_rrll;
        dptr->idunion.mca_rl.dim[4]=
         mptrr->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_rlrr;
        dptr->idunion.mca_rl.dim[5]=
         mptrr->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_rlrl;
        dptr->idunion.mca_rl.dim[6]=
         mptrr->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_rlll;
        dptr->idunion.mca_rl.dim[7]=
         mptrr->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_llrr;
        dptr->idunion.mca_rl.dim[8]=
         mptrr->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_llrl;
        dptr->idunion.mca_rl.dim[9]=
         mptrr->idunion.mca_p2p1.n12*mptrl->idunion.mca_p2p1.n1234;

/* rl- ... */
        dptr->idunion.mca_rl.dim[10]=
         dptr->idunion.mca_rl.n12_rl*(dptr->idunion.mca_rl.n12_rl+1)/2;
        dptr->idunion.mca_rl.dim[11]=
         dptr->idunion.mca_rl.n12_rl*mptrl->idunion.mca_p2p1.n12;
        dptr->idunion.mca_rl.dim[12]=
         dptr->idunion.mca_rl.n12_rl*mptrr->idunion.mca_p2p1.n1234;
        dptr->idunion.mca_rl.dim[13]=
         dptr->idunion.mca_rl.n12_rl*dptr->idunion.mca_rl.n1234_rrrl;
        dptr->idunion.mca_rl.dim[14]=
         dptr->idunion.mca_rl.n12_rl*dptr->idunion.mca_rl.n1234_rrll;
        dptr->idunion.mca_rl.dim[15]=
         dptr->idunion.mca_rl.n12_rl*dptr->idunion.mca_rl.n1234_rlrr;
        dptr->idunion.mca_rl.dim[16]=
         dptr->idunion.mca_rl.n12_rl*dptr->idunion.mca_rl.n1234_rlrl;
        dptr->idunion.mca_rl.dim[17]=
         dptr->idunion.mca_rl.n12_rl*dptr->idunion.mca_rl.n1234_rlll;
        dptr->idunion.mca_rl.dim[18]=
         dptr->idunion.mca_rl.n12_rl*dptr->idunion.mca_rl.n1234_llrr;
        dptr->idunion.mca_rl.dim[19]=
         dptr->idunion.mca_rl.n12_rl*dptr->idunion.mca_rl.n1234_llrl;
        dptr->idunion.mca_rl.dim[20]=
         dptr->idunion.mca_rl.n12_rl*mptrl->idunion.mca_p2p1.n1234;

/* ll- ... */
        dptr->idunion.mca_rl.dim[21]=
         mptrl->idunion.mca_p2p1.n12*mptrr->idunion.mca_p2p1.n1234;
        dptr->idunion.mca_rl.dim[22]=
         mptrl->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_rrrl;
        dptr->idunion.mca_rl.dim[23]=
         mptrl->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_rrll;
        dptr->idunion.mca_rl.dim[24]=
         mptrl->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_rlrr;
        dptr->idunion.mca_rl.dim[25]=
         mptrl->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_rlrl;
        dptr->idunion.mca_rl.dim[26]=
         mptrl->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_rlll;
        dptr->idunion.mca_rl.dim[27]=
         mptrl->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_llrr;
        dptr->idunion.mca_rl.dim[28]=
         mptrl->idunion.mca_p2p1.n12*dptr->idunion.mca_rl.n1234_llrl;

/* rrrr- ... */
        dptr->idunion.mca_rl.dim[29]=
         mptrr->idunion.mca_p2p1.n1234*dptr->idunion.mca_rl.n1234_rrrl;
        dptr->idunion.mca_rl.dim[30]=
         mptrr->idunion.mca_p2p1.n1234*dptr->idunion.mca_rl.n1234_rrll;
        dptr->idunion.mca_rl.dim[31]=
         mptrr->idunion.mca_p2p1.n1234*dptr->idunion.mca_rl.n1234_rlrr;
        dptr->idunion.mca_rl.dim[32]=
         mptrr->idunion.mca_p2p1.n1234*dptr->idunion.mca_rl.n1234_rlrl;
        dptr->idunion.mca_rl.dim[33]=
         mptrr->idunion.mca_p2p1.n1234*dptr->idunion.mca_rl.n1234_rlll;
        dptr->idunion.mca_rl.dim[34]=
         mptrr->idunion.mca_p2p1.n1234*dptr->idunion.mca_rl.n1234_llrr;
        dptr->idunion.mca_rl.dim[35]=
         mptrr->idunion.mca_p2p1.n1234*dptr->idunion.mca_rl.n1234_llrl;
        dptr->idunion.mca_rl.dim[36]=
         mptrr->idunion.mca_p2p1.n1234*mptrl->idunion.mca_p2p1.n1234;

/* rrrl- ... */
        dptr->idunion.mca_rl.dim[37]=
         dptr->idunion.mca_rl.n1234_rrrl*(dptr->idunion.mca_rl.n1234_rrrl+1)/2;
        dptr->idunion.mca_rl.dim[38]=
         dptr->idunion.mca_rl.n1234_rrrl*dptr->idunion.mca_rl.n1234_rrll;
        dptr->idunion.mca_rl.dim[39]=
         dptr->idunion.mca_rl.n1234_rrrl*dptr->idunion.mca_rl.n1234_rlrr;
        dptr->idunion.mca_rl.dim[40]=
         dptr->idunion.mca_rl.n1234_rrrl*dptr->idunion.mca_rl.n1234_rlrl;
        dptr->idunion.mca_rl.dim[41]=
         dptr->idunion.mca_rl.n1234_rrrl*dptr->idunion.mca_rl.n1234_rlll;
        dptr->idunion.mca_rl.dim[42]=
         dptr->idunion.mca_rl.n1234_rrrl*dptr->idunion.mca_rl.n1234_llrr;
        dptr->idunion.mca_rl.dim[43]=
         dptr->idunion.mca_rl.n1234_rrrl*dptr->idunion.mca_rl.n1234_llrl;
        dptr->idunion.mca_rl.dim[44]=
         dptr->idunion.mca_rl.n1234_rrrl*mptrl->idunion.mca_p2p1.n1234;

/* rrll- ... */
        dptr->idunion.mca_rl.dim[45]=
         dptr->idunion.mca_rl.n1234_rrll*(dptr->idunion.mca_rl.n1234_rrll+1)/2;
        dptr->idunion.mca_rl.dim[46]=
         dptr->idunion.mca_rl.n1234_rrll*dptr->idunion.mca_rl.n1234_rlrr;
        dptr->idunion.mca_rl.dim[47]=
         dptr->idunion.mca_rl.n1234_rrll*dptr->idunion.mca_rl.n1234_rlrl;
        dptr->idunion.mca_rl.dim[48]=
         dptr->idunion.mca_rl.n1234_rrll*dptr->idunion.mca_rl.n1234_rlll;
        dptr->idunion.mca_rl.dim[49]=
         dptr->idunion.mca_rl.n1234_rrll*dptr->idunion.mca_rl.n1234_llrr;
        dptr->idunion.mca_rl.dim[50]=
         dptr->idunion.mca_rl.n1234_rrll*dptr->idunion.mca_rl.n1234_llrl;
        dptr->idunion.mca_rl.dim[51]=
         dptr->idunion.mca_rl.n1234_rrll*mptrl->idunion.mca_p2p1.n1234;

/* rlrr- ... */
        dptr->idunion.mca_rl.dim[52]=
         dptr->idunion.mca_rl.n1234_rlrr*(dptr->idunion.mca_rl.n1234_rlrr+1)/2;
        dptr->idunion.mca_rl.dim[53]=
         dptr->idunion.mca_rl.n1234_rlrr*dptr->idunion.mca_rl.n1234_rlrl;
        dptr->idunion.mca_rl.dim[54]=
         dptr->idunion.mca_rl.n1234_rlrr*dptr->idunion.mca_rl.n1234_rlll;
        dptr->idunion.mca_rl.dim[55]=
         dptr->idunion.mca_rl.n1234_rlrr*dptr->idunion.mca_rl.n1234_llrr;
        dptr->idunion.mca_rl.dim[56]=
         dptr->idunion.mca_rl.n1234_rlrr*dptr->idunion.mca_rl.n1234_llrl;
        dptr->idunion.mca_rl.dim[57]=
         dptr->idunion.mca_rl.n1234_rlrr*mptrl->idunion.mca_p2p1.n1234;

/* rlrl- ... */
        dptr->idunion.mca_rl.dim[58]=
         dptr->idunion.mca_rl.n1234_rlrl*(dptr->idunion.mca_rl.n1234_rlrl+1)/2;
        dptr->idunion.mca_rl.dim[59]=
         dptr->idunion.mca_rl.n1234_rlrl*dptr->idunion.mca_rl.n1234_rlll;
        dptr->idunion.mca_rl.dim[60]=
         dptr->idunion.mca_rl.n1234_rlrl*dptr->idunion.mca_rl.n1234_llrr;
        dptr->idunion.mca_rl.dim[61]=
         dptr->idunion.mca_rl.n1234_rlrl*dptr->idunion.mca_rl.n1234_llrl;
        dptr->idunion.mca_rl.dim[62]=
         dptr->idunion.mca_rl.n1234_rlrl*mptrl->idunion.mca_p2p1.n1234;

/* rlll- ... */
        dptr->idunion.mca_rl.dim[63]=
         dptr->idunion.mca_rl.n1234_rlll*(dptr->idunion.mca_rl.n1234_rlll+1)/2;
        dptr->idunion.mca_rl.dim[64]=
         dptr->idunion.mca_rl.n1234_rlll*dptr->idunion.mca_rl.n1234_llrr;
        dptr->idunion.mca_rl.dim[65]=
         dptr->idunion.mca_rl.n1234_rlll*dptr->idunion.mca_rl.n1234_llrl;
        dptr->idunion.mca_rl.dim[66]=
         dptr->idunion.mca_rl.n1234_rlll*mptrl->idunion.mca_p2p1.n1234;

/* llrr- ... */
        dptr->idunion.mca_rl.dim[67]=
         dptr->idunion.mca_rl.n1234_llrr*(dptr->idunion.mca_rl.n1234_llrr+1)/2;
        dptr->idunion.mca_rl.dim[68]=
         dptr->idunion.mca_rl.n1234_llrr*dptr->idunion.mca_rl.n1234_llrl;
        dptr->idunion.mca_rl.dim[69]=
         dptr->idunion.mca_rl.n1234_llrr*mptrl->idunion.mca_p2p1.n1234;

/* llrl- ... */
        dptr->idunion.mca_rl.dim[70]=
         dptr->idunion.mca_rl.n1234_llrl*(dptr->idunion.mca_rl.n1234_llrl+1)/2;
        dptr->idunion.mca_rl.dim[71]=
         dptr->idunion.mca_rl.n1234_llrl*mptrl->idunion.mca_p2p1.n1234;

        dptr->idunion.mca_rl.s = (_REAL **)
         get2(dptr->idunion.mca_rl.blocks,sizeof(_REAL *));
        dptr->idunion.mca_rl.f = (_REAL **)
         get2(dptr->idunion.mca_rl.blocks,sizeof(_REAL *));
        for(j=0;j<dptr->idunion.mca_rl.blocks;j++) {
         dptr->idunion.mca_rl.s[j] = (_REAL *)
          get2(dptr->idunion.mca_rl.dim[j],sizeof(_REAL));
         dptr->idunion.mca_rl.f[j] = (_REAL *)
          get2(dptr->idunion.mca_rl.dim[j],sizeof(_REAL));

         for(i=0;i<dptr->idunion.mca_rl.dim[j];i++) {
          dptr->idunion.mca_rl.s[j][i] =0.;
          dptr->idunion.mca_rl.f[j][i] =0.;
         }
        }


        dptr->idunion.mca_rl.ns = (int *)
         get2(dptr->idunion.mca_rl.nt12+dptr->idunion.mca_rl.nt1234,sizeof(int ));
        dptr->idunion.mca_rl.i12_2 = (int *)
         get2(dptr->idunion.mca_rl.nt12,sizeof(int ));
        dptr->idunion.mca_rl.i1234_2 = (int *)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(int ));
        dptr->idunion.mca_rl.i1234_4 = (int *)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(int ));
        dptr->idunion.mca_rl.n12 = (int *)
         get2(dptr->idunion.mca_rl.nt12,sizeof(int ));
        dptr->idunion.mca_rl.n1234 = (int *)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(int ));
        dptr->idunion.mca_rl.n12_1 = (int *)
         get2(dptr->idunion.mca_rl.nt12,sizeof(int ));
        dptr->idunion.mca_rl.n12_2 = (int *)
         get2(dptr->idunion.mca_rl.nt12,sizeof(int ));
        dptr->idunion.mca_rl.n1234_1 = (int *)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(int ));
        dptr->idunion.mca_rl.n1234_2 = (int *)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(int ));
        dptr->idunion.mca_rl.n1234_3 = (int *)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(int ));
        dptr->idunion.mca_rl.n1234_4 = (int *)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(int ));

        dptr->idunion.mca_rl.n12[0]=mptrr->idunion.mca_p2p1.n12;
        dptr->idunion.mca_rl.n12[1]=dptr->idunion.mca_rl.n12_rl;
        dptr->idunion.mca_rl.n12[2]=mptrl->idunion.mca_p2p1.n12;
        dptr->idunion.mca_rl.n12_1[0] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n12_2[0] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i12_2[0] = 1;
        dptr->idunion.mca_rl.n12_1[1] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n12_2[1] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i12_2[1] = 0;
        dptr->idunion.mca_rl.n12_1[2] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n12_2[2] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i12_2[2] = 1;

        dptr->idunion.mca_rl.n1234[0]=mptrr->idunion.mca_p2p1.n1234;
        dptr->idunion.mca_rl.n1234[1]=dptr->idunion.mca_rl.n1234_rrrl;
        dptr->idunion.mca_rl.n1234[2]=dptr->idunion.mca_rl.n1234_rrll;
        dptr->idunion.mca_rl.n1234[3]=dptr->idunion.mca_rl.n1234_rlrr;
        dptr->idunion.mca_rl.n1234[4]=dptr->idunion.mca_rl.n1234_rlrl;
        dptr->idunion.mca_rl.n1234[5]=dptr->idunion.mca_rl.n1234_rlll;
        dptr->idunion.mca_rl.n1234[6]=dptr->idunion.mca_rl.n1234_llrr;
        dptr->idunion.mca_rl.n1234[7]=dptr->idunion.mca_rl.n1234_llrl;
        dptr->idunion.mca_rl.n1234[8]=mptrl->idunion.mca_p2p1.n1234;
        dptr->idunion.mca_rl.n1234_1[0] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_2[0] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_3[0] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_4[0] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i1234_2[0] = 1;
        dptr->idunion.mca_rl.i1234_4[0] = 1;
        dptr->idunion.mca_rl.n1234_1[1] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_2[1] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_3[1] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_4[1] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i1234_2[1] = 1;
        dptr->idunion.mca_rl.i1234_4[1] = 0;
        dptr->idunion.mca_rl.n1234_1[2] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_2[2] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_3[2] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_4[2] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i1234_2[2] = 1;
        dptr->idunion.mca_rl.i1234_4[2] = 1;
        dptr->idunion.mca_rl.n1234_1[3] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_2[3] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_3[3] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_4[3] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i1234_2[3] = 0;
        dptr->idunion.mca_rl.i1234_4[3] = 1;
        dptr->idunion.mca_rl.n1234_1[4] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_2[4] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_3[4] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_4[4] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i1234_2[4] = 0;
        dptr->idunion.mca_rl.i1234_4[4] = 0;
        dptr->idunion.mca_rl.n1234_1[5] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_2[5] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_3[5] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_4[5] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i1234_2[5] = 0;
        dptr->idunion.mca_rl.i1234_4[5] = 1;
        dptr->idunion.mca_rl.n1234_1[6] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_2[6] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_3[6] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_4[6] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i1234_2[6] = 1;
        dptr->idunion.mca_rl.i1234_4[6] = 1;
        dptr->idunion.mca_rl.n1234_1[7] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_2[7] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_3[7] = mptrr->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_4[7] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i1234_2[7] = 1;
        dptr->idunion.mca_rl.i1234_4[7] = 0;
        dptr->idunion.mca_rl.n1234_1[8] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_2[8] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_3[8] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.n1234_4[8] = mptrl->idunion.mca_p2p1.n0;
        dptr->idunion.mca_rl.i1234_2[8] = 1;
        dptr->idunion.mca_rl.i1234_4[8] = 1;

        dptr->idunion.mca_rl.ns[0]=0;
        for(i=1;i<dptr->idunion.mca_rl.nt12;i++) {
         dptr->idunion.mca_rl.ns[i]=dptr->idunion.mca_rl.ns[i-1]
                                    +dptr->idunion.mca_rl.n12[i-1];
        }
        for(i=0,j=dptr->idunion.mca_rl.nt12;i<dptr->idunion.mca_rl.nt1234;i++,j++) {
         dptr->idunion.mca_rl.ns[j]=dptr->idunion.mca_rl.ns[j-1]
                                    +dptr->idunion.mca_rl.n1234[i-1];
        }

        dptr->idunion.mca_rl.l12 = (mod12struct **)
         get2(dptr->idunion.mca_rl.nt12,sizeof(mod12struct *));
        for(j=0;j<dptr->idunion.mca_rl.nt12;j++) {
         dptr->idunion.mca_rl.l12[j] = (mod12struct *)
          get2(dptr->idunion.mca_rl.n12[j],sizeof(mod12struct));
        }

        dptr->idunion.mca_rl.l1234 = (mod1234struct **)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(mod1234struct *));
        for(j=0;j<dptr->idunion.mca_rl.nt1234;j++) {
         dptr->idunion.mca_rl.l1234[j] = (mod1234struct *)
          get2(dptr->idunion.mca_rl.n1234[j],sizeof(mod1234struct));
        }

        dptr->idunion.mca_rl.p3 = (modstruct **)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(modstruct *));
        dptr->idunion.mca_rl.p4 = (modstruct **)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(modstruct *));
        dptr->idunion.mca_rl.m0 = (int *)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(int));
        dptr->idunion.mca_rl.alpha = (int *)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(int));
        dptr->idunion.mca_rl.beta = (int *)
         get2(dptr->idunion.mca_rl.nt1234,sizeof(int));
        j=1;
        dptr->idunion.mca_rl.p3[j]=mptrr->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.p4[j]=mptrl->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.m0[j]=dptr->idunion.mca_rl.ns[0];
        dptr->idunion.mca_rl.alpha[j]=1;
        dptr->idunion.mca_rl.beta[j]=0;
        j=2;
        dptr->idunion.mca_rl.p3[j]=mptrl->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.p4[j]=mptrl->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.m0[j]=dptr->idunion.mca_rl.ns[0];
        dptr->idunion.mca_rl.alpha[j]=1;
        dptr->idunion.mca_rl.beta[j]=0;
        j=3;
        dptr->idunion.mca_rl.p3[j]=mptrr->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.p4[j]=mptrl->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.m0[j]=dptr->idunion.mca_rl.ns[1];
        dptr->idunion.mca_rl.alpha[j]=0;
        dptr->idunion.mca_rl.beta[j]=1;
        j=4;
        dptr->idunion.mca_rl.p3[j]=mptrr->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.p4[j]=mptrl->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.m0[j]=dptr->idunion.mca_rl.ns[1];
        dptr->idunion.mca_rl.alpha[j]=0;
        dptr->idunion.mca_rl.beta[j]=1;
        j=5;
        dptr->idunion.mca_rl.p3[j]=mptrl->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.p4[j]=mptrl->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.m0[j]=dptr->idunion.mca_rl.ns[1];
        dptr->idunion.mca_rl.alpha[j]=0;
        dptr->idunion.mca_rl.beta[j]=1;
        j=6;
        dptr->idunion.mca_rl.p3[j]=mptrr->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.p4[j]=mptrr->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.m0[j]=dptr->idunion.mca_rl.ns[2];
        dptr->idunion.mca_rl.alpha[j]=1;
        dptr->idunion.mca_rl.beta[j]=0;
        j=7;
        dptr->idunion.mca_rl.p3[j]=mptrr->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.p4[j]=mptrl->idunion.mca_p2p1.m;
        dptr->idunion.mca_rl.m0[j]=dptr->idunion.mca_rl.ns[2];
        dptr->idunion.mca_rl.alpha[j]=1;
        dptr->idunion.mca_rl.beta[j]=0;

        dptr->idunion.mca_rl.Fm = rot_matrix
         (0,mptrr->idunion.mca_p2p1.n0-1,0,mptrl->idunion.mca_p2p1.n0-1);
        bptrr = mptrr->idunion.mca_p2p1.mca->idunion.mcarig.bead;
        bptrl = mptrl->idunion.mca_p2p1.mca->idunion.mcarig.bead;
        dptr->idunion.mca_rl.Fb = rot_matrix
         (0,bptrr->idunion.bead.nbond-1,0,bptrl->idunion.bead.nbond-1);

/*
        dptr->idunion.mca_rl.zav = 
          (bptrr->idunion.bead.nbond*bptrr->idunion.bead.zav +
           bptrl->idunion.bead.nbond*bptrl->idunion.bead.zav )/
          (bptrr->idunion.bead.nbond + bptrl->idunion.bead.nbond);
        bavr=bptrr->idunion.bead.bav*bptrr->idunion.bead.bav;
        bavl=bptrl->idunion.bead.bav*bptrl->idunion.bead.bav;
        dptr->idunion.mca_rl.bav = l_sqrt(
          (bptrr->idunion.bead.nbond*bavr + bptrl->idunion.bead.nbond*bavl )/
          (bptrr->idunion.bead.nbond + bptrl->idunion.bead.nbond) );

        bavr = bptrr->idunion.bead.zr*bptrr->idunion.bead.zav/bptrr->idunion.bead.bav;
        dptr->idunion.mca_rl.zr = bavr*
         dptr->idunion.mca_rl.bav/dptr->idunion.mca_rl.zav;
*/
        dptr->idunion.mca_rl.zr = bptrr->idunion.bead.zr;
        dptr->idunion.mca_rl.zav = bptrr->idunion.bead.zav;
        dptr->idunion.mca_rl.bav = bptrr->idunion.bead.bav;

        printf("setup_mca_rl: n0:      %d\n",dptr->idunion.mca_rl.n0);
        printf("setup_mca_rl: n12:     %d\n",dptr->idunion.mca_rl.n12);
        printf("setup_mca_rl: n1234:   %d\n",dptr->idunion.mca_rl.n1234);
        printf("setup_mca_rl: zr:      %f\n",dptr->idunion.mca_rl.zr);
        printf("setup_mca_rl: zav:     %f\n",dptr->idunion.mca_rl.zav);
        printf("setup_mca_rl: bav:     %f\n",dptr->idunion.mca_rl.bav);
}

/************************************************************************/
do_first_rl(idstruct *objptr,idstruct *serptr)
{
   _REAL *box=serptr->idunion.ser.set->boxc;
   idstruct *mcar,*mcal,*mptrr,*mptrl,*bptrr,*bptrl,*bptr;

   int nt12=objptr->idunion.mca_rl.nt12;
   int nt1234=objptr->idunion.mca_rl.nt1234;
   int *ns=objptr->idunion.mca_rl.ns;
   int *n12=objptr->idunion.mca_rl.n12;
   int *n1234=objptr->idunion.mca_rl.n1234;
   int *n12_1=objptr->idunion.mca_rl.n12_1;
   int *n12_2=objptr->idunion.mca_rl.n12_2;
   int *n1234_1=objptr->idunion.mca_rl.n1234_1;
   int *n1234_2=objptr->idunion.mca_rl.n1234_2;
   int *n1234_3=objptr->idunion.mca_rl.n1234_3;
   int *n1234_4=objptr->idunion.mca_rl.n1234_4;
   int *m0=objptr->idunion.mca_rl.m0;
   int *alpha=objptr->idunion.mca_rl.alpha;
   int *beta=objptr->idunion.mca_rl.beta;
   mod12struct **l12=objptr->idunion.mca_rl.l12;
   mod1234struct **l1234=objptr->idunion.mca_rl.l1234;

   int i1,i2,i3,i4,j1,j2,j3,j4,j1_0,j2_0,j3_0,j4_0,
       diag,block,i12_2,j12_2,i1234_2,j1234_2,i1234_4,j1234_4,
       i,j,k,l,m,ib,jb,ibt,jbt,kk,m1,m2,nbondmr,nbondml,nmoder,nmodel,overlap;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,r,zi,zj,zr,cutsq,
         boxm1[3],x1,y1,z1,x2,y2,z2,dx,dy,dz,scp;
   bondstruct *lbr,*lbl,*lb1,*lb2;
   modstruct *pr,*pl,*p1,*p2,*p3,*p4;
   mod12struct *l12_rr,*l12_ll;
   mod1234struct *l1234_rrrr,*l1234_llll;
   _REAL xm,ym,zm,f[5][5],Ct[3],csum,D[3][3],b,c;

   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }
   ++objptr->idunion.mca_rl.cases;

   mcar=objptr->idunion.mca_rl.mcar;
   mcal=objptr->idunion.mca_rl.mcal;
   mptrr=mcar->idunion.mca_p2p1.mca;
   mptrl=mcal->idunion.mca_p2p1.mca;
   bptrr=mptrr->idunion.mcarig.bead;
   bptrl=mptrl->idunion.mcarig.bead;
   bptr=bptrr;
   zr=objptr->idunion.mca_rl.zr;
   cutsq=bptr->idunion.bead.cutsq;
   nmoder=mcar->idunion.mca_p2p1.n0;
   nmodel=mcal->idunion.mca_p2p1.n0;
   nbondmr=bptrr->idunion.bead.nbond;
   nbondml=bptrl->idunion.bead.nbond;
   nt12 = objptr->idunion.mca_rl.nt12;
   nt1234 = objptr->idunion.mca_rl.nt1234;

   b=l_sqrt(3./8.);
   c=0.5;

   pr=mcar->idunion.mca_p2p1.m;
   pl=mcal->idunion.mca_p2p1.m;
   l12_rr=mcar->idunion.mca_p2p1.l12;
/* rr */
   j=0;
   for (i1=0,i=0;i1<n12_1[j];i1++) {
    for (i2=i1;i2<n12_2[j];i2++,i++) {
     for(k=0;k<5;k++) {
      l12[j][i].phi[k]=l12_rr->phi[k];
      for(kk=0;kk<3;kk++) {
       l12[j][i].dphidm1[k][kk]=l12_rr->dphidm1[k][kk];
       l12[j][i].dphidm2[k][kk]=l12_rr->dphidm2[k][kk];
      }
     }
     l12_rr++;
    }
   }

/* rl */
   j=1;
   for (i1=0,p1=pr,i=0;i1<n12_1[j];i1++,p1++) {
    for (i2=0,p2=pl;i2<n12_2[j];i2++,p2++,i++) {
     scp = p1->m[0]*p2->m[0]+p1->m[1]*p2->m[1]+p1->m[2]*p2->m[2];
     l12[j][i].phi[0]=c*(3.*p1->m[2]*p2->m[2] - scp);
     l12[j][i].dphidm1[0][0]=-c*p2->m[0];
     l12[j][i].dphidm1[0][1]=-c*p2->m[1];
     l12[j][i].dphidm1[0][2]= 2.*c*p2->m[2];
     l12[j][i].dphidm2[0][0]=-c*p1->m[0];
     l12[j][i].dphidm2[0][1]=-c*p1->m[1];
     l12[j][i].dphidm2[0][2]= 2.*c*p1->m[2];
     l12[j][i].phi[1]=b*(p1->m[0]*p2->m[0] - p1->m[1]*p2->m[1]);
     l12[j][i].dphidm1[1][0]= b*p2->m[0];
     l12[j][i].dphidm1[1][1]=-b*p2->m[1];
     l12[j][i].dphidm1[1][2]= 0.;
     l12[j][i].dphidm2[1][0]= b*p1->m[0];
     l12[j][i].dphidm2[1][1]=-b*p1->m[1];
     l12[j][i].dphidm2[1][2]= 0.;
     l12[j][i].phi[2]=b*(p1->m[0]*p2->m[1] + p1->m[1]*p2->m[0]);
     l12[j][i].dphidm1[2][0]= b*p2->m[1];
     l12[j][i].dphidm1[2][1]= b*p2->m[0];
     l12[j][i].dphidm1[2][2]= 0.;
     l12[j][i].dphidm2[2][0]= b*p1->m[1];
     l12[j][i].dphidm2[2][1]= b*p1->m[0];
     l12[j][i].dphidm2[2][2]= 0.;
     l12[j][i].phi[3]=b*(p1->m[0]*p2->m[2] + p1->m[2]*p2->m[0]);
     l12[j][i].dphidm1[3][0]= b*p2->m[2];
     l12[j][i].dphidm1[3][1]= 0.;
     l12[j][i].dphidm1[3][2]= b*p2->m[0];
     l12[j][i].dphidm2[3][0]= b*p1->m[2];
     l12[j][i].dphidm2[3][1]= 0.;
     l12[j][i].dphidm2[3][2]= b*p1->m[0];
     l12[j][i].phi[4]=b*(p1->m[1]*p2->m[2] + p1->m[2]*p2->m[1]);
     l12[j][i].dphidm1[4][0]= 0.;
     l12[j][i].dphidm1[4][1]= b*p2->m[2];
     l12[j][i].dphidm1[4][2]= b*p2->m[1];
     l12[j][i].dphidm2[4][0]= 0.;
     l12[j][i].dphidm2[4][1]= b*p1->m[2];
     l12[j][i].dphidm2[4][2]= b*p1->m[1];
    }
   }

/* ll */
   j=2;
   l12_ll=mcal->idunion.mca_p2p1.l12;
   for (i1=0,i=0;i1<n12_1[j];i1++) {
    for (i2=i1;i2<n12_2[j];i2++,i++) {
     for(k=0;k<5;k++) {
      l12[j][i].phi[k]=l12_ll->phi[k];
      for(kk=0;kk<3;kk++) {
       l12[j][i].dphidm1[k][kk]=l12_ll->dphidm1[k][kk];
       l12[j][i].dphidm2[k][kk]=l12_ll->dphidm2[k][kk];
      }
     }
     l12_ll++;
    }
   }

/* rrrr */
   l1234_rrrr=mcar->idunion.mca_p2p1.l1234;
   j=0;
   for (i1=0,i=0;i1<n1234_1[j];i1++) {
    for (i2=i1;i2<n1234_2[j];i2++) {
     for (i3=0;i3<n1234_3[j];i3++) {
      for (i4=i3;i4<n1234_4[j];i4++,i++) {
       for(k=0;k<5;k++) {
        l1234[j][i].phi[k]=l1234_rrrr->phi[k];
        for(kk=0;kk<3;kk++) {
         l1234[j][i].dphidm1[k][kk]=l1234_rrrr->dphidm1[k][kk];
         l1234[j][i].dphidm2[k][kk]=l1234_rrrr->dphidm2[k][kk];
         l1234[j][i].dphidm3[k][kk]=l1234_rrrr->dphidm3[k][kk];
         l1234[j][i].dphidm4[k][kk]=l1234_rrrr->dphidm4[k][kk];
        }
       }
       l1234_rrrr++;
      }
     }
    }
   }

/* mix */
   for(j=1;j<8;j++) {
    p3=objptr->idunion.mca_rl.p3[j];
    p4=objptr->idunion.mca_rl.p4[j];
    i1234_2=objptr->idunion.mca_rl.i1234_2[j];
    i1234_4=objptr->idunion.mca_rl.i1234_4[j];
    j12_2=j/3;

   for (i1=0,i=0;i1<n1234_1[j];i1++) {
    for (i2=i1234_2*i1;i2<n1234_2[j];i2++) {
     for (i3=0;i3<n1234_3[j];i3++) {
      for (i4=i1234_4*i3;i4<n1234_4[j];i4++,i++) {
       m=m0[j]+alpha[j]*((2*n1234_2[j]-i1+1)*i1/2+(i2-i1))+beta[j]*(n1234_2[j]*i1+i2);
       scp=p3[i3].m[0]*p4[i4].m[0]+p3[i3].m[1]*p4[i4].m[1]+p3[i3].m[2]*p4[i4].m[2];
       for(k=0;k<5;k++) {
        l1234[j][i].phi[k]=l12[j12_2][m].phi[k]*scp;
        for(kk=0;kk<3;kk++) {
         l1234[j][i].dphidm1[k][kk]=l12[j12_2][m].dphidm1[k][kk]*scp;
         l1234[j][i].dphidm2[k][kk]=l12[j12_2][m].dphidm2[k][kk]*scp;
         l1234[j][i].dphidm3[k][kk]=l12[j12_2][m].phi[k]*p4[i4].m[kk];
         l1234[j][i].dphidm4[k][kk]=l12[j12_2][m].phi[k]*p3[i3].m[kk];
        }
       }
      }
     }
    }
   }

   } /* enddo j */

/* llll */
   l1234_llll=mcal->idunion.mca_p2p1.l1234;
   j=8;
   for (i1=0,i=0;i1<n1234_1[j];i1++) {
    for (i2=i1;i2<n1234_2[j];i2++) {
     for (i3=0;i3<n1234_3[j];i3++) {
      for (i4=i3;i4<n1234_4[j];i4++,i++) {
       for(k=0;k<5;k++) {
        l1234[j][i].phi[k]=l1234_llll->phi[k];
        for(kk=0;kk<3;kk++) {
         l1234[j][i].dphidm1[k][kk]=l1234_llll->dphidm1[k][kk];
         l1234[j][i].dphidm2[k][kk]=l1234_llll->dphidm2[k][kk];
         l1234[j][i].dphidm3[k][kk]=l1234_llll->dphidm3[k][kk];
         l1234[j][i].dphidm4[k][kk]=l1234_llll->dphidm4[k][kk];
        }
       }
       l1234_llll++;
      }
     }
    }
   }

   lbr=bptrr->idunion.bead.l;
   lbl=bptrl->idunion.bead.l;
   for (ib=0,lb1=lbr;ib<nbondmr;ib++,lb1++) {
    for (jb=0,lb2=lbl;jb<nbondml;jb++,lb2++) {
#include "fbb.c"
     for(m1=0;m1<3;m1++) {
      for(m2=m1;m2<3;m2++) {
       f[m2][m1]=f[m1][m2];
      }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_rl.Fb[ib][jb].f[m1][m2]=f[m1][m2];
/*
       printf("Fb(%d,%d,%d,%d)= %f\n",
        ib,jb,m1,m2,objptr->idunion.mca_rl.Fb[ib][jb].f[m1][m2]);
*/
      }
     }
    }
   }

   for (i=0;i<nmoder;i++) {
    for (j=0;j<nmodel;j++) {
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }
     for (ib=0;ib<nbondmr;ib++) {
      for (jb=0;jb<nbondml;jb++) {
       for(m1=0;m1<3;m1++) {
        for(m2=0;m2<3;m2++) {
         f[m1][m2]+=
          mptrr->idunion.mcarig.evec[ib+1][i+1] *
          objptr->idunion.mca_rl.Fb[ib][jb].f[m1][m2] *
          mptrl->idunion.mcarig.evec[jb+1][j+1];
        }
       }
      }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_rl.Fm[i][j].f[m1][m2]=f[m1][m2];
/*
       printf("Fm(%d,%d,%d,%d)= %f\n",
        i,j,m1,m2,objptr->idunion.mca_rl.Fm[i][j].f[m1][m2]);
*/
      }
     }
    }
   }

   for(ib=0;ib<nt12;ib++) {
    for(jb=ib;jb<nt12;jb++) {
     block=objptr->idunion.mca_rl.block[ib][jb];

   if( block>=0 ) {

   i12_2=objptr->idunion.mca_rl.i12_2[ib];
   j12_2=objptr->idunion.mca_rl.i12_2[jb];
   diag = objptr->idunion.mca_rl.diag[ib][jb];


   for (i1=0,i=0,kk=0;i1<n12_1[ib];i1++) {
    for (i2=i12_2*i1;i2<n12_2[ib];i2++,i++,kk++) {

    j1_0 = i1*diag;
    j2_0 = i2*diag;

   for (j1=j1_0,j=0;j1<n12_1[jb];j1++,j2_0=j12_2*j1) {
    for (j2=j2_0;j2<n12_2[jb];j2++,j++) {

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
      {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l12[ib][i].dphidm1[m1][k] *
          objptr->idunion.mca_rl.Fm[i1][j1].f[k][l]  *
          l12[jb][j].dphidm1[m2][l];
         f[m1][m2]+=
          l12[ib][i].dphidm1[m1][k] *
          objptr->idunion.mca_rl.Fm[i1][j2].f[k][l]  *
          l12[jb][j].dphidm2[m2][l];
         f[m1][m2]+=
          l12[ib][i].dphidm2[m1][k] *
          objptr->idunion.mca_rl.Fm[i2][j1].f[k][l]  *
          l12[jb][j].dphidm1[m2][l];
         f[m1][m2]+=
          l12[ib][i].dphidm2[m1][k] *
          objptr->idunion.mca_rl.Fm[i2][j2].f[k][l]  *
          l12[jb][j].dphidm2[m2][l];
        }
       }
  
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     objptr->idunion.mca_rl.f[block][kk] += scp/((_REAL)5);

     scp=0.;
     scp+=l12_rr[i].phi[0]*l12_ll[j].phi[0];
     scp+=2.*l12_rr[i].phi[1]*l12_ll[j].phi[1];
     scp+=2.*l12_rr[i].phi[2]*l12_ll[i].phi[2];
     scp+=2.*l12_rr[i].phi[3]*l12_ll[i].phi[3];
     scp+=2.*l12_rr[i].phi[4]*l12_ll[i].phi[4];
     objptr->idunion.mca_rl.s[block][kk] += scp/((_REAL)5);

    }
   }
    }
   }
    } /* endif block */
    } /* enddo jb */
   } /* enddo ib */

   for(ib=0,ibt=0;ib<nt12;ib++,ibt++) {
    for(jb=0,jbt=nt12;jb<nt1234;jb++,jbt++) {

     block=objptr->idunion.mca_rl.block[ibt][jbt];

   if( block>=0 ) {

   i12_2=objptr->idunion.mca_rl.i12_2[ib];
   j1234_2=objptr->idunion.mca_rl.i1234_2[jb];
   j1234_4=objptr->idunion.mca_rl.i1234_4[jb];
    
   for (i1=0,i=0,kk=0;i1<n12_1[ib];i1++) {
    for (i2=i12_2*i1;i2<n12_2[ib];i2++,i++) {

   for (j1=0,j=0;j1<n1234_1[jb];j1++) {
    for (j2=j1234_2*j1;j2<n1234_2[jb];j2++) {
     for (j3=0;j3<n1234_3[jb];j3++) {
      for (j4=j1234_4*j3;j4<n1234_4[jb];j4++,j++,kk++) {

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l12[ib][i].dphidm1[m1][k] *
          objptr->idunion.mca_rl.Fm[i1][j1].f[k][l]  *
          l1234[jb][j].dphidm1[m2][l];
         f[m1][m2]+=
          l12[ib][i].dphidm1[m1][k] *
          objptr->idunion.mca_rl.Fm[i1][j2].f[k][l]  *
          l1234[jb][j].dphidm2[m2][l];
         f[m1][m2]+=
          l12[ib][i].dphidm1[m1][k] *
          objptr->idunion.mca_rl.Fm[i1][j3].f[k][l]  *
          l1234[jb][j].dphidm3[m2][l];
         f[m1][m2]+=
          l12[ib][i].dphidm1[m1][k] *
          objptr->idunion.mca_rl.Fm[i1][j4].f[k][l]  *
          l1234[jb][j].dphidm4[m2][l];
         f[m1][m2]+=
          l12[ib][i].dphidm2[m1][k] *
          objptr->idunion.mca_rl.Fm[i2][j1].f[k][l]  *
          l1234[jb][j].dphidm1[m2][l];
         f[m1][m2]+=
          l12[ib][i].dphidm2[m1][k] *
          objptr->idunion.mca_rl.Fm[i2][j2].f[k][l]  *
          l1234[jb][j].dphidm2[m2][l];
         f[m1][m2]+=
          l12[ib][i].dphidm2[m1][k] *
          objptr->idunion.mca_rl.Fm[i2][j3].f[k][l]  *
          l1234[jb][j].dphidm3[m2][l];
         f[m1][m2]+=
          l12[ib][i].dphidm2[m1][k] *
          objptr->idunion.mca_rl.Fm[i2][j4].f[k][l]  *
          l1234[jb][j].dphidm4[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     objptr->idunion.mca_rl.f[block][kk] += scp/((_REAL)5);

     scp=0.;
     scp+=l12[ib][i].phi[0]*l1234[jb][j].phi[0];
     scp+=2.*l12[ib][i].phi[1]*l1234[jb][j].phi[1];
     scp+=2.*l12[ib][i].phi[2]*l1234[jb][j].phi[2];
     scp+=2.*l12[ib][i].phi[3]*l1234[jb][j].phi[3];
     scp+=2.*l12[ib][i].phi[4]*l1234[jb][j].phi[4];
     objptr->idunion.mca_rl.s[block][kk] += scp/((_REAL)5);

      }
     }
    }
   }
    }
   }

    } /* endif block */
    } /* enddo jb */
   } /* enddo ib */

   for(ib=0,ibt=nt12;ib<nt1234;ib++,ibt++) {
    for(jb=ib,jbt=ibt;jb<nt1234;jb++,jbt++) {

     block=objptr->idunion.mca_rl.block[ibt][jbt];

   if( block>=0 ) {

   i1234_2=objptr->idunion.mca_rl.i1234_2[ib];
   i1234_4=objptr->idunion.mca_rl.i1234_4[ib];
   j1234_2=objptr->idunion.mca_rl.i1234_2[jb];
   j1234_4=objptr->idunion.mca_rl.i1234_4[jb];
   diag=objptr->idunion.mca_rl.diag[ibt][jbt];

   for (i1=0,i=0,kk=0;i1<n1234_1[ib];i1++) {
    for (i2=i1234_2*i1;i2<n1234_2[ib];i2++) {
     for (i3=0;i3<n1234_3[ib];i3++) {
      for (i4=i1234_4*i3;i4<n1234_4[ib];i4++,i++) {
       j1_0=diag*i1;
       j2_0=diag*i2;
       j3_0=diag*i3;
       j4_0=diag*i4;

   for (j1=j1_0,j=0;j1<n1234_1[jb];j1++,j2_0=j1234_2*j1) {
    for (j2=j2_0;j2<n1234_2[jb];j2++) {
     for (j3=j3_0;j3<n1234_2[jb];j3++,j4_0=j1234_4*j3) {
      for (j4=j4_0;j4<n1234_4[jb];j4++,j++,kk++) {

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l1234[ib][i].dphidm1[m1][k] *
          objptr->idunion.mca_rl.Fm[i1][j1].f[k][l]  *
          l1234[jb][j].dphidm1[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm1[m1][k] *
          objptr->idunion.mca_rl.Fm[i1][j2].f[k][l]  *
          l1234[jb][j].dphidm2[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm1[m1][k] *
          objptr->idunion.mca_rl.Fm[i1][j3].f[k][l]  *
          l1234[jb][j].dphidm3[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm1[m1][k] *
          objptr->idunion.mca_rl.Fm[i1][j4].f[k][l]  *
          l1234[jb][j].dphidm4[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm2[m1][k] *
          objptr->idunion.mca_rl.Fm[i2][j1].f[k][l]  *
          l1234[jb][j].dphidm1[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm2[m1][k] *
          objptr->idunion.mca_rl.Fm[i2][j2].f[k][l]  *
          l1234[jb][j].dphidm2[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm2[m1][k] *
          objptr->idunion.mca_rl.Fm[i2][j3].f[k][l]  *
          l1234[jb][j].dphidm3[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm2[m1][k] *
          objptr->idunion.mca_rl.Fm[i2][j4].f[k][l]  *
          l1234[jb][j].dphidm4[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm3[m1][k] *
          objptr->idunion.mca_rl.Fm[i3][j1].f[k][l]  *
          l1234[jb][j].dphidm1[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm3[m1][k] *
          objptr->idunion.mca_rl.Fm[i3][j2].f[k][l]  *
          l1234[jb][j].dphidm2[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm3[m1][k] *
          objptr->idunion.mca_rl.Fm[i3][j3].f[k][l]  *
          l1234[jb][j].dphidm3[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm3[m1][k] *
          objptr->idunion.mca_rl.Fm[i3][j4].f[k][l]  *
          l1234[jb][j].dphidm4[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm4[m1][k] *
          objptr->idunion.mca_rl.Fm[i4][j1].f[k][l]  *
          l1234[jb][j].dphidm1[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm4[m1][k] *
          objptr->idunion.mca_rl.Fm[i4][j2].f[k][l]  *
          l1234[jb][j].dphidm2[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm4[m1][k] *
          objptr->idunion.mca_rl.Fm[i4][j3].f[k][l]  *
          l1234[jb][j].dphidm3[m2][l];
         f[m1][m2]+=
          l1234[ib][i].dphidm4[m1][k] *
          objptr->idunion.mca_rl.Fm[i4][j4].f[k][l]  *
          l1234[jb][j].dphidm4[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     objptr->idunion.mca_rl.f[block][kk] += scp/((_REAL)5);

     scp=0.;
     scp+=l1234[ib][i].phi[0]*l1234[jb][j].phi[0];
     scp+=2.*l1234[ib][i].phi[1]*l1234[jb][j].phi[1];
     scp+=2.*l1234[ib][i].phi[2]*l1234[jb][j].phi[2];
     scp+=2.*l1234[ib][i].phi[3]*l1234[jb][j].phi[3];
     scp+=2.*l1234[ib][i].phi[4]*l1234[jb][j].phi[4];
     objptr->idunion.mca_rl.s[block][kk] += scp/((_REAL)5);

      }
     }
    }
   }
      }
     }
    }
   }

    } /* endif block */
    } /* enddo jb */
   } /* enddo ib */

}

/************************************************************************/
void
write_first_rl(idstruct *dptr)
{
	int	i,j,i0,j0,i2,j2,kk,k,l,m,ib,jb,ibt,jbt,block,diag;
	int	states = dptr->idunion.mca_rl.cases;
	FILE 	*file = dptr->idunion.mca_rl.file,*dump;
        idstruct *mcar,*mcal,*mptrr,*mptrl,*bptrr,*bptrl;

        int nt12=dptr->idunion.mca_rl.nt12;
        int nt1234=dptr->idunion.mca_rl.nt1234;
        int *ns=dptr->idunion.mca_rl.ns;
        int *n12=dptr->idunion.mca_rl.n12;
        int *n1234=dptr->idunion.mca_rl.n1234;
	_REAL fact,tol,factns,sum,sum2;
        int rank;

        int nn,nf;
        _REAL **S;
        _REAL **F,**Fm1,**tmp,csum,ctmp;
        _REAL *ui,*v,*w;
#ifdef IMSL
        d_complex *wi,*vi;
#else
        _REAL *vi;
#endif

        mcar=dptr->idunion.mca_rl.mcar;
        mcal=dptr->idunion.mca_rl.mcal;
        mptrr=mcar->idunion.mca_p2p1.mca;
        bptrr=mptrr->idunion.mcarig.bead;
	fact = 1./ ((_REAL)states);
        factns = bptrr->idunion.bead.factns;

        printf("write_first_rl: number of states: %d\n",states);

        nn =   mcar->idunion.mca_p2p1.n12
              +dptr->idunion.mca_rl.n12_rl
              +mcal->idunion.mca_p2p1.n12
              +mcar->idunion.mca_p2p1.n1234
              +dptr->idunion.mca_rl.n1234_rrrl
              +dptr->idunion.mca_rl.n1234_rrll
              +dptr->idunion.mca_rl.n1234_rlrr
              +dptr->idunion.mca_rl.n1234_rlrl
              +dptr->idunion.mca_rl.n1234_rlll
              +dptr->idunion.mca_rl.n1234_llrr
              +dptr->idunion.mca_rl.n1234_llrl
              +mcal->idunion.mca_p2p1.n1234;
        F = matrix(1,nn,1,nn);
        S = matrix(1,nn,1,nn);

           ib=0; jb=0;
           ibt=0; jbt=0;
           i0=ns[ibt];
           j0=ns[jbt];
           for (i=1,kk=0; i<=n12[ib]; i++) {
            for (j=i; j<=n12[jb]; j++,kk++) {
             i2=i+i0;
             j2=j+j0;
             F[i2][j2]=mcar->idunion.mca_p2p1.f_I_I[kk];
             F[j2][i2]=F[i2][j2];
             S[i2][j2]=mcar->idunion.mca_p2p1.s_I_I[kk];
             S[j2][i2]=S[i2][j2];
            }
           }
           ib=0; jb=0;
           ibt=0; jbt=3;
           i0=ns[ibt];
           j0=ns[jbt];
           for (i=1,kk=0; i<=n12[ib]; i++) {
            for (j=1; j<=n1234[jb]; j++,kk++) {
             i2=i+i0;
             j2=j+j0;
             F[i2][j2]=mcar->idunion.mca_p2p1.f_I_II[kk];
             F[j2][i2]=F[i2][j2];
             S[i2][j2]=mcar->idunion.mca_p2p1.s_I_II[kk];
             S[j2][i2]=S[i2][j2];
            }
           }
           ib=2; jb=2;
           ibt=2; jbt=2;
           i0=ns[ibt];
           j0=ns[jbt];
           for (i=1,kk=0; i<=n12[ib]; i++) {
            for (j=i; j<=n12[jb]; j++,kk++) {
             i2=i+i0;
             j2=j+j0;
             F[i2][j2]=mcal->idunion.mca_p2p1.f_I_I[kk];
             F[j2][i2]=F[i2][j2];
             S[i2][j2]=mcal->idunion.mca_p2p1.s_I_I[kk];
             S[j2][i2]=S[i2][j2];
            }
           }
           ib=2; jb=8;
           ibt=2; jbt=11;
           i0=ns[ibt];
           j0=ns[jbt];
           for (i=1,kk=0; i<=n12[ib]; i++) {
            for (j=1; j<=n1234[jb]; j++,kk++) {
             i2=i+i0;
             j2=j+j0;
             F[i2][j2]=mcal->idunion.mca_p2p1.f_I_II[kk];
             F[j2][i2]=F[i2][j2];
             S[i2][j2]=mcal->idunion.mca_p2p1.s_I_II[kk];
             S[j2][i2]=S[i2][j2];
            }
           }
           ib=0; jb=0;
           ibt=3; jbt=3;
           i0=ns[ibt];
           j0=ns[jbt];
           for (i=1,kk=0; i<=n1234[ib]; i++) {
            for (j=i; j<=n1234[jb]; j++,kk++) {
             i2=i+i0;
             j2=j+j0;
             F[i2][j2]=mcar->idunion.mca_p2p1.f_II_II[kk];
             F[j2][i2]=F[i2][j2];
             S[i2][j2]=mcar->idunion.mca_p2p1.s_II_II[kk];
             S[j2][i2]=S[i2][j2];
            }
           }
           ib=8; jb=8;
           ibt=11; jbt=11;
           i0=ns[ibt];
           j0=ns[jbt];
           for (i=1,kk=0; i<=n1234[ib]; i++) {
            for (j=i; j<=n1234[jb]; j++,kk++) {
             i2=i+i0;
             j2=j+j0;
             F[i2][j2]=mcal->idunion.mca_p2p1.f_II_II[kk];
             F[j2][i2]=F[i2][j2];
             S[i2][j2]=mcal->idunion.mca_p2p1.s_II_II[kk];
             S[j2][i2]=S[i2][j2];
            }
           }

        for (ib=0; ib<nt12; ib++) {
         for (jb=ib; jb<nt12; jb++) {
          block=dptr->idunion.mca_rl.block[ib][jb];
          diag=dptr->idunion.mca_rl.diag[ib][jb];

          if(block>=0) {

           i0=ns[ib];
           j0=ns[jb];
           for (i=1,kk=0; i<=n12[ib]; i++) {
            for (j=diag*i; j<=n12[jb]; j++,kk++) {
             i2=i0+i;
             j2=j0+j;
             F[i2][j2]=dptr->idunion.mca_rl.f[block][kk]*fact;
             F[j2][i2]=F[i2][j2];
             S[i2][j2]=dptr->idunion.mca_rl.s[block][kk]*fact;
             S[j2][i2]=S[i2][j2];
            }
           }
          }
         }
        }
        for (ib=0,ibt=0; ib<nt12; ib++,ibt++) {
         for (jb=0,jbt=nt12; jb<nt1234; jb++,jbt++) {
          block=dptr->idunion.mca_rl.block[ibt][jbt];
          diag=dptr->idunion.mca_rl.diag[ibt][jbt];

          if(block>=0) {

           i0=ns[ibt];
           j0=ns[jbt];
           for (i=1,kk=0; i<=n12[ib]; i++) {
            for (j=diag*i; j<=n1234[jb]; j++,kk++) {
             i2=i0+i;
             j2=j0+j;
             F[i2][j2]=dptr->idunion.mca_rl.f[block][kk]*fact;
             F[j2][i2]=F[i2][j2];
             S[i2][j2]=dptr->idunion.mca_rl.s[block][kk]*fact;
             S[j2][i2]=S[i2][j2];
            }
           }
          }
         }
        }
        for (ib=0,ibt=nt12; ib<nt1234; ib++,ibt++) {
         for (jb=ib,jbt=ibt; jb<nt1234; jb++,jbt++) {
          block=dptr->idunion.mca_rl.block[ibt][jbt];
          diag=dptr->idunion.mca_rl.diag[ibt][jbt];

          if(block>=0) {

           i0=ns[ibt];
           j0=ns[jbt];
           for (i=1,kk=0; i<=n1234[ib]; i++) {
            for (j=diag*i; j<=n1234[jb]; j++,kk++) {
             i2=i0+i;
             j2=j0+j;
             F[i2][j2]=dptr->idunion.mca_rl.f[block][kk]*fact;
             F[j2][i2]=F[i2][j2];
             S[i2][j2]=dptr->idunion.mca_rl.s[block][kk]*fact;
             S[j2][i2]=S[i2][j2];
            }
           }
          }
         }
        }

/*
*/
        fprintf(file,"S,F:\n");
        for(i=1;i<=nn;i++)
         for(j=1;j<=nn;j++)
          fprintf(file,"%d %d %e %e\n",i,j,S[i][j],F[i][j]);
        fflush(file);

        for(i=1;i<=nn;i++) {
         if(S[i][i]<ZERO) {
          printf("write_first_rl:\n");
          printf("S element is wrong: S(%d,%d) = %e\n",i,i,S[i][i]);
         }
        }

        fflush(stdout);
        Fm1 = matrix(1,nn,1,nn);
#ifdef IMSL
        tol=TOLSVD;
        ui = (_REAL *)get2(nn*nn,sizeof(_REAL));
        if(dptr->idunion.mca_rl.rigid==0) {
         for(i=1,k=0;i<=nn;i++)
          for(j=1;j<=nn;j++,k++)
           ui[k]=S[i][j];
        } else {
         for(i=1,k=0;i<=nn;i++)
          for(j=1;j<=nn;j++,k++)
           ui[k]=F[i][j];
        }
        w = imsl_d_lin_svd_gen(nn,nn,ui,
               IMSL_RANK,tol,&rank,
               IMSL_INVERSE,&v,0);
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          Fm1[i][j]=v[k];
         }
        }
        free(v);
        free(w);
#else
        if(dptr->idunion.mca_rl.rigid==0) {
         rank=invert(S,nn,nn,Fm1);
        } else {
         rank=invert(F,nn,nn,Fm1);
        }
#endif
        fprintf(file,"write_first_rl: rank of F= %d\n",rank); fflush(file);

        tmp = matrix(1,nn,1,nn);

/*
        fprintf(file,"Fm1*F:\n");
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=Fm1[i][k]*F[k][j];
          }
          fprintf(file,"%d %d %e\n",i,j,csum);
         }
        }
        fflush(file);
*/

        if(dptr->idunion.mca_rl.rigid==0) {

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=Fm1[i][k]*F[k][j];
          }
          tmp[i][j]=csum;
         }
        }

        } else {

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=Fm1[i][k]*S[k][j];
          }
          tmp[i][j]=csum;
         }
        }

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=S[i][k]*tmp[k][j];
          }
          Fm1[i][j]=csum;
         }
        }
        for(i=1;i<=nn;i++) {
         tol=1./S[i][i];
         csum=tol*Fm1[i][i];
/*         csum=RCmul(factns,csum); */
         fprintf(file,"tau(%d)= %e\n",i,csum);
        }
        fflush(file);

        }

        w=vector(1,nn);
#ifdef IMSL
        for(i=1,l=0;i<=nn;i++)
         for(j=1;j<=nn;j++,l++)
          ui[l]=tmp[i][j];
        wi = imsl_d_eig_gen(nn,ui,IMSL_VECTORS,&vi,0);
        if(dptr->idunion.mca_p2p1.rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1]=1./wi[i].re;
         }
        } else {
         for(i=0;i<nn;i++) {
          w[i+1]=wi[i].re;
         }
        }
        free(ui);
        free(wi);
#else
        vi=(_REAL *)get2(nn*nn,sizeof(_REAL));
#ifdef EISPACK
        eigvec_ev(nn,tmp,vi,w);
        if(dptr->idunion.mca_p2p1.rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1]=1./w[i+1];
         }
        }
#else
        eigvec_nr(tmp,nn,w,vi);
#endif
#endif

        fprintf(file,"Fm1*S eigenvalues:\n");
        nf=0;
        for(i=1;i<=nn;i++) {
         fprintf(file,"w(%d)= %e\n",i,w[i]);
         if(Nabs(w[i])>=ZERO) {
          nf++;
         }
        }
        fflush(file);

        dptr->idunion.mca_rl.nf = nf;
        nf=nn;
        dptr->idunion.mca_rl.eval= (_REAL *)get2(nf,sizeof(_REAL));
        dptr->idunion.mca_rl.evec= matrix(1,nn,1,nf);

        nf=0;
        for(i=1;i<=nn;i++) {
         dptr->idunion.mca_rl.eval[i-1]=1./w[i];
        }
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
#ifdef IMSL
          dptr->idunion.mca_rl.evec[i][j]=vi[k].re;
#else
          dptr->idunion.mca_rl.evec[i][j]=vi[k];
#endif
         }
        }
/*
        free_vector(w,1,nn);
        free(vi);
*/
/*
        for(i=1,k=0;i<=nn;i++,k++) {
         if(Nabs(wi[i-1].re)>=ZERO)
         {
          nf++;
          dptr->idunion.mca_rl.eval[nf-1]=1./wi[i-1].re;
          vt=&vi[k];
          for(j=1;j<=nn;j++) {
           dptr->idunion.mca_rl.evec[j][nf]=(*vt).re;
           vt+=nn;
          }
         }
        }
*/

/*
*/
        fprintf(file,"evec:\n");
        for(j=1;j<=5;j++) {
         for(i=1;i<=nn;i++) {
          fprintf(file,"%d %d %e\n",i,j,dptr->idunion.mca_rl.evec[i][j]);
         }
        }

        nf=nn;
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_rl.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_rl.evec[k][i]*tmp[k][j];
          }
          Fm1[i][j]=sum;
         }
        }
/*
        fprintf(file,"Ct*S*C:\n");
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          fprintf(file,"%d %d %e\n",i,j,Fm1[i][j]);
         }
        }
        fflush(file);
*/

/*
        for(i=1,k=0;i<=nf;i++) {
         for(j=1;j<=nf;j++,k++) {
          ui[k]=Fm1[i][j];
         }
        }
        w = imsl_d_eig_sym(nf,ui,IMSL_VECTORS,&v,0);
        for(i=1,k=0;i<=nf;i++) {
         if(Nabs(w[i-1]) > ZERO)
         {
          w[i-1]=1./l_sqrt(w[i-1]);
          for(j=1;j<=nf;j++,k++) {
           tmp[i][j]=v[k];
          }
         } else {
          w[i-1]=1.;
          for(j=1;j<=nf;j++,k++) {
           tmp[i][j]=0.;
          }
          tmp[i][i]=1.;
         }
        }
        for(i=1,m=0;i<=nf;i++) {
         for(j=1;j<=nf;j++,m++) {
          for(k=1,sum=0.;k<=nf;k++) {
           for(l=1,sum2=0.;l<=nf;l++) {
            sum2+=tmp[k][l]*w[l-1]*tmp[j][l];
           }
           sum+=dptr->idunion.mca_rl.evec[i][k]*sum2;
          }
          v[m]=sum;
         }
        }
        for(i=1,k=0;i<=nf;i++) {
         for(j=1;j<=nf;j++,k++) {
          dptr->idunion.mca_rl.evec[i][j]=v[k];
         }
        }
*/

        for(j=1,k=0;j<=nn;j++) {
         if(Nabs(Fm1[j][j]) > ZERO) {
          k++;
          sum=Fm1[j][j];
         } else {
          sum=1.;
         }
         sum=l_sqrt(sum);
         for(i=1;i<=nn;i++) {
          dptr->idunion.mca_rl.evec[i][j]/=sum;
         }
        }
        if(k!=nf) {
         printf("write_first: missed zero in Ct*S*C.\n");
        }

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_rl.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_rl.evec[k][i]*tmp[k][j];
          }
          S[i][j]=sum;
         }
        }
/*
        fprintf(file,"Ct*S1*C:\n");
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          fprintf(file,"%d %d %e\n",i,j,S[i][j]);
         }
        }
        fflush(file);
*/

        free_matrix(F,1,nn,1,nn);
        free_matrix(S,1,nn,1,nn);
        free_matrix(tmp,1,nn,1,nn);
}

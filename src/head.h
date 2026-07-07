#ifndef _HEAD_H
#define _HEAD_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*
*/
#define DOUBLE

#ifdef DOUBLE
#define _REAL		double
#define l_sqrt		sqrt
#define l_acos		acos
#define l_asin		asin
#define l_cos		cos
#define l_sin		sin
#else
#define _REAL		float
#define l_sqrt		sqrtf
#define l_acos		acos
#define l_asin		asin
#define l_cos		cos
#define l_sin		sin
#endif

#define MAXIMUM(A, B)  (((A) > (B) ? (A) : (B)))
#define MINIMUM(A, B)  (((A) < (B) ? (A) : (B)))
#define Nabs(dumx) ((dumx) >= 0 ? (dumx) : -(dumx))
#define Nint(dumx)  ( ((dumx) >=0.0) ? (int)((dumx) + 0.5) : (int)((dumx) - 0.5) )
#define SIGN(a,b) ((b) > 0 ? Nabs(a) : -Nabs(a))
#define rdelta(x,y) ( (x==y) ? 1. : 0. )
#define idelta(x,y) ( (x==y) ? 1 : 0 )

typedef struct rdiad {
        _REAL sq,xx,xy,xz,yx,yy,yz,zx,zy,zz;
} rdiadstruct;

#include "complex.h"
#include "nr.h"
#include "nrutil.h"
#ifdef IMSL
#include <imsl.h>
#endif

#include  "prm.h"

#define TRUE		1
#define FALSE		0

#define _PI 	3.141592653589793
#define TWOPI	6.28318530717958
#define HALFPI	1.570796326794895
#define TODEG 	57.295779513082379	/* (180.0 / PI) */
#define TORAD 	.017453292519943	/* (PI / 180.0) */

#define AVOG   6.02205
#define BOLTZ  8.31441

#define TOLSVD  1.e-10
#define ZERO    1.e-08
 
#define TOKLEN		81
#define ATOMTOKLEN	 4
 
#define PARM             1
#define STREAM           2
#define CRD              3
#define RST              4
#define PDB              5
#define MOVIE            6
#define XYZ              7
#define ATOM             8
#define TABLE            9
#define BOND            10
#define ANGLE           11
#define DIHE            12
#define GR              13
#define MOL             14
#define APE             15
#define STDD            16
#define APE_TO_DEC      17
#define RIG             18
#define ORDERB          19
#define DENSITY1        20
#define DISCOVER2       21
#define DISCOVER3       22
#define GYR             23
#define BEAD            24
#define MCARIG          25
#define MCA_P2P1        26
#define ACFF            27
#define STDF            28
#define GROUP           29
#define GROUPCENTER     30
#define GROUPCOM        31
#define GROUPGYR        32
#define ASA             33
#define MCA_RL          34
#define ORAC3_DUMP      35
#define ORAC4_DUMP      36
#define MCA_P2P1L       37
#define MCA_P1          38
#define STATIC          39
#define PUCKER          40
#define NUCLEIC         41
#define RMS             42
#define ROTNE_PRAGER    43
#define OTTINGER        44
#define ZR              45
#define OSEEN           46
#define DLPOLY          47
#define MCA_P2P1MM      48
#define ORIENT          49
#define STREAMW         50
#define MCA_P1L         51
#define CCFF            52
#define CHARMM          53
#define PERS            54
#define ASAF            55
#define MAP             56
#define MCA_P1LL        57
#define MCA_P2          58
#define MCA_P2LL        59
#define CBMC            60
#define MUCA            61
#define MCA_DIHE        62
#define HBOND           63
#define MCA_P1ML        64
#define MCA_DOT         65
#define NWCHEM          66
#define NWCHEM_RST      67
#define POS             68
#define AXSF            69
#define AMBER           70

/* bit-wise flags for idunion.group */

#define COMFLAG         0x1     /* 1st bit */
#define GEOMFLAG        0x2     /* 2nd bit */
#define GYRFLAG         0x4     /* 3rd bit */

typedef struct name {
	struct name	*next;
	char		*name;
} namestruct;
 
typedef struct set {
	struct set	*next;
        parmstruct      *prm;
	_REAL		*crd;
        int             max;
	_REAL		*boxc,*boxm1;
	double		Etot,Epot,Ekin,Press,*Enerd;
	_REAL		dir[3],P2;
	_REAL		rotdir[3][3],**cart_to_cell,**cell_to_cart;
        float           *crd_s;
} setstruct;

typedef struct setr {
	_REAL		*crd;
	_REAL		*boxc;
} setstructr;

typedef struct select {
	int	item;
	int	code;
} selectstruct;

#define SEL 1
#define SELTILL -1
#define SELTILLEND -2
#define SELDONE -3
 
typedef struct file {
	struct file	*next;
	char		f_name[TOKLEN];
	FILE		*f_file;
} filestruct;
 
typedef struct print {
	struct print	*next;
	int		p_type;
	int		p_subtype;
	int		p_nitems;
	char		*p_value,*p_value1,*p_value2,*p_value3,*p_value4;
} printstruct;

typedef struct stats {
        int     ct,weight;
        _REAL   avg, avg2, max, min;
} statstruct;

typedef struct Cstats {
        int     ct;
        _REAL   x, y;
} Cstatstruct;

typedef struct bond {
        int id1,id2;
        _REAL zeta1,zeta2;
        _REAL x1,y1,z1,x2,y2,z2;
        _REAL x,y,z,sq,q0,q1,q2,q3;
        int *offset_l,*doffset;
        int *offset_l2,*doffset2;
//        fcomplex D1[3],dD1[3][3];
//        fcomplex D2[5],dD2[5][3];
//        fcomplex D4[9],dD4[9][3];
} bondstruct;

typedef struct dihe {
        int id[4];
        _REAL zeta[4],*r[4];
        _REAL x[4],y[4],z[4],fx[4],fy[4],fz[4];
        _REAL base[21],dbase_dc[21];
} dihestruct;

typedef struct rot2 {
        /* fcomplex dD2[5][3]; */
        _REAL D2[5];
} rot2struct;

typedef struct rot1 {
        _REAL f[5][3];
} rot1struct;

typedef struct rot {
        _REAL f[3][3];
} rotstruct;

typedef struct mod {
        _REAL mod,m[3];
        fcomplex D2[5];
} modstruct;

typedef struct mod1 {
        _REAL phi[3],dphidm1[3][3];
} mod1struct;

typedef struct dot {
        _REAL phi,dphidm1[3],dphidm2[3],mod;
} dotstruct;

typedef struct mod12 {
        _REAL phi[5],dphidm1[5][3],dphidm2[5][3],mod;
} mod12struct;

typedef struct mod123 {
        _REAL phi[3],dphidm1[3][3],dphidm2[3][3],
                     dphidm3[3][3];
} mod123struct;

typedef struct mod1234 {
        _REAL phi[5],dphidm1[5][3],dphidm2[5][3],
                     dphidm3[5][3],dphidm4[5][3],mod;
} mod1234struct;

typedef struct mod123456 {
        _REAL phi[5],dphidm1[5][3],dphidm2[5][3],
                     dphidm3[5][3],dphidm4[5][3],
                     dphidm5[5][3],dphidm6[5][3];
} mod123456struct;

typedef struct id {
	struct id	*next;
	char		id_name[TOKLEN];
	int		id_type;
	int		id_subtype;
	int		id_setup;
	union {
		struct {
                        int             fmt,conv;
			int		ifile;
			int		iset;
			int		fileset;
			int		nfiles;
			namestruct	*names;
			char		*name;
			FILE		*fileptr;
			int		skip, win, winct, type;
                        _REAL           w;
                        _REAL           w2[9];
		} serw;

		struct {
                        parmstruct      *prm;
                        int             fmt,conv,reconstr,skipa,skipb,iweight;
			int		ifile;
			int		iset;
			int		fileset;
			int		nfiles;
			int		nstates;
                        int             nproc,shift;
                        int             ibox,scale;
			namestruct	*names;
			char		*name;
			FILE		*fileptr;
			int		skip, win, winct;
			int		atom_selected;
			setstruct	*set;
                        struct id       *w;
		} ser;

                struct {
                        char            *name;
                        parmstruct      *prm;
                } parm;

                struct {
                        int             natoms;
                        int             nsets;
                        setstruct       *sets;
                } stat;

                struct {
                        int             natoms;
                        int             fmt,conv;
                        int             blank;
                        char            name[TOKLEN];
                        FILE            *file;
                        int             mod;
                        selectstruct    *select;
                        int             nsel;
                        int             modsel;
                        int             selmod;
                        setstruct       *set;
                        setstructr      *setr;
                        parmstruct      *prm;
                        int             option;
			_REAL           weight;
                } crd;

                struct {
                        int             nitems;
                        int             append;
                        int             wrap;
                        char            format[TOKLEN];
                        int             itemsline;
                        int             sofar;
                        char            linedelim;
                        int             nfiles;
                        char            filename[TOKLEN];
                        filestruct      *files;
                        filestruct      *curfile;
                        FILE            *fptr;
                        printstruct     *print;
                } tab;

                struct {
			FILE		*file;
			char		fname[TOKLEN];
                        int             type;
                        int             nitems,nlist;
                        int             nboxes;
                        _REAL           scale,min,max;
                        _REAL           *boxes,casesw;
                        int             *list,jump;
                        _REAL           *r1,*r2,*r3,*r4;
                        int             cases,callcase;
                        int             nelem,nmax;
                        int             glrate;
                        int             *statesp,*states;
                        int             **ratemat;
                        _REAL           *value;
                        int             *gauche;
                        union {
                                struct {
                                        _REAL           value;
                                        statstruct      stat;
                                } scalar;
                        } resunion;
                } work;

                struct {
			FILE		*file;
			char		fname[TOKLEN];
                        int             type;
                        struct id       *grp1,*grp2;
                        int             nboxes,nboxesa,nboxesb;
                        _REAL           scalea,mina,maxa;
                        _REAL           scaleb,minb,maxb;
                        int             *boxes;
                        int             cases;
                } orient;

                struct {
			FILE		*file;
                        parmstruct      *prm;
			char		fname[TOKLEN];
                        int             type;
                        int             nboxes,nboxesa,nboxesb;
                        _REAL           scalea,mina,maxa;
                        _REAL           scaleb,minb,maxb;
                        _REAL           *boxes;
                        int             jump;
                        _REAL           *r1,*r2,*r3,*r4;
                        int             cases,callcase;
                        _REAL           *value;
                        union {
                                struct {
                                        _REAL           value;
                                        statstruct      stat;
                                } scalar;
                        } resunion;
                } work2;

                struct {
			FILE		*file;
                        parmstruct      *prm;
			char		fname[TOKLEN];
                        int             type;
                        int             nlist;
                        _REAL           **r1,**r2;
                        int             cases;
                        _REAL           *u21,*r21,*umat,gyr;
                        int             udim;
                } work3;

                struct {
                        parmstruct      *prm;
                        int             nbond,nbead,hydrod;
                        _REAL           sigma,factns,stick_slip;
                        _REAL           *zt,zav,zr,bav,bavt,bavt2,lambda,cutsq,screen;
                        _REAL           **r1,**r2,*z1,*z2;
                        bondstruct      *l;
                        // rdiadstruct     *g,**M1,**F1;
                        _REAL           **r,*z;
                } bead;

                struct {
			FILE		*file;
			char		fname[TOKLEN];
                        struct id       *bead;
                        int             type;
                        _REAL           Dxx,Dyy,Dzz,qdir[4];
                        _REAL           eigenstates1[4][4],eigenvalues1[4];
                        _REAL           eigenstates2[6][6],eigenvalues2[6];
                } rig;

                struct {
			FILE		*file;
                        struct id       *bead;
			char		fname[TOKLEN];
                        int             type,dump,rigid,sel_modes,sel_bonds;
                        int             cases;
                        int             nf,n1,n12,n2,dim1,dim12,dim2,nn;
                        _REAL           *s1,*f1,*s12,*f12,*s2,*f2,casesw;
                        rotstruct       *st1,*ft1,*st12,*ft12,*st2,*ft2;
                        _REAL           **Rm1;
                        _REAL           **evec,*eval;
                        _REAL           *f0;
                        _REAL           **evec0,*eval0;
                        mod123struct    *l123;
                        rotstruct       **Fb;
                        int             tens,list,nlist,*j1,*j2,*j3;
                } mcarig;

                struct {
			FILE		*file;
                        struct id       *mca;
                        struct id       *bead;
			char		fname[TOKLEN];
                        int             type,dump,rigid;
                        int             cases;
                        int             m1,m2,n0,n12,n1234,n123456,nf,
                                        dim_I_I,dim_I_II,dim_II_II,
                                        dim_I_III,dim_II_III,dim_III_III;
                        _REAL           *s_I_I,*s_I_II,*s_II_II,
                                        *s_I_III,*s_II_III,*s_III_III,
                                        *S1,*eval,**evec,casesw;
                        _REAL           *f_I_I,*f_I_II,*f_II_II,
                                        *f_I_III,*f_II_III,*f_III_III;
                        rotstruct       **Fb,**Fm;
                        modstruct       *m;
                        mod12struct     *l12;
                        mod1234struct   *l1234;
                        mod123456struct   *l123456;
                } mca_p2p1;

                struct {
			FILE		*file;
			char		fname[TOKLEN];
                        parmstruct      *prm;
                        struct id       *model;
                        int             type,typo,dump,printj;
                        _REAL           **r1,**r2,**r3,**r4;
                        int             nbond,nbase,nbase2,neval,cases;
                        _REAL           **proj,**proj2,**proj3,**proj4,**proj5,**proj6,
                                        *modsq,*rtmp,*btmp,*tmp,casesw,
                                        **relax_Sz,**relax_Sxy,**relax_Iz_Sz,
                                        **nOe_Iz_Sz,*Ssq;
                        fcomplex        **projc;
                        int             nt,nzero,nomega;
                        _REAL           domega,dtime,T,*gam1,*gam2,*gNgH,*csa,gamH,
                                        **omegaN,**omegaH,omega0[200],
                                        epsilon,tzero;
                        _REAL           *tau_p2,*p2t,*workft;
                } acff;

                struct {
                        int             type;
                        char            fname1[TOKLEN],fname2[TOKLEN];
                        FILE            *file1,*file2;
                        _REAL           max,min,maxz,minz,scale,scalez;
                        _REAL           rho;
                        _REAL           *boxav,smallside;
                        int             nboxes,nboxesz;
                        _REAL           *boxes,*boxespa,*boxespe,casesw;
                        _REAL           *listi,*listj;
                        int             jump;
                        int             win, winct;
                        int             cases,grid;
                        _REAL           *rmin,*rmax;
                } db;   /* distribution */

                struct {
                        setstruct       *set;
                        int             natoms, maxatom, flag;
                        int             *Atomlist, *AtomXlist;
                        _REAL           total_mass;
                        _REAL           geomcenter[3];
                        _REAL           masscenter[3],q[4];
                        _REAL           momin, radgyr;
                        _REAL           xmin,xmax,ymin,ymax,zmin,zmax;
                        _REAL           *vdw;
                } group;

               struct {
                        FILE            *file;
			char		fname[TOKLEN];
                        int             npuckers;
                        int             npoints;
                        _REAL           **points;       /* pointers to xyz */
                        _REAL           *scratch;       /* copy of xyz */
                        _REAL           *results;
                        statstruct      *stats;         /* alt: stat/Cstat */
                        int             *boxes,nboxes; 
                        _REAL           min,max,scale;
                        int             cases;
                } puck;

                struct {
                        struct id       *grp;
                        int             natoms;
                        _REAL           *refcrd, *refcrdctr, *crd2, *crd2a;
                        _REAL           refcenter[3];
                        _REAL           rms;
                        statstruct      stat;
                        int             fit, it0, ok, refstatic;
                } rms;

                struct {
                        _REAL           f,res,probe;
                        int             nmax,close,output;
                        struct id       *group,*groupc;
                        FILE            *fpdb;
                } asa;

                struct {
                        _REAL           fsurf,res,probe,cell;
                        int             nmax,output;
                        FILE            *fpdb;
                        _REAL           *vdw,*vdw2,*sum2,*diff2,*sigma;
                } asaf;

                struct {
                        int Nx,Ny,Nz,Mx,My,Mz,NN,ncell,mapsize,Mstar;
                        int *offset,*Map,*Head,*List;
                        _REAL cell,cellx,celly,cellz,cellxm1,cellym1,cellzm1,
                              rcut,box[3];
                        setstruct *set;
                } map;

                struct {
                        FILE            *file;
                        parmstruct      *prm;
                        char            fname[TOKLEN];
                        int             center,offset,distmax,cases,nn;
                        _REAL           *distr;
                        int             *count;
                } pers;

                struct {
                        char            filename[TOKLEN];
                        FILE            *table;
                        FILE            *list;
                        int             *BondHAt1,*BondHAt2, Nbonh;
                        int             *donor1, *donor2, *acceptor,
                                        ndonor, nacceptor,
                                        solv1;
                        setstruct       *set;
                        _REAL           cut, angle;
                        statstruct      *distat, *angstat;
                        int             *stat;
                        int             *AtomRes;
                } hb;   /* hydrogen bond */

	} idunion;
 
} idstruct;


#ifdef NUCLEIC_DEF
/* static char     *NucleicNames[] = {"CYT ", "GUA ", "THY ", "ADE ", "URA "}; */
static char     *NucleicNames[] = {"DC  ", "DG  ", "DT  ", "DA  ", "U   ",
                                   "C   ", "G   ", "A   ", 
                                   "C3  ", "G3  ", "A3  ", 
                                   "C5  ", "G5  ", "A5  ", 
                                   "DC3 ", "DG3 ", "DT3 ", "DA3 ", "U3  ",
                                   "DC5 ", "DG5 ", "DT5 ", "DA5 ", "U5  "};
#define N_NUCLEIC      24
#endif


#ifdef AMINO_DEF
static char     *AminoNames[] = {
        "ALA ", "GLY ", "SER ", "THR ", "LEU ", "ILE ", "VAL ", "ASN ",
        "GLN ", "ARG ", "HID ", "HIE ", "HIP ", "TRP ", "PHE ", "TYR ",
        "GLU ", "ASP ", "LYS ", "PRO ", "CYS ", "CYX ", "MET "};
#define N_AMINO 23
#endif

#ifdef SINGLE_DEF
static char     *SingleNames[] = { "HB  ", "HE  ", "OME "};
#define N_SINGLE        3
#endif

extern void     gettok(), inerr();
extern void *get2(int,int);
extern char     *get(), *syserr();
extern setstruct       *getstatset();
extern char tok[TOKLEN];
extern char prmtop[TOKLEN];

#endif

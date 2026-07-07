/*
 *  main.c - main() + higher level input/setup routines
 */

#include "head.h"
 
void usage(char *arg0)
{
	fprintf (stderr, "Usage: %s [-i analin] [-o analout]\n",
			arg0);
	fprintf (stderr, "       %s < analin > analout\n");
	exit(1);
}

/***********************************
	global vars
************************************/
 
// char prmtop[TOKLEN], tok[TOKLEN];
char tok[TOKLEN];
// char		analin[TOKLEN], analout[TOKLEN], prmtop[TOKLEN];
char		analin[TOKLEN], analout[TOKLEN];
statstruct      nullstat = {0.0, 0.0, 0.0};

int first=0;

/***********************************
	external vars
************************************/
 
/****	stream.c   ****/
 
extern idstruct	*firstser;	/* default crd set for measurements */
extern idstruct	*firstserw;
 
/****	util.c	****/
 
// extern char	tok[];
char	tok[];
extern char	line[];

extern idstruct	*startid();
 
 
/***********************************
	forward declarations
************************************/
 
void	init(), init2(), init3(),readin(), readin2(), files_in(), files_out(),
        declare(), declare2(),output(),readin3(),declare3();
 
/***********************************************************************
 							MAIN()
************************************************************************/
 
int main(int argc, char *argv[])
{
        init(argc, argv);
        doit();
/*
*/
        init2();
        /* doit2(); */
/*
*/
        init3();
        if(first>0) {
         doit3();
        }
	return(0);
}
 
/***********************************************************************
 							INIT()
************************************************************************/
 
/*
 * init() - get filenames if any from command line and read in analin
 */
 
void
init(int argc, char *argv[])
{
	fileargs(argc, argv);
	readin();
 
}

void
init2()
{
	readin2();
}
 
void
init3()
{
	readin3();
}
 
 
/***********************************************************************
 							READIN()
************************************************************************/
 
/*
 *  readin() - read each section of analin, setting up idstructs
 */
 
void
readin()
{
	line[0] = '\n';
        files_in();
        files_out();
        declare();
	output();
 
	/*
	 *  if last token is "STOP," quit here 
	 *	use scanf() rather than gettok() to avoid error if eof
	 */
 
/*
	scanf(" %s", tok);
	if (!strcmp(tok, "STOP")) {
		printf("STOP\n");
		exit(0);
	}
*/

	/*
	 *  in case it will be a long run, make sure user can see
	 *	the interpretation of input
	 */

	fflush(stdout);
}

void
readin2()
{
 printf("%s\n",tok);
        declare2();
 
/*
	scanf(" %s", tok);
	if (!strcmp(tok, "STOP")) {
		printf("STOP\n");
		exit(0);
	}
*/
	fflush(stdout);
}
 
void
readin3()
{
 printf("%s\n",tok);
        declare3();
 
	scanf(" %s", tok);
	if (!strcmp(tok, "STOP")) {
		printf("STOP\n");
		exit(0);
	}
	fflush(stdout);
}

/***********************************************************************
 							FILES_IN()
************************************************************************/
 
/*
 *  files_in() - read & set up input file structs: STREAM
 */
 
void
files_in()
{
	int	nser = 0;
	int	nserw = 0;
 
	/*
	 *  consider analin input
	 */

	gettok();
	if (strcmp(tok, "FILES_IN") )
		inerr("beginning: expected FILES_IN","");
 
	for (;;) {
		gettok();
		if (!strcmp(tok, "PARM")) {
		 gettok();
		 if (!strcmp(tok, "FORM")) {
		  gettok();
		  if (!strcmp(tok, "AMBER")) {
            		rdprm();
		  } else if (!strcmp(tok, "PDB")) {
            		rdprm_pdb();
		  } else if (!strcmp(tok, "DISCOVER3")) {
            		rdprm_his3();
		  } else if (!strcmp(tok, "DISCOVER2")) {
            		rdprm_his2();
                  } else {
                   printf("main: parm form required.\n");
                   exit(1);
                  }
                 } else {
                  printf("main: parm form required.\n");
                  exit(1);
                 }
		} else if (!strcmp(tok, "STREAMW")) {
			nserw++;
            		setupstreamw();
		} else if (!strcmp(tok, "STREAM")) {
			nser++;
            		setupstream();
		} else if (!strcmp(tok, "STATIC")) {
            		setupstat();
		} else if (!strcmp(tok, "FILES_OUT") ) {
			break;
		} else
			inerr("FILES_IN: unexpected token","");
	}
	if (!nser) 
		inerr("Must specify at least 1 STREAM file in FILES_IN","");
}
 
 
/***********************************************************************
 							FILES_OUT()
************************************************************************/
 
/*
 *  files_out() - read & set up output file structs
 */
 
void
files_out()
{
 
	for(;;) {
		gettok();
		if (!strcmp(tok, "DECLARE")) {
			return;
		} else if (!strcmp(tok, "MCARIG")) {
			init_mcarig();			
		} else if (!strcmp(tok, "MCA_P2P1")) {
			init_mca_p2p1();
		} else if (!strcmp(tok, "CRD")) {
			initcoord();
		} else if (!strcmp(tok, "BOND")) {
			initbond();
		} else if (!strcmp(tok, "ANGLE")) {
			initang();
		} else if (!strcmp(tok, "DIHE")) {
			initdihe();
		} else if (!strcmp(tok, "PUCKER")) {
			initpucker();
		} else if (!strcmp(tok, "GR")) {
			initgr();
		} else if (!strcmp(tok, "ORDERB")) {
			initorderb();
		} else if (!strcmp(tok, "PERS")) {
			initpers();
		} else if (!strcmp(tok, "ORIENT")) {
			initorient();
		} else if (!strcmp(tok, "ACFF")) {
			init_acff();			
		} else if (!strcmp(tok, "DENSITY1")) {
			initdensity1();
                } else if (!strcmp(tok, "HBOND")) {
                        inithbond();
		} else if (!strcmp(tok, "TABLE")) {
			inittable();
		} else {
			inerr("unexpected token in FILES_OUT","");
		}
	}
}
 
 
/***********************************************************************
 							DECLARE()
************************************************************************/
 
/*
 * declare() - set up structs for measurements
 */
 
void
declare()
{
	for (;;) {
 
		gettok();
		if (!strcmp(tok, "OUTPUT")) {
			return;			
		} else if (!strcmp(tok, "BOND")) {
			setupbond(2);			
		} else if (!strcmp(tok, "ANGLE")) {
			setupang(3);			
		} else if (!strcmp(tok, "DIHE")) {
			setupdihe();			
		} else if (!strcmp(tok, "PUCKER")) {
			pucker();			
		} else if (!strcmp(tok, "RMS")) {
			rms();			
                } else if (!strcmp(tok, "GROUP")) {
                        group();
                } else if (!strcmp(tok, "HBOND")) {
                        setuphb();
                } else if (!strcmp(tok, "ASA")) {
                        setup_asa();
                } else if (!strcmp(tok, "ASAF")) {
                        setup_asaf();
		} else if (!strcmp(tok, "GR")) {
			setupgr();			
		} else if (!strcmp(tok, "ORDERB")) {
			setuporderb();			
		} else if (!strcmp(tok, "PERS")) {
			setuppers();			
		} else if (!strcmp(tok, "MAP")) {
			init_map();			
		} else if (!strcmp(tok, "ORIENT")) {
			setuporient();			
		} else if (!strcmp(tok, "BEADS")) {
			setup_bead_fric();			
		} else if (!strcmp(tok, "BONDS")) {
			setup_bead_bond();			
		} else if (!strcmp(tok, "RIG")) {
			setup_rig();
		} else if (!strcmp(tok, "MCARIG_M1")) {
			setup_mcarig_m1();
		} else if (!strcmp(tok, "ACFF")) {
			setup_acff();			
		} else if (!strcmp(tok, "GYR")) {
			setupgyr();			
		} else if (!strcmp(tok, "DENSITY1")) {
			setupdensity1();
		} else {
			inerr("DECLARE: unexpected token","");
		}
	}
}

void
declare2()
{
	for (;;) {
 
		gettok();
/*
		if (!strcmp(tok, "END")) {
			return;			
*/
		if (!strcmp(tok, "DECLARE3")) {
			return;			
		} else if (!strcmp(tok, "ACFF")) {
                        first++;
			setup_acff();			
		} else {
			inerr("DECLARE2: unexpected token","");
		}
	}
}
 
void
declare3()
{
	for (;;) {
 
		gettok();
		if (!strcmp(tok, "END")) {
			return;			
		} else if (!strcmp(tok, "MCA_P2P1")) {
                        first++;
			setup_mca_p2p1();
		} else {
			inerr("DECLARE3: unexpected token","");
		}
	}
}
 
 
/***********************************************************************
 							OUTPUT()
************************************************************************/
 
/*
 * output() - set up previously declared output structs, coordinating
 *	with DECLAREd things if necc
 */
 
void
output()
{
	for(;;) {
		gettok();
 
/*
*/
		if (!strcmp(tok, "END")) {
			return;			
		} else if (!strcmp(tok, "DECLARE2")) {
			return;			
		} else if (!strcmp(tok, "CRD")) {
			setupcoord();
		} else if (!strcmp(tok, "TABLE")) {
			setuptab();
		} else
			inerr("DO: unexpected token","");
	}
}

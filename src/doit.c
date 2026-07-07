#include "head.h"

/*
 *  doit.c - run the analysis on the input stream(s)
 */

/**** 	global   ****/

int	it;	/* iteration */
int	itw;
int	lastpass = 0;

/****	extern   ****/

extern idstruct	*idlist;

/***********************************************************************
                                                        DOIT()
************************************************************************/

/*
 * doit() - iterate: reading in streams, measuring, outputting
 */

void doit()
{
  void printsumm();
	/*
	 * process streams then print summary
	 */

/*
	idstruct	*idptr;
	printf("PROCESSING ORDER: types\n");
	for (idptr=idlist; idptr!=NULL; idptr=idptr->next) 
		printf("%d ", idptr->id_type);
	printf("\n");
	fflush(stdout);
*/
	for(it=0; procq(); it++,itw++)
		;
	lastpass++;
	printsumm();
}

void doit2()
{
	for(it=0; procq2(); it++)
		;
}

void doit3()
{
	for(it=0; procq3(); it++)
		;
}

/***********************************************************************
                                                        PROCQ()
************************************************************************/

/*
 *  procq() - traverse idlist, taking appropriate actions: reading streams,
 *	performing measurements, outputting results. If no stream yields a
 *	set, return done=0
 */

procq()
{
	idstruct	*idptr,*serptr,*mapptr;
	int		someser = 0, someserw=0;
        int j;

	for (idptr=idlist; idptr!=NULL; idptr=idptr->next) {
/*
		printf("TYPE %d\n", idptr->id_type); fflush(stdout);
*/
		switch (idptr->id_type) {
			case STATIC:
			case PARM:
				break;
			case STREAMW:
				someserw = nextserw(idptr);
				break;
			case STREAM:
				someser = nextser(idptr);
                                serptr = idptr;
				if(someser) 
				{
				 if(idptr->idunion.ser.reconstr) {
                                  sputinbox(serptr);
                                 }
/*
				 inertdir(serptr);
                                 dipol(serptr);
*/
                                }
				break;
			case MAP:
                                if (someser) {
                                        setup_map(idptr);
                                        mapptr=idptr;
                                }
                                break;
			case GROUP:
                                if (someser)
                                        dogrp(idptr);
                                break;
			case BOND:
				if (someser)
					distribbond(idptr,serptr);
				else
					writebond(idptr);
				break;
			case ANGLE:
				if (someser)
					distribang(idptr,serptr);
				else
					writeang(idptr);
				break;
			case DIHE:
				if (someser)
					distribdihe(idptr,serptr);
				else
                                 if(idptr->idunion.work.type==1)
					writetorch(idptr);
                                 else
					writedihe(idptr);
				break;
                        case PUCKER:
                                if (someser)
                                        calcpucker(idptr);
                                else
                                        writepucker(idptr);
                                break;
                        case RMS:
                                if (someser)
                                        fitrms(idptr);
                                break;
                        case HBOND:
                                if (someser)
                                        // hbond(idptr);
					hbond(idptr,serptr);
                                else
                                        closehb(idptr);
                                break;
			case GR:
				if (someser)
                                 if(idptr->idunion.db.grid==0)
                                 {
					gr(idptr,serptr);
                                 }
                                 else
                                 {
					gr_grid(idptr,serptr);
                                 }
				else
					writegr(idptr);
				break;
			case ASA:
				if (someser)
					do_asa(idptr);
				break;
			case ASAF:
				if (someser)
					surf_mol(idptr,serptr,mapptr);
				break;
			case ORDERB:
				if (someser)
					distriborderb(idptr,serptr);
				else
					writeorderb(idptr);
				break;
			case PERS:
				if (someser)
					do_pers_flory(idptr,serptr);
				else
					write_pers_flory(idptr);
				break;
			case ORIENT:
				if (someser)
					distriborient(idptr);
				else
					writeorient(idptr);
				break;
			case GYR:
				if (someser)
					gyr(idptr,serptr);
				break;
			case TABLE:
				if (someser)
					writetab(idptr);
				else continue;
/*
					closetab(idptr);
*/
				break;
			case BEAD:
				break;
			case RIG:
				if (someser) {
                                 rotdiff(idptr,serptr);
                                } else break;
                                break;
			case MCARIG:
				if (someser) {
                                 if(idptr->idunion.mcarig.dump>=0)
				 {
                                  if(idptr->idunion.mcarig.list==0)
				  {
                                   if(idptr->idunion.mcarig.tens==0)
				   {
				 	do_first_m1(idptr,serptr);
				   }
                                   else
				   {
					do_first_m1_tens(idptr,serptr);
				   }
				  }
                                  else
				  {
				 	do_first_m1_list(idptr,serptr);
				  }
				 }
				 else
				 {
					write_first(idptr);
					someser=0;
				 }
				} else {
/*
*/
					write_first(idptr);
/*
                    serptr->idunion.ser.nstates=1000000;
*/
/*
					write_first_tril(idptr);
*/
                                }
				break;
			case MCA_P2P1:
				break;
			case ACFF:
				break;
			case DENSITY1:
				if (someser)
					distribdensity1(idptr,serptr);
				else
					writedensity1(idptr);
				break;
			case CRD:
				if (someser)
					writecrd(idptr,serptr);
				else
					break;
				break;
			default:
				printf("Program error: unexpected id_type %d\n",
						idptr->id_type);
				exit(1);
		}
	}
	return(someser);
}


/***********************************************************************/
procq2()
{
	idstruct	*idptr,*serptr;
	int		someser = 0;
	int		someserw = 0;
        int j;


	for (idptr=idlist; idptr!=NULL; idptr=idptr->next) {
/*
		printf("TYPE %d\n", idptr->id_type); fflush(stdout);
*/
		switch (idptr->id_type) {
			case PARM:
				break;
			case STREAMW:
				someserw = nextserw(idptr);
				break;
			case STREAM:
				someser = nextser(idptr);
                                serptr = idptr;
				if(someser)
		                {
				 if(idptr->idunion.ser.reconstr) 
				 {
                                  sputinbox(serptr);
				 }
//				 inertdir(serptr);
				}
				break;
			case BOND:
				break;
			case ANGLE:
				break;
			case DIHE:
				break;
			case GR:
				break;
			case ORDERB:
				break;
			case PERS:
				break;
			case ORIENT:
				break;
			case GYR:
				break;
			case TABLE:
				break;
			case BEAD:
				break;
			case MCARIG:
				break;
			case MCA_P2P1:
				break;
			case ACFF:
				break;
			case DENSITY1:
				break;
			case CRD:
				break;
			default:
				printf("Program error: unexpected id_type %d\n",
						idptr->id_type);
				exit(1);
		}
	}
	return(someser);
}

/***********************************************************************/
procq3()
{
	idstruct	*idptr,*serptr;
	int		someser = 0;
	int		someserw = 0;
        int j;


	for (idptr=idlist; idptr!=NULL; idptr=idptr->next) {
/*
		printf("TYPE %d\n", idptr->id_type); fflush(stdout);
*/
		switch (idptr->id_type) {
			case PARM:
				break;
			case STREAMW:
				if (it==0)
				{
                                 rewind(idptr->idunion.serw.fileptr);
				}
				someserw = nextserw(idptr);
				break;
			case STREAM:
				if (it==0)
				{
                                 rewind(idptr->idunion.ser.fileptr);
				}
				someser = nextser(idptr);
                                serptr = idptr;
				if(someser)
				{
				 if(idptr->idunion.ser.reconstr) 
				 {
                                  sputinbox(serptr);
				 }
//				 inertdir(serptr);
				}
				break;
			case BOND:
				break;
			case ANGLE:
				break;
			case DIHE:
				break;
			case GR:
				break;
			case ORDERB:
				break;
			case PERS:
				break;
			case ORIENT:
				break;
			case GYR:
				break;
			case TABLE:
				break;
			case BEAD:
				break;
			case MCARIG:
                                break;
			case MCA_P2P1:
				if (someser) {
					do_first_p2p1(idptr,serptr);
				} else {
					write_first_p2p1(idptr);
                                }
				break;
			case ACFF:
				if (someser) {
                                 if (it==0) adjust_acff(idptr);
                                 if(idptr->idunion.acff.model->id_type==MCA_P2P1) 
				 {
				  do_acff2_stress(idptr,serptr);
                                 } 
				} else {
                                  write_acff_stress(idptr);
                                }
				break;
			case DENSITY1:
				break;
			case CRD:
				break;
			default:
				printf("Program error: unexpected id_type %d\n",
						idptr->id_type);
				exit(1);
		}
	}
	return(someser);
}
/***********************************************************************
                                                        PRINTSUMM()
************************************************************************/

/*
 *  printsumm() - print summary statistics on measurements
 */

void printsumm()
{
	idstruct	*idptr;

	printf("SUMMARY\n");
	for (idptr=idlist; idptr!=NULL; idptr=idptr->next) {
/*
		printf("TYPE %d\n", idptr->id_type); fflush(stdout);
*/
		switch (idptr->id_type) {
		case STATIC:
		case HBOND:
		case PARM:
		  break;
		case STREAMW:
		  break;
		case STREAM:
		  break;
		case BOND:
                  printf("--Distance %s: ", idptr->id_name);
                  printstat(&idptr->idunion.work.resunion.scalar.stat, 0.0, 0);
		  break;
		case ANGLE:
                  printf("--Angle %s: ", idptr->id_name);
                  printstat(&idptr->idunion.work.resunion.scalar.stat, 0.0, 1);
		  break;
		case GYR:
		  break;
		case TABLE:
		  break;
		case DIHE:
                  printf("--Dihedral %s: ", idptr->id_name);
                  printstat(&idptr->idunion.work.resunion.scalar.stat, 0.0, 1);
		  break;
                case PUCKER:
                  printpucker(idptr);
                  break;
                case RMS:
                  printf("--RMS %s: ", idptr->id_name);
                  printstat(&idptr->idunion.rms.stat, 0.0, 0);
                  break;
		case GR:
		  break;
		case ORDERB:
                  printf("--Orderb %s: ", idptr->id_name);
                  printstat(&idptr->idunion.work.resunion.scalar.stat, 0.0, 0);
		  break;
		case PERS:
		  break;
		case ORIENT:
		  break;
		case BEAD:
		  break;
		case DENSITY1:
		  break;
		case CRD:
		  break;
		case MCARIG:
		  break;
		case MCA_P2P1:
		  break;
		case RIG:
		  break;
		case ACFF:
		  break;
		case CCFF:
		  break;
		case GROUP:
		  break;
		case ASA:
		  break;
		default:
		  printf("Program error: unexpected id_type %d\n", 
								idptr->id_type);
		  exit(1);
		}
	}
}


printpucker(idptr)
idstruct	*idptr;
{
	int		npuck = idptr->idunion.puck.npuckers;
	int		nresults = idptr->idunion.puck.npoints - 3;
	statstruct	*stats = idptr->idunion.puck.stats;
	int		i, j;

	printf("--Pucker %s:\n", idptr->id_name);

	for (i=0; i<npuck; i++) {
		printf("  pucker # %d:\n", i+1);
		for (j=0; j<nresults/2; j++) {
			printf("    amplitude %d: ", j+1);
			printstat(stats, 0.0, 0);
			stats++;
			printf("    angle     %d: ", j+1);
			printCstat((Cstatstruct *) stats, 0.0);
			stats++;
		}
		if (nresults % 2) {
			printf("    amplitude %d: ", j+1);
			printstat(stats, 0.0, 0);
			stats++;
		}
	}
}

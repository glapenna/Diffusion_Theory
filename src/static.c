
/************************************************************************
*                            AMBER                                     **
*                                                                      **
* Copyright (c) 1986, 1991 Regents of the University of California     **
*                       All Rights Reserved.                           ** 
*                                                                      **
*  This software provided pursuant to a license agreement containing   **
*  restrictions on its disclosure, duplication, and use. This software **
*  contains confidential and proprietary information, and may not be   **
*  extracted or distributed, in whole or in part, for any purpose      **
*  whatsoever, without the express written permission of the authors.  **
*  This notice, and the associated author list, must be attached to    **
*  all copies, or extracts, of this software. Any additional           **
*  restrictions set forth in the license agreement also apply to this  **
*  software.                                                           **
*************************************************************************/
#include "head.h"

extern int      errno;

extern idstruct		*startid(), *getstruct();
// extern char		tok[],prmtop[];
// char		tok[],prmtop[];
extern parmstruct 	*defprm;
extern char		*get();
extern FILE		*genopen();

/***********************************************************************
 							SETUPSTAT()
************************************************************************/

/* 
 * setupstat() - allocate space & read in an id's worth of STATIC 
 *	coordinate set(s)
 */

void
setupstat()
{
	idstruct	*statptr, *idptr;
	setstruct	*set, *lastset;
	FILE		*file;
	parmstruct	*prm;
	char		id[TOKLEN], tmpname[TOKLEN+3];
	_REAL		*crd;
	int		crdspace, i, ifmt, nfiles = 0, box = 0, natbox;

	lastset = NULL;

	/*
	 * get stat record in idlist & start filling in
	 */

	statptr = startid(STATIC, id);
        statptr->idunion.stat.natoms = 0;

	/*
	 *  see if parm id 
	 */

	gettok();
	if ((idptr = getstruct(tok)) != NULL) {
		if (idptr->id_type == PARM) {
			/* TODO: err if not null & not parm? */
			prm = idptr->idunion.parm.prm;
		}
	} else {
		prm = defprm;
		printf("Using default parm %s for STAT %s\n", prmtop, id);
	}

	/*
	 *  see if ATOM &/or BOX stmt, determining size of crd array(s)
	 */

	for (;;) {
		if (!strcmp(tok, "ATOM")) {
        		if (statptr->idunion.stat.natoms)
				inerr("STAT: 2 ATOM statements?", "");
        		statptr->idunion.stat.natoms = checkatoms(prm, 0);
		} else if (!strcmp(tok, "BOX")) {
			if (box)
				inerr("STAT: 2 BOX statements?", "");
			box++;
			gettok();
		} else
			break;
	}

        if (!statptr->idunion.stat.natoms)
        	statptr->idunion.stat.natoms = checkatoms(prm, 0);

	natbox = statptr->idunion.stat.natoms;
	crdspace = 3 * natbox * sizeof(_REAL);

	/*
	 *  get file names, starting w/ current token, & read coords
	 */

        statptr->idunion.stat.nsets = 0;
	for (;;) {

		/*
		 * check for done 
		 */

		if (tok[0] == ';') {
			if (!nfiles)
				inerr ("STATIC: no files named, id=", id);
			statptr->id_setup++;
			return;
		}
		nfiles++;
		
		/*
		 * check that tok isn't keyword, open file & get it ready
		 */

		if (iskeywd())
			inerr("STAT: expected file name, got keyword: ", tok);
		/* allow room for .Z extension */
		strcpy(tmpname, tok);
		file = genopen(tmpname, "stat");
		ifmt = fmtnstart(&file, tmpname, statptr->idunion.stat.natoms);

		/*
		 * read in coordinate set(s) for this file
		 */

		for(statptr->idunion.stat.nsets=0; ;
					statptr->idunion.stat.nsets++) {

			/*
			 * allocate space, queue it on stat struct
			 *	and read a set
			 */

			set = (setstruct *) get (sizeof(setstruct));
			set->crd = (_REAL *) get(crdspace);
			set->boxc = (_REAL *) get(3);
			set->max = statptr->idunion.stat.natoms;
			set->prm = prm;
			if (lastset == NULL)
				statptr->idunion.stat.sets = set;
			else
				lastset->next = set;
			lastset = set;
			set->next = NULL;
			crd = set->crd;
			for (i=0; i<natbox; i++) {
				int	x = 3 * i;
				if (scanfmt(ifmt, file, &crd[x])  ) {
					if (i == 0)
						break;
					printf("Error: %s %s %s %d %s %d\n",
					    "non-real or unexpected EOF in", 
					    tmpname, "in", 
					    statptr->idunion.stat.nsets, 
					    "set near atom", i);
					exit(1);
				}
			}
			printf("stat: read in %d atoms\n", i);
			if (ifmt == RST  ||  ifmt == PDB)
				break;
			if (i == 0)
				break;
		}
		genclose(file, tmpname);
		gettok();
	}
}


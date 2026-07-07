#include "head.h"

// extern char	tok[];
// char	tok[];

/************************************************************************
								MOD()
*************************************************************************/

void
mod(int *modx, selectstruct *select)
{
	if (*modx)
		inerr("MOD: only 1 MOD allowed","");

	gettok();
	if (!tokint(modx))
		inerr("MOD: expected a number", "");
	if (*modx < 2)
		inerr("MOD: expect a number > 1","");
	if (select != NULL) 
		*modx = *modx * -1;	/* flag to apply mod to select */
}


/************************************************************************
								SELECT()
*************************************************************************/

#define NSELECT 10

selectstruct *
selectx()
{
	selectstruct 	*sel, *asel, nulsel;
	int 		i, sn = 0, getx = 1, dash = 0;

	nulsel.item = 0;
	nulsel.code = 0;

	sel = (selectstruct *) get(NSELECT*sizeof(selectstruct));
	for (i=0; i<NSELECT; i++)
		sel[i] = nulsel;

	gettok();
	if (tok[0] != '(')
		inerr("SELECT: expected '('","");

	for (;;) {
		gettok();

		/*
		 * wrap up if necc
		 */

		if (tok[0] == ')') {
			if (sn == 0) {
				if (dash)
					inerr("SELECT: all sets specified","");
				inerr("SELECT: nothing specified","");
			}
			if (dash)
				sel[sn].code = SELTILLEND;
			else
				sel[sn].code = SELDONE;
			return(sel);
		} 

		/*
		 * grow array if necc
		 */

		if (!((sn+1)%NSELECT)) {
			getx++;
			asel = (selectstruct *) 
					get(getx*NSELECT*sizeof(selectstruct));
			for (i=0; i<sn; i++) 
				asel[i] = sel[i];
			for (i=sn; i<getx*NSELECT; i++)
				asel[i] = nulsel;
			free((char *)sel);
			sel = asel;
		}

		/*
		 * token
		 */

		if (tok[0] == '-') {
			dash = 1;
			sel[sn].code = SELTILL;	/* tolerating double dash */
		} else if (tok[0] == ',') {
			if (dash)
				inerr("SELECT: '-,' not allowed","");
		} else {		/* better be a # */
			if (!tokint(&sel[sn].item))
				inerr("SELECT: expected a number", "");
			if (sel[sn].item < 1)
				inerr("SELECT: numbers must be >= 1","");
			if (!sel[sn].code)
				sel[sn].code = SEL;
			sn++;
		}
	}
}

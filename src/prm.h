#define MAXDIHK 7
#define MAXDIHKI 4

typedef struct parm {
	int 	Nmol,Natom,Ntypes,Nbond,Nangle,Ndihe,Ndihei,Natexcl,Ncharge,
                Nbondp,Nanglep,Ndihep,Ndiheip,Nsit,Nres,
                Nsolu,Nions,Nsolv,Natom_solu,Natom_ions,Natom_solv,
                Nsolu_tot,Nions_tot,Nsolv_tot;
        int     *atomtype,*Ibond,*Iangle,*Idihe,*Idihei,*Nexcl,*Iexcl;
        int     *excltemp;
	_REAL   *Masses,*Charge,*Kbond,*Kangle,*Kdihe,*Kdihei,
                *ljAB,*fric;
	_REAL	sig0r,eps0r,mass0r;
        char    *AtomNames, *ResNames, *ResNamesa;
        int     *Swapat;
        int     Nmov,Nrest,Nener,Nmolt;
        int     *Imov_at,*Iat_mov;
        double  *Ener;
        int     *ResNum;
        int     *Ipres,*iend;
} parmstruct;

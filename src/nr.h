/* typedef struct FCOMPLEX {_REAL r,i;} fcomplex; */
typedef struct IMMENSE {unsigned long l,r;} immense;
typedef struct GREAT {unsigned short l,c,r;} great;

#define TRADITIONAL 1   /* but leave LINT_ARGS and ANSI undefined */
#ifdef LINT_ARGS
	void eigvec_nr(_REAL **, int, _REAL *, _REAL *);
	int invert(_REAL **, int, int, _REAL **);
	void adi(double **, double **, double **, double **, 
		double **, double **, double **, double **, 
		int, int, double, double, double);
	void amoeba(_REAL **, _REAL *, int, _REAL, _REAL (*)(), int *);
	void anneal(_REAL *, _REAL *, int *, int);
	void avevar(_REAL *, int, _REAL *, _REAL *);
	void balanc(_REAL **, int);
	void bcucof(_REAL *, _REAL *, _REAL *, _REAL *, _REAL, 
		_REAL, _REAL **);
	void bcuint(_REAL *, _REAL *, _REAL *, _REAL *, _REAL, 
		_REAL, _REAL, _REAL, _REAL, _REAL, _REAL *, 
		_REAL *, _REAL *);
	_REAL bessi(int, _REAL);
	_REAL bessi0(_REAL);
	_REAL bessi1(_REAL);
	_REAL bessj(int, _REAL);
	_REAL bessj0(_REAL);
	_REAL bessj1(_REAL);
	_REAL bessk(int, _REAL);
	_REAL bessk0(_REAL);
	_REAL bessk1(_REAL);
	_REAL bessy(int, _REAL);
	_REAL bessy0(_REAL);
	_REAL bessy1(_REAL);
	_REAL beta(_REAL, _REAL);
	_REAL betacf(_REAL, _REAL, _REAL);
	_REAL betai(_REAL, _REAL, _REAL);
	_REAL bico(int, int);
	void bksub(int, int, int, int, int, _REAL ***);
	_REAL bnldev(_REAL, int, int *);
	_REAL brent(_REAL, _REAL, _REAL, _REAL (*)(), _REAL, _REAL *);
	void bsstep(_REAL *, _REAL *, int, _REAL *, _REAL,
		_REAL, _REAL *, _REAL *, _REAL *, void (*)());
	void caldat(long, int *, int *, int *);
	_REAL cel(_REAL, _REAL, _REAL, _REAL);
	void chder(_REAL, _REAL, _REAL *, _REAL *, int);
	_REAL chebev(_REAL, _REAL, _REAL *, int, _REAL);
	void chebft(_REAL, _REAL, _REAL *, int, _REAL (*)());
	void chebpc(_REAL *, _REAL *, int);
	void chint(_REAL, _REAL, _REAL *, _REAL *, int);
	void chsone(_REAL *, _REAL *, int, int, _REAL *, 
		_REAL *, _REAL *);
	void chstwo(_REAL *, _REAL *, int, int, _REAL *, 
		_REAL *, _REAL *);
	void cntab1(int **, int, int, _REAL *, _REAL *, 
		_REAL *, _REAL *, _REAL *);
	void cntab2(int **, int, int, _REAL *, _REAL *, 
		_REAL *, _REAL *, _REAL *, _REAL *, _REAL *, _REAL *);
	void convlv(_REAL *, int, _REAL *, int, int, _REAL *);
	void correl(_REAL *, _REAL *, int, _REAL *);
	void cosft(_REAL *, int, int);
	void covsrt(_REAL **, int, int *, int);
	void crank(int, _REAL *, _REAL *);
	_REAL dbrent(_REAL, _REAL, _REAL, _REAL (*)(), _REAL (*)(),
		_REAL, _REAL *);
	void ddpoly(_REAL *, int, _REAL, _REAL *, int);
	void des(immense, immense, int *, int, immense *);
	void ks(immense, int, great *);
	void cyfun(unsigned long, great, unsigned long *);
	_REAL df1dim(_REAL);
	void dfpmin(_REAL *, int, _REAL, int *, _REAL *, _REAL (*)(),
		void (*)());
	void difeq(int, int, int, int, int, int, int *, int, 
		_REAL **, _REAL **);
	void dlinmin(_REAL *, _REAL *, int, _REAL *, _REAL (*)(), void (*)());
	void eclass(int *, int, int *, int *, int);
	void eclazz(int *, int, int (*)());
	void eigsrt(_REAL *, _REAL **, int);
	_REAL el2(_REAL, _REAL, _REAL, _REAL);
	void elmhes(_REAL **, int);
	_REAL erf(_REAL);
	_REAL erfc(_REAL);
	_REAL erfcc(_REAL);
	void eulsum(_REAL *, _REAL, int, _REAL *);
	_REAL evlmem(_REAL, _REAL *, int, _REAL);
	_REAL expdev(int *);
	_REAL f1dim(_REAL);
	_REAL factln(int);
	_REAL factrl(int);
	void fgauss(_REAL, _REAL *, _REAL *, _REAL *, int);
	void fit(_REAL *, _REAL *, int, _REAL *, int, _REAL *, 
		_REAL *, _REAL *, _REAL *, _REAL *, _REAL *);
	void fixrts(_REAL *, int);
	void fleg(_REAL, _REAL *, int);
	void flmoon(int, int, long *, _REAL *);
	void four1(_REAL *, int, int);
	void fourn(_REAL *, int *, int, int);
	void fpoly(_REAL, _REAL *, int);
	void frprmn(_REAL *, int, _REAL, int *, _REAL *, _REAL (*)(),
		void (*)());
	void ftest(_REAL *, int, _REAL *, int, _REAL *, _REAL *);
	_REAL gamdev(int, int *);
	_REAL gammln(_REAL);
	_REAL gammp(_REAL, _REAL);
	_REAL gammq(_REAL, _REAL);
	_REAL gasdev(int *);
	void gauleg(double, double, double *, double *, int);
	void gaussj(_REAL **, int, _REAL **, int);
	void gcf(_REAL *, _REAL, _REAL, _REAL *);
	_REAL golden(_REAL, _REAL, _REAL, _REAL (*)(), _REAL, _REAL *);
	void gser(_REAL *, _REAL, _REAL, _REAL *);
	void hqr(_REAL **, int, _REAL *, _REAL *);
	void hunt(_REAL *, int, _REAL, int *);
	void indexx(int, _REAL *, int *);
	int irbit1(unsigned long int *);
	int irbit2(unsigned long int *);
	void jacobi(_REAL **, int, _REAL *, _REAL **, int *);
	long julday(int, int, int);
	void kendl1(_REAL *, _REAL *, int, _REAL *, _REAL *, _REAL *);
	void kendl2(_REAL **, int, int, _REAL *, _REAL *, _REAL *);
	void ksone(_REAL *, int, _REAL (*)(), _REAL *, _REAL *);
	void kstwo(_REAL *, int, _REAL *, int, _REAL *, _REAL *);
	void laguer(fcomplex *, int, fcomplex *, _REAL, int);
	void lfit(_REAL *, _REAL *, _REAL *, int, _REAL *, int, 
		int *, int, _REAL **, _REAL *, void (*)());
	void linmin(_REAL *, _REAL *, int, _REAL *, _REAL (*)());
	void locate(_REAL *, int, _REAL, int *);
	void lubksb(_REAL **, int, int *, _REAL *);
	void ludcmp(_REAL **, int, int *, _REAL *);
	void mdian1(_REAL *, int, _REAL *);
	void mdian2(_REAL *, int, _REAL *);
	void medfit(_REAL *, _REAL *, int, _REAL *, _REAL *, _REAL *);
	void memcof(_REAL *, int, int, _REAL *, _REAL *);
	_REAL midexp(_REAL (*)(), _REAL, _REAL, int);
	_REAL midinf(_REAL (*)(), _REAL, _REAL, int);
	_REAL midpnt(_REAL (*)(), _REAL, _REAL, int);
	_REAL midsql(_REAL (*)(), _REAL, _REAL, int);
	_REAL midsqu(_REAL (*)(), _REAL, _REAL, int);
	void mmid(_REAL *, _REAL *, int, _REAL, _REAL, int, _REAL *,
		void (*)());
	void mnbrak(_REAL *, _REAL *, _REAL *, _REAL *, 
		_REAL *, _REAL *, _REAL (*)());
	void mnewt(int, _REAL *, int, _REAL, _REAL);
	void moment(_REAL *, int, _REAL *, _REAL *, _REAL *, 
		_REAL *, _REAL *, _REAL *);
	void mprove(_REAL **, _REAL **, int, int *, _REAL *, _REAL *);
	void mrqcof(_REAL *, _REAL *, _REAL *, int, _REAL *, int, 
		int *, int, _REAL **, _REAL *, _REAL *, void (*)());
	void mrqmin(_REAL *, _REAL *, _REAL *, int, _REAL *, int, 
		int *, int, _REAL **, _REAL **, _REAL *, void (*)(),
		_REAL *);
	void odeint(_REAL *, int, _REAL, _REAL, _REAL, _REAL, 
		_REAL, int *, int *, void (*)(), void (*)());
	void pcshft(_REAL, _REAL, _REAL *, int);
	void pearsn(_REAL *, _REAL *, int, _REAL *, _REAL *, _REAL *);
	void piksr2(int, _REAL *, _REAL *);
	void piksrt(int, _REAL *);
	void pinvs(int, int, int, int, int, int, _REAL ***, _REAL **);
	_REAL plgndr(int, int, _REAL);
	_REAL poidev(_REAL, int *);
	void polcoe(_REAL *, _REAL *, int, _REAL *);
	void polcof(_REAL *, _REAL *, int, _REAL *);
	void poldiv(_REAL *, int, _REAL *, int, _REAL *, _REAL *);
	void polin2(_REAL *, _REAL *, _REAL **, int, int, _REAL, 
		_REAL, _REAL *, _REAL *);
	void polint(_REAL *, _REAL *, int, _REAL, _REAL *, _REAL *);
	void powell(_REAL *, _REAL **, int, _REAL, int *, _REAL *, _REAL (*)());
	void predic(_REAL *, int, _REAL *, int, _REAL *, int);
	_REAL probks(_REAL);
	void pzextr(int, _REAL, _REAL *, _REAL *, _REAL *, int, int);
	void qcksrt(int, _REAL *);
	_REAL qgaus(_REAL (*)(_REAL), _REAL, _REAL);
	_REAL qromb(_REAL (*)(_REAL), _REAL, _REAL);
	_REAL qromo(_REAL (*)(_REAL), _REAL, _REAL, _REAL (*)());
	void qroot(_REAL *, int, _REAL *, _REAL *, _REAL);
	_REAL qsimp(_REAL (*)(_REAL), _REAL, _REAL);
	_REAL qtrap(_REAL (*)(_REAL), _REAL, _REAL);
	_REAL quad3d(_REAL (*)(_REAL), _REAL, _REAL);
	_REAL ran0(int *);
	_REAL ran1(int *);
	_REAL ran2(long *);
	_REAL ran3(int *);
	_REAL ran4(int *);
	void rank(int, int *, int *);
	void ratint(_REAL *, _REAL *, int, _REAL, _REAL *, _REAL *);
	void realft(_REAL *, int, int);
	void red(int, int, int, int, int, int, int, int, int, int, 
		int, _REAL ***, _REAL **);
	void rk4(_REAL *, _REAL *, int, _REAL, _REAL, _REAL *, void (*)());
	void rkdumb(_REAL *, int, _REAL, _REAL, int, void (*)());
	void rkqc(_REAL *, _REAL *, int, _REAL *, _REAL, _REAL, 
		_REAL *, _REAL *, _REAL *, void (*)());
	_REAL rofunc(_REAL);
	_REAL rtbis(_REAL (*)(), _REAL, _REAL, _REAL);
	_REAL rtflsp(_REAL (*)(), _REAL, _REAL, _REAL);
	_REAL rtnewt(void (*)(), _REAL, _REAL, _REAL);
	_REAL rtsafe(void (*)(), _REAL, _REAL, _REAL);
	_REAL rtsec(_REAL (*)(), _REAL, _REAL, _REAL);
	void rzextr(int, _REAL, _REAL *, _REAL *, _REAL *, int, int);
	void scrsho(_REAL (*)());
	void shell(int, _REAL *);
	void shoot(int, _REAL *, _REAL *, int, _REAL, _REAL, _REAL, 
		_REAL, _REAL, _REAL *, _REAL *);
	void shootf(int, _REAL *, _REAL *, _REAL *, _REAL *, int, int, 
		_REAL, _REAL, _REAL, _REAL, _REAL, _REAL, _REAL *, 
		_REAL *, _REAL *);
	void simp1(_REAL **, int, int *, int, int, int *, _REAL *);
	void simp2(_REAL **, int, int *, int, int *, int, _REAL *);
	void simp3(_REAL **, int, int, int, int);
	void simplx(_REAL **, int, int, int, int, int, 
		int *, int *, int *);
	void sinft(_REAL *, int);
	void smooft(_REAL *, int, _REAL);
	void sncndn(_REAL, _REAL, _REAL *, _REAL *, _REAL *);
	void solvde(int, _REAL, _REAL, _REAL *, int *, int, int, 
		int, _REAL **, _REAL ***, _REAL **);
	void sor(double **, double **, double **, double **, double **, 
		double **, double **, int, double);
	void sort(int, _REAL *);
	void sort2(int, _REAL *, _REAL *);
	void sort3(int, _REAL *, _REAL *, _REAL *);
	void sparse(_REAL *, int, _REAL *, _REAL *);
	void spctrm(FILE *, _REAL *, int, int, int);
	void spear(_REAL *, _REAL *, int, _REAL *, 
		_REAL *, _REAL *, _REAL *, _REAL *);
	void splie2(_REAL *, _REAL *, _REAL **, int, int, _REAL **);
	void splin2(_REAL *, _REAL *, _REAL **, _REAL **, int, int, 
		_REAL, _REAL, _REAL *);
	void spline(_REAL *, _REAL *, int, _REAL, _REAL, _REAL *);
	void splint(_REAL *, _REAL *, _REAL *, int, _REAL, _REAL *);
	void svbksb(_REAL **, _REAL *, _REAL **, int, int,
		_REAL *, _REAL *);
	void svdcmp(_REAL **, int, int, _REAL *, _REAL **);
	void svdfit(_REAL *, _REAL *, _REAL *, int, _REAL *, int, 
		_REAL **, _REAL **, _REAL *, _REAL *, void (*)());
	void svdvar(_REAL **, int, _REAL *, _REAL **);
	void toeplz(_REAL *, _REAL *, _REAL *, int);
	void tptest(_REAL *, _REAL *, int, _REAL *, _REAL *);
	void tqli(_REAL *, _REAL *, int, _REAL **);
	_REAL trapzd(_REAL (*)(_REAL), _REAL, _REAL, int);
	void tred2(_REAL **, int, _REAL *, _REAL *);
	void tridag(_REAL *, _REAL *, _REAL *, _REAL *, _REAL *, int);
	void ttest(_REAL *, int, _REAL *, int, _REAL *, _REAL *);
	void tutest(_REAL *, int, _REAL *, int, _REAL *, _REAL *);
	void twofft(_REAL *, _REAL *, _REAL *, _REAL *, int);
	void vander(_REAL *, _REAL *, _REAL *, int);
	int zbrac(_REAL (*)(), _REAL *, _REAL *);
	void zbrak(_REAL (*)(), _REAL, _REAL, int, _REAL *, _REAL *, int *);
	_REAL zbrent(_REAL (*)(), _REAL, _REAL, _REAL);
	void zroots(fcomplex *, int, fcomplex *, int);
#endif
#ifdef ANSI
	void eigvec_nr(_REAL **a, int n, _REAL *wr, _REAL *z);
	int invert(_REAL **a, int m, int n, _REAL **b);
	void  adi(double **a, double **b, double **c, double **d, double **e,
		double **f, double **g, double **u, int jmax, int k,
		double alpha, double beta, double eps);
	void  amoeba(_REAL **p, _REAL *y, int ndim, _REAL ftol,
		_REAL (*funk)(_REAL *), int *iter);
	void  anneal(_REAL *x, _REAL *y, int *iorder, int ncity);
	void  avevar(_REAL *data, int n, _REAL *ave, _REAL *svar);
	void  balanc(_REAL **a, int n);
	void  bcucof(_REAL *y, _REAL *y1, _REAL *y2, _REAL *y12, _REAL d1,
		_REAL d2, _REAL **c);
	void  bcuint(_REAL *y, _REAL *y1, _REAL *y2, _REAL *y12, _REAL x1l, 
		_REAL x1u, _REAL x2l, _REAL x2u, _REAL x1, _REAL x2, 
		_REAL *ansy, _REAL *ansy1, _REAL *ansy2);
	_REAL bessi(int n, _REAL x);
	_REAL bessi0(_REAL x);
	_REAL bessi1(_REAL x);
	_REAL bessj(int n, _REAL x);
	_REAL bessj0(_REAL x);
	_REAL bessj1(_REAL x);
	_REAL bessk(int n, _REAL x);
	_REAL bessk0(_REAL x);
	_REAL bessk1(_REAL x);
	_REAL bessy(int n, _REAL x);
	_REAL bessy0(_REAL x);
	_REAL bessy1(_REAL x);
	_REAL beta(_REAL z, _REAL w);
	_REAL betacf(_REAL a, _REAL b, _REAL x);
	_REAL betai(_REAL a, _REAL b, _REAL x);
	_REAL bico(int n, int k);
	void  bksub(int ne, int nb, int jf, int k1, int k2, _REAL ***c);
	_REAL bnldev(_REAL pp, int n, int *idum);
	_REAL brent(_REAL ax, _REAL bx, _REAL cx, _REAL (*f)(_REAL), _REAL tol,
		_REAL *xmin);
	void  bsstep(_REAL *y, _REAL *dydx, int nv, _REAL *xx, _REAL htry,
		_REAL eps, _REAL *yscal, _REAL *hdid, _REAL *hnext, 
		void (*derivs)(_REAL,_REAL *,_REAL *));
	void  caldat(long julian, int *mm, int *id, int *iyyy);
	_REAL cel(_REAL qqc, _REAL pp, _REAL aa, _REAL bb);
	void  chder(_REAL a, _REAL b, _REAL *c, _REAL *cder, int n);
	_REAL chebev(_REAL a, _REAL b, _REAL *c, int m, _REAL x);
	void  chebft(_REAL a, _REAL b, _REAL *c, int n, _REAL (*func)(_REAL));
	void  chebpc(_REAL *c, _REAL *d, int n);
	void  chint(_REAL a, _REAL b, _REAL *c, _REAL *cint, int n);
	void  chsone(_REAL *bins, _REAL *ebins, int nbins, int knstrn, 
		_REAL *df, _REAL *chsq, _REAL *prob);
	void  chstwo(_REAL *bins1, _REAL *bins2, int nbins, int knstrn, 
		_REAL *df, _REAL *chsq, _REAL *prob);
	void  cntab1(int **nn, int n1, int nj, _REAL *chisq, _REAL *df, 
		_REAL *prob, _REAL *cramrv, _REAL *ccc);
	void  cntab2(int **nn, int ni, int nj, _REAL *h, _REAL *hx, _REAL *hy, 
		_REAL *hygx, _REAL *hxgy, _REAL *uygx, _REAL *uxgy,
		_REAL *uxy);
	void  convlv(_REAL *data, int n, _REAL *respns, int m, int isign, 
		_REAL *ans);
	void  correl(_REAL *data1, _REAL *data2, int n, _REAL *ans);
	void  cosft(_REAL *y, int n, int isign);
	void  covsrt(_REAL **covar, int ma, int *lista, int mfit);
	void  crank(int n, _REAL *w, _REAL *s);
	_REAL dbrent(_REAL ax, _REAL bx, _REAL cx, _REAL (*f)(_REAL),
		_REAL (*df)(_REAL), _REAL tol, _REAL *xmin);
	void  ddpoly(_REAL *c, int nc, _REAL x, _REAL *pd, int nd);
	void  des(immense inp, immense key, int *newkey, int isw, immense *out);
	void  ks(immense key, int n, great *kn);
	void  cyfun(unsigned long ir, great k, unsigned long *iout);
	_REAL df1dim(_REAL x);
	void  dfpmin(_REAL *p, int n, _REAL ftol, int *iter, _REAL *fret, 
		_REAL (*func)(_REAL *), void (*dfunc)(_REAL *,_REAL *));
	void  difeq(int k, int k1, int k2, int jsf, int is1, int isf, 
		int *indexv, int ne, _REAL **s, _REAL **y);
	void  dlinmin(_REAL *p, _REAL *xi, int n, _REAL *fret,
		_REAL (*func)(_REAL *), void (*dfunc)(_REAL *,_REAL *));
	void  eclass(int *nf, int n, int *lista, int *listb, int m);
	void  eclazz(int *nf, int n, int (*equiv)(int,int));
	void  eigsrt(_REAL *d, _REAL **v, int n);
	_REAL el2(_REAL x, _REAL qqc, _REAL aa, _REAL bb);
	void  elmhes(_REAL **a, int n);
	_REAL erf(_REAL x);
	_REAL erfc(_REAL x);
	_REAL erfcc(_REAL x);
	void  eulsum(_REAL *sum, _REAL term, int jterm, _REAL *wksp);
	_REAL evlmem(_REAL fdt, _REAL *cof, int m, _REAL pm);
	_REAL expdev(int *idum);
	_REAL f1dim(_REAL x);
	_REAL factln(int n);
	_REAL factrl(int n);
	void  fgauss(_REAL x, _REAL *a, _REAL *y, _REAL *dyda, int na);
	void  fit(_REAL *x, _REAL *y, int ndata, _REAL *sig, int mwt, _REAL *a, 
		_REAL *b, _REAL *siga, _REAL *sigb, _REAL *chi2, _REAL *q);
	void  fixrts(_REAL *d, int npoles);
	void  fleg(_REAL x, _REAL *pl, int nl);
	void  flmoon(int n, int nph, long *jd, _REAL *frac);
	void  four1(_REAL *data, int nn, int isign);
	void  fourn(_REAL *data, int *nn, int ndim, int isign);
	void  fpoly(_REAL x, _REAL *p, int np);
	void  frprmn(_REAL *p, int n, _REAL ftol, int *iter, _REAL *fret, 
		_REAL (*func)(_REAL *), void (*dfunc)(_REAL *,_REAL *));
	void  ftest(_REAL *data1, int n1, _REAL *data2, int n2, _REAL *f, 
		_REAL *prob);
	_REAL gamdev(int ia, int *idum);
	_REAL gammln(_REAL xx);
	_REAL gammp(_REAL a, _REAL x);
	_REAL gammq(_REAL a, _REAL x);
	_REAL gasdev(int *idum);
	void  gauleg(double x1, double x2, double *x, double *w, int n);
	void  gaussj(_REAL **a, int n, _REAL **b, int m);
	void  gcf(_REAL *gammcf, _REAL a, _REAL x, _REAL *gln);
	_REAL golden(_REAL ax, _REAL bx, _REAL cx, _REAL (*f)(_REAL), _REAL tol, 
		_REAL *xmin);
	void  gser(_REAL *gamser, _REAL a, _REAL x, _REAL *gln);
	void  hqr(_REAL **a, int n, _REAL *wr, _REAL *wi);
	void  hunt(_REAL *xx, int n, _REAL x, int *jlo);
	void  indexx(int n, _REAL *arrin, int *indx);
	int   irbit1(unsigned long int *iseed);
	int   irbit2(unsigned long int *iseed);
	void  jacobi(_REAL **a, int n, _REAL *d, _REAL **v, int *nrot);
	long  julday(int mm, int id, int iyyy);
	void  kendl1(_REAL *data1, _REAL *data2, int n, _REAL *tau, _REAL *z,
		_REAL *prob);
	void  kendl2(_REAL **tab, int i, int j, _REAL *tau, _REAL *z,
		_REAL *prob);
	void  ksone(_REAL *data, int n, _REAL (*func)(_REAL), _REAL *d,
		_REAL *prob);
	void  kstwo(_REAL *data1, int n1, _REAL *data2, int n2, _REAL *d,
		_REAL *prob);
	void  laguer(fcomplex *a, int m, fcomplex *x, _REAL eps, int polish);
	void  lfit(_REAL *x, _REAL *y, _REAL *sig, int ndata, _REAL *a, int ma, 
		int *lista, int mfit, _REAL **covar, _REAL *chisq,
		void (*funcs)(_REAL,_REAL *,int));
	void  linmin(_REAL *p, _REAL *xi, int n, _REAL *fret, _REAL (*func)(_REAL));
	void  locate(_REAL *xx, int n, _REAL x, int *j);
	void  lubksb(_REAL **a, int n, int *indx, _REAL *b);
	void  ludcmp(_REAL **a, int n, int *indx, _REAL *d);
	void  mdian1(_REAL *x, int n, _REAL *xmed);
	void  mdian2(_REAL *x, int n, _REAL *xmed);
	void  medfit(_REAL *x, _REAL *y, int ndata, _REAL *a, _REAL *b,
		_REAL *abdev);
	void  memcof(_REAL *data, int n, int m, _REAL *pm, _REAL *cof);
	_REAL midexp(_REAL (*funk)(_REAL), _REAL aa, _REAL bb, int n);
	_REAL midinf(_REAL (*funk)(_REAL), _REAL aa, _REAL bb, int n);
	_REAL midpnt(_REAL (*func)(_REAL), _REAL a, _REAL b, int n);
	_REAL midsql(_REAL (*funk)(_REAL), _REAL aa, _REAL bb, int n);
	_REAL midsqu(_REAL (*funk)(_REAL), _REAL aa, _REAL bb, int n);
	void  mmid(_REAL *y, _REAL *dydx, int nvar, _REAL xs, _REAL htot,
		int nstep, _REAL *yout,
		void (*derivs)(_REAL,_REAL *,_REAL *));
	void  mnbrak(_REAL *ax, _REAL *bx, _REAL *cx, _REAL *fa, _REAL *fb, 
		_REAL *fc, _REAL (*func)(_REAL));
	void  mnewt(int ntrial, _REAL *x, int n, _REAL tolx, _REAL tolf);
	void  moment(_REAL *data, int n, _REAL *ave, _REAL *adev, _REAL *sdev, 
		_REAL *svar, _REAL *skew, _REAL *curt);
	void  mprove(_REAL **a, _REAL **alud, int n, int *indx, _REAL *b, 
		_REAL *x);
	void  mrqcof(_REAL *x, _REAL *y, _REAL *sig, int ndata, _REAL *a, int ma, 
		int *lista, int mfit, _REAL **alpha, _REAL *beta, _REAL
		*chisq, void (*funcs)(_REAL,_REAL *,_REAL *,_REAL *,int));
	void  mrqmin(_REAL *x, _REAL *y, _REAL *sig, int ndata, _REAL *a,
		int ma, int *lista, int mfit, _REAL **covar, _REAL **alpha, 
		_REAL *chisq, void (*funcs)(_REAL,_REAL *,_REAL *,_REAL *,
		int),_REAL *alamda);
	void  odeint(_REAL *ystart, int nvar, _REAL x1, _REAL x2, _REAL eps,
		_REAL h1, _REAL hmin, int *nok, int *nbad,
		void (*derivs)(_REAL,_REAL *,_REAL *),
		void  (*rkqc)(_REAL *,_REAL *,int,_REAL *,_REAL,_REAL,_REAL
		*,_REAL *,_REAL *,void (*)(_REAL,_REAL *,_REAL *)));
	void  pcshft(_REAL a, _REAL b, _REAL *d, int n);
	void  pearsn(_REAL *x, _REAL *y, int n, _REAL *r, _REAL *prob, _REAL *z);
	void  piksr2(int n, _REAL *arr, _REAL *brr);
	void  piksrt(int n, _REAL *arr);
	void  pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
		_REAL ***c, _REAL **s);
	_REAL plgndr(int l, int m, _REAL x);
	_REAL poidev(_REAL xm, int *idum);
	void  polcoe(_REAL *x, _REAL *y, int n, _REAL *cof);
	void  polcof(_REAL *xa, _REAL *ya, int n, _REAL *cof);
	void  poldiv(_REAL *u, int n, _REAL *v, int nv, _REAL *q, _REAL *r);
	void  polin2(_REAL *x1a, _REAL *x2a, _REAL **ya, int m, int n, _REAL x1, 
		_REAL x2, _REAL *y, _REAL *dy);
	void  polint(_REAL *xa, _REAL *ya, int n, _REAL x, _REAL *y, _REAL *dy);
	void  powell(_REAL *p, _REAL **xi, int n, _REAL ftol, int *iter,
		_REAL *fret, _REAL (*func)(_REAL *));
	void  predic(_REAL *data, int ndata, _REAL *d, int npoles, 
		_REAL *future, int nfut);
	_REAL probks(_REAL alam);
	void  pzextr(int iest, _REAL xest, _REAL *yest, _REAL *yz, _REAL *dy,
		int nv, int nuse);
	void  qcksrt(int n, _REAL *arr);
	_REAL qgaus(_REAL (*func)(_REAL), _REAL a, _REAL b);
	_REAL qromb(_REAL (*func)(_REAL), _REAL a, _REAL b);
	_REAL qromo(_REAL (*func)(_REAL), _REAL a, _REAL b,
		_REAL (*choose)(_REAL (*)(_REAL),_REAL,_REAL,int));
	void  qroot(_REAL *p, int n, _REAL *b, _REAL *c, _REAL eps);
	_REAL qsimp(_REAL (*func)(_REAL), _REAL a, _REAL b);
	_REAL qtrap(_REAL (*func)(_REAL), _REAL a, _REAL b);
	_REAL quad3d(_REAL (*func)(_REAL,_REAL,_REAL), _REAL x1, _REAL x2);
	_REAL ran0(int *idum);
	_REAL ran1(int *idum);
	_REAL ran2(long *idum);
	_REAL ran3(int *idum);
	_REAL ran4(int *idum);
	void  rank(int n, int *indx, int *irank);
	void  ratint(_REAL *xa, _REAL *ya, int n, _REAL x, _REAL *y, _REAL *dy);
	void  realft(_REAL *data, int n, int isign);
	void  red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
		int ic1, int jc1, int jcf, int kc, _REAL ***c, _REAL **s);
	void  rk4(_REAL *y, _REAL *dydx, int n, _REAL x, _REAL h, _REAL *yout,
		void (*derivs)(_REAL,_REAL *,_REAL *));
	void  rkdumb(_REAL *vstart, int nvar, _REAL x1, _REAL x2, int nstep, 
		void (*derivs)(_REAL,_REAL *,_REAL *));
	void  rkqc(_REAL *y, _REAL *dydx, int n, _REAL *x, _REAL htry, 
		_REAL eps, _REAL *yscal, _REAL *hdid, _REAL *hnext, 
		void (*derivs)(_REAL,_REAL *,_REAL *));
	_REAL rofunc(_REAL b);
	_REAL rtbis(_REAL (*func)(_REAL), _REAL x1, _REAL x2, _REAL xacc);
	_REAL rtflsp(_REAL (*func)(_REAL), _REAL x1, _REAL x2, _REAL xacc);
	_REAL rtnewt(void (*funcd)(_REAL,_REAL *,_REAL *), _REAL x1, _REAL x2,
		_REAL xacc);
	_REAL rtsafe(void (*funcd)(_REAL,_REAL *,_REAL *), _REAL x1, _REAL x2,
		_REAL xacc);
	_REAL rtsec(_REAL (*func)(_REAL), _REAL x1, _REAL x2, _REAL xacc);
	void  rzextr(int iest, _REAL xest, _REAL *yest, _REAL *yz, _REAL *dy,
		int nv, int nuse);
	void  scrsho(_REAL (*fx)(_REAL));
	void  shell(int n, _REAL *arr);
	void  shoot(int nvar, _REAL *v, _REAL *delv, int n2, _REAL x1, _REAL x2,
		_REAL eps, _REAL h1, _REAL hmin, _REAL *f, _REAL *dv);
	void  shootf(int nvar, _REAL *v1, _REAL *v2, _REAL *delv1, _REAL *delv2,
		int n1, int n2, _REAL x1, _REAL x2, _REAL xf, _REAL eps, 
		_REAL h1, _REAL hmin, _REAL *f, _REAL *dv1, _REAL *dv2);
	void  simp1(_REAL **a, int mm, int *ll, int nll, int iabf, int *kp,
		_REAL *bmax);
	void  simp2(_REAL **a, int n, int *l2, int nl2, int *ip, int kp,
		_REAL *q1);
	void  simp3(_REAL **a, int i1, int k1, int ip, int kp);
	void  simplx(_REAL **a, int m, int n, int m1, int m2, int m3, 
		int *icase, int *izrov, int *iposv);
	void  sinft(_REAL *y, int n);
	void  smooft(_REAL *y, int n, _REAL pts);
	void  sncndn(_REAL uu, _REAL emmc, _REAL *sn, _REAL *cn, _REAL *dn);
	void  solvde(int itmax, _REAL conv, _REAL slowc, _REAL *scalv,
		int *indexv, int ne, int nb, int m, _REAL **y, _REAL ***c,
		_REAL **s);
	void  sor(double **a, double **b, double **c, double **d, double **e, 
		double **f, double **u, int jmax, double rjac);
	void  sort(int n, _REAL *ra);
	void  sort2(int n, _REAL *ra, _REAL *rb);
	void  sort3(int n, _REAL *ra, _REAL *rb, _REAL *rc);
	void  sparse(_REAL *b, int n, _REAL *x, _REAL *rsq);
	void  spctrm(FILE *fp, _REAL *p, int m, int k, int ovrlap);
	void  spear(_REAL *data1, _REAL *data2, int n, _REAL *d, _REAL *zd,
		_REAL *probd, _REAL *rs, _REAL *probrs);
	void  splie2(_REAL *x1a, _REAL *x2a, _REAL **ya, int m, int n,
		_REAL **y2a);
	void  splin2(_REAL *x1a, _REAL *x2a, _REAL **ya, _REAL **y2a, int m,
		int n, _REAL x1, _REAL x2, _REAL *y);
	void  spline(_REAL *x, _REAL *y, int n, _REAL yp1, _REAL ypn, _REAL *y2);
	void  splint(_REAL *xa, _REAL *ya, _REAL *y2a, int n, _REAL x, _REAL *y);
	void  svbksb(_REAL **u, _REAL *w, _REAL **v, int m, int n, _REAL *b,
		_REAL *x);
	void  svdcmp(_REAL **a, int m, int n, _REAL *w, _REAL **v);
	void  svdfit(_REAL *x, _REAL *y, _REAL *sig, int ndata, _REAL *a, 
		int ma, _REAL **u, _REAL **v, _REAL *w, _REAL *chisq,
		void (*funcs)(_REAL,_REAL *,int));
	void  svdvar(_REAL **v, int ma, _REAL *w, _REAL **cvm);
	void  toeplz(_REAL *r, _REAL *x, _REAL *y, int n);
	void  tptest(_REAL *data1, _REAL *data2, int n, _REAL *t, _REAL *prob);
	void  tqli(_REAL *d, _REAL *e, int n, _REAL **z);
	_REAL trapzd(_REAL (*func)(_REAL), _REAL a, _REAL b, int n);
	void  tred2(_REAL **a, int n, _REAL *d, _REAL *e);
	void  tridag(_REAL *a, _REAL *b, _REAL *c, _REAL *r, _REAL *u, int n);
	void  ttest(_REAL *data1, int n1, _REAL *data2, int n2, _REAL *t,
		_REAL *prob);
	void  tutest(_REAL *data1, int n1, _REAL *data2, int n2, _REAL *t,
		_REAL *prob);
	void  twofft(_REAL *data1, _REAL *data2, _REAL *fft1, _REAL *fft2,
		int n);
	void  vander(_REAL *x, _REAL *w, _REAL *q, int n);
	int   zbrac(_REAL (*func)(_REAL), _REAL *x1, _REAL *x2);
	void  zbrak(_REAL (*fx)(_REAL), _REAL x1, _REAL x2, int n, _REAL *xb1,
		_REAL *xb2, int *nb);
	_REAL zbrent(_REAL (*func)(_REAL), _REAL x1, _REAL x2, _REAL tol);
	void  zroots(fcomplex *a, int m, fcomplex *roots, int polish);
#endif
#ifdef TRADITIONAL
	void eigvec_nr();
	int invert();
	void adi();
	void amoeba();
	void anneal();
	void avevar();
	void balanc();
	void bcucof();
	void bcuint();
	_REAL bessi();
	_REAL bessi0();
	_REAL bessi1();
	_REAL bessj();
	_REAL bessj0();
	_REAL bessj1();
	_REAL bessk();
	_REAL bessk0();
	_REAL bessk1();
	_REAL bessy();
	_REAL bessy0();
	_REAL bessy1();
	_REAL beta();
	_REAL betacf();
	_REAL betai();
	_REAL bico();
	void bksub();
	_REAL bnldev();
	_REAL brent();
	void bsstep();
	void caldat();
	_REAL cel();
	void chder();
	_REAL chebev();
	void chebft();
	void chebpc();
	void chint();
	void chsone();
	void chstwo();
	void cntab1();
	void cntab2();
	void convlv();
	void correl();
	void cosft();
	void covsrt();
	void crank();
	_REAL dbrent();
	void ddpoly();
	void des();
	void ks();
	void cyfun();
	_REAL df1dim();
	void dfpmin();
	void dlinmin();
	void difeq();
	void eclass();
	void eclazz();
	void eigsrt();
	_REAL el2();
	void elmhes();
	_REAL erf();
	_REAL erfc();
	_REAL erfcc();
	void eulsum();
	_REAL evlmem();
	_REAL expdev();
	_REAL f1dim();
	_REAL factln();
	_REAL factrl();
	void fgauss();
	void fit();
	void fixrts();
	void fleg();
	void flmoon();
	void four1();
	void fourn();
	void fpoly();
	void frprmn();
	void ftest();
	_REAL gamdev();
	_REAL gammln();
	_REAL gammp();
	_REAL gammq();
	_REAL gasdev();
	void gauleg();
	void gaussj();
	void gcf();
	_REAL golden();
	void gser();
	void hqr();
	void hunt();
	void indexx();
	int irbit1();
	int irbit2();
	void jacobi();
	long julday();
	void kendl1();
	void kendl2();
	void ksone();
	void kstwo();
	void laguer();
	void lfit();
	void linmin();
	void locate();
	void lubksb();
	void ludcmp();
	void mdian1();
	void mdian2();
	void medfit();
	void memcof();
	_REAL midexp();
	_REAL midinf();
	_REAL midpnt();
	_REAL midsql();
	_REAL midsqu();
	void mmid();
	void mnbrak();
	void mnewt();
	void moment();
	void mprove();
	void mrqcof();
	void mrqmin();
	void odeint();
	void pcshft();
	void pearsn();
	void piksr2();
	void piksrt();
	void pinvs();
	_REAL plgndr();
	_REAL poidev();
	void polcoe();
	void polcof();
	void poldiv();
	void polin2();
	void polint();
	void powell();
	void predic();
	_REAL probks();
	void pzextr();
	void qcksrt();
	_REAL qgaus();
	_REAL qromb();
	_REAL qromo();
	void qroot();
	_REAL qsimp();
	_REAL qtrap();
	_REAL quad3d();
	_REAL ran0();
	_REAL ran1();
	_REAL ran2();
	_REAL ran3();
	_REAL ran4();
	void rank();
	void ratint();
	void realft();
	void red();
	void rk4();
	void rkdumb();
	void rkqc();
	_REAL rofunc();
	_REAL rtbis();
	_REAL rtflsp();
	_REAL rtnewt();
	_REAL rtsafe();
	_REAL rtsec();
	void rzextr();
	void scrsho();
	void shell();
	void shoot();
	void shootf();
	void simp1();
	void simp2();
	void simp3();
	void simplx();
	void sinft();
	void smooft();
	void sncndn();
	void solvde();
	void sor();
	void sort();
	void sort2();
	void sort3();
	void sparse();
	void spctrm();
	void spear();
	void splie2();
	void splin2();
	void spline();
	void splint();
	void svbksb();
	void svdcmp();
	void svdfit();
	void svdvar();
	void toeplz();
	void tptest();
	void tqli();
	_REAL trapzd();
	void tred2();
	void tridag();
	void ttest();
	void tutest();
	void twofft();
	void vander();
	int zbrac();
	void zbrak();
	_REAL zbrent();
	void zroots();
#endif

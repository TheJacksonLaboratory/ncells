#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP dropout_with_centering (SEXP _x, SEXP _muvec, SEXP _p, SEXP _m1, SEXP _n, SEXP _pi1, SEXP _alpha, SEXP _logFC, SEXP _pp, SEXP _marginalsd, SEXP _marginalmean, SEXP _keep)
{
    double *x = REAL(_x);
    double *muvec = REAL(_muvec);
    double logFC = REAL(_logFC)[0];
    double alpha = REAL(_alpha)[0];
    double pp = REAL(_pp)[0];
    double marginalsd = REAL(_marginalsd)[0];
    double marginalmean = REAL(_marginalmean)[0];
    int p = INTEGER(_p)[0];
    int m1 = INTEGER(_m1)[0];
    int n = INTEGER(_n)[0];
    int pi1 = INTEGER(_pi1)[0];
    int *keep = INTEGER(_keep);

    int count = 0;
    double *s = Calloc(p, double);

    GetRNGstate();

    int i, j; // x = [x11, x12, ..., x1n, x21, x22, ..., x2n, ..., xp1, xp2, ..., xpn]
    for (i = 0; i < p; i++) {
	for (j = 0; j < n; j++) {
	    double xij = 0.0;
	    int ij = i*n + j;
	    if (i < m1 && j < pi1)
		xij = x[ij] + muvec[i] + logFC;
	    else
		xij = x[ij] + muvec[i];

	    double y0 = pnorm(xij, marginalmean, marginalsd, 1, 0);
	    double y1 = pnorm(xij, marginalmean + logFC, marginalsd, 1, 0);
	    if (alpha*log((1-pp)*y0 + pp*y1) < log(unif_rand())) {
		x[ij] = 0.0;
		count++;
	    } else {
		x[ij] = xij;
		keep[i] += 1;
	    }
	    s[i] += x[ij];
	}
    }

    for (i = 0; i < p; i++)
	for (j = 0; j < n; j++)
	    x[i*n + j] -= s[i]/n;

    /* Rprintf("zero proportion = %lf\n", (double) count/(p*n)); */

    PutRNGstate();

    Free(s);

    return R_NilValue;
}



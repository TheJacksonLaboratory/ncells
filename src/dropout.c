#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP dropout_with_centering (SEXP _x, SEXP _mu, SEXP _p, SEXP _p1, SEXP _n, SEXP _n1, SEXP _alpha, SEXP _logFC, SEXP _pp, SEXP _sd, SEXP _keep)
{
    double *x = REAL(_x);
    double *mu = REAL(_mu);
    double logFC = REAL(_logFC)[0];
    double alpha = REAL(_alpha)[0];
    double pp = REAL(_pp)[0];
    double sd = REAL(_sd)[0];
    int p = INTEGER(_p)[0];
    int p1 = INTEGER(_p1)[0];
    int n = INTEGER(_n)[0];
    int n1 = INTEGER(_n1)[0];
    int *keep = INTEGER(_keep);

    int count = 0;
    double *s = Calloc(p, double);

    GetRNGstate();

    int i, j; // x = [x11, x12, ..., x1n, x21, x22, ..., x2n, ..., xp1, xp2, ..., xpn]
    for (i = 0; i < p; i++) {
	for (j = 0; j < n; j++) {
	    double xij = 0.0;
	    int ij = i*n + j;
	    if (i < p1 && j < n1)
		xij = x[ij] + mu[i] + logFC;
	    else
		xij = x[ij] + mu[i];

	    double y0 = pnorm(xij, 0, sd, 1, 0);
	    double y1 = pnorm(xij, logFC, sd, 1, 0);
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


SEXP dropout_with_centering_cdf (SEXP _x, SEXP _mu, SEXP _p, SEXP _p1, SEXP _n, SEXP _n1, SEXP _alpha, SEXP _logFC, SEXP _cdf, SEXP _ncdf, SEXP _keep)
{
    double *x = REAL(_x);
    double *mu = REAL(_mu);
    double logFC = REAL(_logFC)[0];
    double alpha = REAL(_alpha)[0];
    double *cdf = REAL(_cdf);
    int ncdf = INTEGER(_ncdf)[0];
    int p = INTEGER(_p)[0];
    int p1 = INTEGER(_p1)[0];
    int n = INTEGER(_n)[0];
    int n1 = INTEGER(_n1)[0];
    int *keep = INTEGER(_keep);

    int pn = p*n;

    GetRNGstate();

    int i, j; // x = [x11, x12, ..., x1n, x21, x22, ..., x2n, ..., xp1, xp2, ..., xpn]
    for (i = 0; i < p; i++) {
	for (j = 0; j < n; j++) {
	    int ij = i*n + j;
	    if (i < p1 && j < n1)
		x[ij] = x[ij] + mu[i] + logFC;
	    else
		x[ij] = x[ij] + mu[i];
	}
    }

    int *ord = Calloc(pn, int);
    R_orderVector1(ord, pn, _x, 1, 0);

    /* Rprintf("ord[0] = %d\n", ord[0]); */

    int *val = Calloc(pn, int);
    j = 0;
    for (i = 0; i < pn; i++) {
	while (cdf[j] < x[ord[i]] && j < ncdf)
	    j += 1;
	val[ord[i]] = j - 1;
    }

    /* Rprintf("val[0] = %d\n", val[0]);     */

    int count = 0;
    double *s = Calloc(p, double);

    for (i = 0; i < p; i++) {
	for (j = 0; j < n; j++) {
	    int ij = i*n + j;
	    if (alpha*log((double)val[ij]/ncdf) < log(unif_rand())) {
		x[ij] = 0.0;
		count++;
	    } else {
		keep[i] += 1;
	    }
	    s[i] += x[ij];
	}
    }

    for (i = 0; i < p; i++)
	for (j = 0; j < n; j++)
	    x[i*n + j] -= s[i]/n;

    /* Rprintf("zero proportion = %lf\n", (double) count/pn); */

    PutRNGstate();

    Free(val);
    Free(ord);
    Free(s);

    return R_NilValue;
}



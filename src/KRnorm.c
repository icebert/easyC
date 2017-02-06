#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))



int KRnorm(float *A, int n, float *x)
{
    int   i, k;
    float Delta  = 3;
    float delta  = 0.1;
    float tol    = 1e-6;
    float g      = 0.9;
    float etamax = 0.1;
    float eta    = etamax;
    float stop_tol = tol * 0.5;
    float rt       = tol * tol;
    float *e       = (float *)malloc(n * sizeof(float));
    float *v       = (float *)malloc(n * sizeof(float));
    float *rk      = (float *)malloc(n * sizeof(float));
    float *y       = (float *)malloc(n * sizeof(float));
    float *Z       = (float *)malloc(n * sizeof(float));
    float *p       = (float *)malloc(n * sizeof(float));
    float *w       = (float *)malloc(n * sizeof(float));
    float *ap      = (float *)malloc(n * sizeof(float));
    float *ynew    = (float *)malloc(n * sizeof(float));
    float rho_km1, rho_km2, rout, rold, rat, eta_o, res_norm, innertol;
    float minynew, maxynew;
    float alpha, beta, gamma, tmp;

    for (i=0; i<n; i++)
    {
        e[i] = 1.0;
        x[i] = 1.0;
    }

    cblas_sgemv(CblasRowMajor, CblasNoTrans, n, n, 1, A, n, x, 1, 0, v, 1);
    for (i=0; i<n; i++)
    {
        v[i] *= x[i];
        rk[i] = 1.0 - v[i];
    }
    rho_km1 = cblas_sdot(n, rk, 1, rk, 1);
    rout    = rho_km1;
    rold    = rout;

    while (rout > rt)
    {
        k = 0;
        cblas_scopy(n, e, 1, y, 1);
        innertol = max(eta*eta*rout, rt);

        while (rho_km1 > innertol)
        {
            k++;
            if (k == 1)
            {
                for (i=0; i<n; i++)
                    if (v[i] == 0)
                        Z[i] = INFINITY;
                    else
                        Z[i] = rk[i] / v[i];
                cblas_scopy(n, Z, 1, p, 1);
                rho_km1 = cblas_sdot(n, rk, 1, Z, 1);
            }
            else
            {
                beta = rho_km1 / rho_km2;
                cblas_sscal(n, beta, p, 1);
                for (i=0; i<n; i++) p[i] += Z[i];
            }

            for (i=0; i<n; i++) w[i] = x[i] * p[i];
            cblas_sgemv(CblasRowMajor, CblasNoTrans, n, n, 1, A, n, w, 1, 0, ap, 1);
            for (i=0; i<n; i++) w[i] = (x[i] * ap[i]) + (v[i] * p[i]);
            tmp = cblas_sdot(n, p, 1, w, 1);
            if (tmp != 0)
                alpha = rho_km1 / tmp;
            else
                alpha = INFINITY;
            cblas_scopy(n, p, 1, ap, 1);
            cblas_sscal(n, alpha, ap, 1);

            for (i=0; i<n; i++) ynew[i] = y[i] + ap[i];
            minynew = ynew[0];
            maxynew = ynew[0];
            for (i=1; i<n; i++)
            {
                if (ynew[i] < minynew) minynew = ynew[i];
                if (ynew[i] > maxynew) maxynew = ynew[i];
            }
            if (minynew <= delta)
            {
                if (delta == 0) break;
                for (i=0; i<n; i++)
                    if (ap[i] < 0)
                    {
                        gamma = (delta - y[i]) / ap[i];
                        break;
                    }
                for (i++; i<n ;i++)
                    if (ap[i] < 0)
                    {
                        tmp = (delta - y[i]) / ap[i];
                        if (tmp < gamma) gamma = tmp;
                    }
                cblas_saxpy(n, gamma, ap, 1, y, 1);
                break;
            }
            if (maxynew >= Delta)
            {
                for (i=0; i<n; i++)
                    if (ynew[i] > Delta)
                    {
                        if (ap[i] == 0)
                            gamma = INFINITY;
                        else
                            gamma = (Delta - y[i]) / ap[i];
                        break;
                    }
                for (i++; i<n ;i++)
                    if (ynew[i] > Delta)
                    {
                        if (ap[i] == 0)
                            tmp = INFINITY;
                        else
                            tmp = (Delta - y[i]) / ap[i];
                        if (tmp < gamma) gamma = tmp;
                    }
                cblas_saxpy(n, gamma, ap, 1, y, 1);
                break;
            }
            cblas_scopy(n, ynew, 1, y, 1);
            cblas_saxpy(n, -alpha, w, 1, rk, 1);
            rho_km2 = rho_km1;
            for (i=0; i<n; i++)
                if (v[i] == 0)
                    Z[i] = INFINITY;
                else
                    Z[i] = rk[i] / v[i];
            rho_km1 = cblas_sdot(n, rk, 1, Z, 1);
        }

        for (i=0; i<n; i++) x[i] *= y[i];
        cblas_sgemv(CblasRowMajor, CblasNoTrans, n, n, 1, A, n, x, 1, 0, v, 1);
        for (i=0; i<n; i++)
        {
            v[i] *= x[i];
            rk[i] = 1.0 - v[i];
        }
        rho_km1  = cblas_sdot(n, rk, 1, rk, 1);
        rout     = rho_km1;
        rat      = rout / rold;
        rold     = rout;
        res_norm = sqrt(rout);
        eta_o    = eta;
        eta      = g * rat;

        if (g * eta_o * eta_o > 0.1) eta = max(eta, g * eta_o * eta_o);
        eta = max(min(eta, etamax), stop_tol / res_norm);
    }


    free(e);
    free(v);
    free(rk);
    free(y);
    free(Z);
    free(p);
    free(w);
    free(ap);
    free(ynew);

    return 0;
}




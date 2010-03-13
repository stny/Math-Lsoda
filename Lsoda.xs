#ifdef __cplusplus
extern "C" {
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"
#include <math.h>
  typedef int (*deriv_func)(double t, double* x, double* dx);
  typedef int (*jac_func)(double t, double* x, double** mat);
  typedef int (*lsoda_func_ptr) (const int* neq, const double* t, double* y, double* ydot);
  typedef int (*lsoda_jacobi_ptr) (const int* neq, const double* t, double* y, const int* ml, const int* mu, double* pd, const int* nrowpd);
  void dlsoda_(lsoda_func_ptr f, const int * neq, double* y,
               double* t, double const * tout, int const* itol, double const* rtol,
               double const* atol, int const* itask, int* istate,
               int const* iopt, double* rwork, int const* lrw, int* iwork,
               int const* liw, lsoda_jacobi_ptr jac, int const* jt);
#ifdef __cplusplus
}
#endif

typedef PerlIO * OutputStream;

#define XS_STATE(type, x)                           \
  INT2PTR(type, SvROK(x) ? SvIV(SvRV(x)) : SvIV(x))

#define XS_STRUCT2OBJ(sv, class, obj)           \
  if (obj == NULL) {                            \
    sv_setsv(sv, &PL_sv_undef);                 \
  } else {                                      \
    sv_setref_pv(sv, class, (void *) obj);      \
  }

static SV* func = (SV*)NULL;

void call_func(SV* func_name, double* t, SV* x, SV* y)
{
  dSP;
  ENTER;
  SAVETMPS;
  PUSHMARK(SP);
  XPUSHs(sv_2mortal(newSVnv(*t)));
  XPUSHs(sv_2mortal(newSVsv(x)));
  XPUSHs(sv_2mortal(newSVsv(y)));
  PUTBACK;
  call_sv(func_name, G_ARRAY);
  SPAGAIN;
  PUTBACK;
  FREETMPS;
  LEAVE;
}

int lsoda(int* dim, double* t, double* x, double* y)
{
  int size = *dim, i;
  SV* svx[size];
  SV* svy[size];
  for(i=0;i<size;i++){
    svx[i] = sv_2mortal(newSVnv(x[i]));
    svy[i] = sv_2mortal(newSVnv(y[i]));
  }
  SV* svx1 = newRV_noinc((SV*)av_make(size, svx));
  SV* svy1 = newRV_noinc((SV*)av_make(size, svy));

  call_func(func, t, svx1, svy1);
  AV* avy = (AV*)SvRV(svy1);
  for(i=0;i<size;i++){
    SV** pv = av_fetch(avy,i,0);
    y[i] = SvNV(*pv);
  }
}

int required_size(double start, double end, double dt)
{
  if(dt!=0){
    int size = (int)ceil( ((end - start)/dt) ) + 1;
    if(size > 0) return size;
  }
  return -1;
}

int
solve(SV* func_name, AV* y, double t, double tout, double dt, AV* rtol, AV* atol, OutputStream stream)
{
  func = func_name;
  int dim = av_len(y) + 1;
  int i, step;
  int maxvalue = 16 > (dim+9) ? 16 : (dim+9);
  int lrw = 22 + dim * maxvalue;
  int liw = 20 + dim;
  int itol = 2;
  int itask = 1;
  int istate = 1;
  int iopt = 0;
  int jt = 2;
  FILE *fp = PerlIO_findFILE(stream);
  int maxStep = required_size(t, tout, dt);
  double t1 = t + dt;
  double ay[dim], artol[dim], aatol[dim], rwork[lrw];
  int iwork[liw];
  for(i=0;i<dim;i++){
    SV** pv = av_fetch(y,i,0);
    ay[i] = SvNV(*pv);
    pv = av_fetch(rtol,i,0);
    artol[i] = SvNV(*pv);
    pv= av_fetch(atol,i,0);
    aatol[i] = SvNV(*pv);
  }

  fprintf(fp, "%g", t);
  for(i=0;i<dim;i++){
    fprintf(fp, "\t%g", ay[i]);
  }
  fprintf(fp, "\n");


  for(step=1;step<maxStep;step++){
    dlsoda_(lsoda, &dim, ay, &t, &t1, &itol, artol, aatol, &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, NULL, &jt);

    fprintf(fp, "%g", t);
    for(i=0;i<dim;i++){
    fprintf(fp, "\t%g", ay[i]);
    }
    fprintf(fp, "\n");

    t1 = t1 + dt;
    switch (istate) {
    case 2: /* successful exit */
      break;
    case -1:
      fprintf(stderr, "[WARNING]: excess work done at t =%14.6e (perhaps wrong JT).\n", t);
      istate = 1;
      return istate;
    case -2:
      fprintf(stderr, "[ERROR]: excess accuracy requested at t =%14.6e (tolerances too small).\n", t);
      return istate;
    case -3:
      fprintf(stderr, "[ERROR]: illegal input detected at t =%14.6e (see printed message).\n", t);
      return istate;
    case -4:
      fprintf(stderr, "[ERROR]: repeated error test failures at t =%14.6e (check all inputs).\n", t);
      return istate;
    case -5:
      fprintf(stderr, "[ERROR]: repeated convergence failures at t =%14.6e (perhaps bad Jacobian supplied or wrong choice of JT or tolerances).\n", t);
      return istate;
    case -6:
      fprintf(stderr, "[ERROR]: error weight became zero at t =%14.6e (Solution component i vanished, and ATOL or ATOL(i) = 0.)\n", t);
      return istate;
    case -7:
      fprintf(stderr, "[ERROR]: work space insufficient at t =%14.6e\n", t);
      return istate;
    default:
      fprintf(stderr, "[ERROR]: unknown error at t =%14.6e\n", t);
      return istate;
    }
  }
  return istate;
}


MODULE = Math::Lsoda  PACKAGE = Math::Lsoda

int
solve(func_name, y, t, tout, dt, rtol, atol, stream)
SV* func_name;
AV* y;
double t;
double tout;
double dt;
AV* rtol;
AV* atol;
OutputStream stream;

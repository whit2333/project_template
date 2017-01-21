#ifndef qcdnum2c_HH
#define qcdnum2c_HH 1
#include "TMath.h"
#ifdef __cplusplus
extern"C" {
#endif

   double func_(int * ipdf, double * x); // a test
   void   qcinit_(int *, char *);
   void   gxmake_(double *,int *,int *,int *,int *,int *);
   void   gqmake_(double *,double *,int *,int *,int *);
   void   fillwt_(int *,int *,int *,int *); 
   void   fillwc_(void (f)(double *, int *, int *, int *, int *, int *),int *,int *,int *); 
   void   setord_(int *);
   void   setalf_(double *, double * );
   void   evolfg_(int *,double (f)(int *,double *), double *,int * ,double *);
   void   setcbt_(int *,int *,int *, int *);
   int    iqfrmq_(double *);
   double fvalxq_(int *,int *, double *, double *, int *);
   double fpdfxq_(int *,double *, double *, double *, int *);
   double fsnsxq_(int *,int *, double *, double *, int *);
   double fsumxq_(int *,double *, double *, double *, int *);

   void   booktab_(double *, int *, int *, int *);
   void   myweight_(double *, int *, int *, int *, int *, int *);
   void   mysinglet_(double *, int *, int *, int *, int *, int *);

   void makewrs_(double *, int *,double (f)(double *, double *, int *), double (g)(double *, double *, int *), double (h)(double *),int *);
   void makewta_(double *, int *,double (f)(double *, double *, int *), double (g)(double *));
   void makewtd_(double *, int *,double (f)(double *, double *, int *), double (g)(double *));

   double achi_(double * );

   double pnsg2r_(double * , double * , int * );
   double pnsg2s_(double * , double * , int * );
   double pnsg2d_(double * , double * , int * );
   double pnsg2a_(double * , double * , int * );
   double pdum_(double * , double * , int * );

   double ptggr_(double * , double * , int * );
   double ptggs_(double * , double * , int * );
   double ptggd_(double * , double * , int * );
   double ptgga_(double * , double * , int * );
   double ptqga_(double * , double * , int * );


#ifdef __cplusplus
}
#endif

#endif


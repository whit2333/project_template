#include "qcdnum2c.h"
#include <iostream>

double func_(int * ipdf, double * x){
   // Function called from fortran

   //std::cout << "calling func_ with ipdf = " << ipdf[0] << " and x = " << x[0] << std::endl;
   if(ipdf[0] == 0 ) return x[0]*0.5; 
   if(ipdf[1] == 1 ) return x[0]*0.9; 
   return 0.0;
}


void myweight_(double * w, int * nw, int * nwords, int * idpij, int * mxord, int * idum) {

   printf("MADE IT HERE!!!!!\n");
//      dimension idpij(7,3)
   *nwords = 0;
//C--   S/R name for error message      
//      call setUmsg('MyWeight')
//C--   Max perturbative order
//      mxord = 1
   *mxord = 1;
//C--   Partition
//      itypes(1) = 2
//      itypes(2) = 2 
//      itypes(3) = 0
//      itypes(4) = 0
   int itypes[4];
   itypes[0] = 2;
   itypes[1] = 0;
   itypes[2] = 0;
   itypes[3] = 0;
//      call BookTab(w,nw,itypes,nwords)
   booktab_(w,nw,&itypes[0],nwords);
//C--   Not enough space       
//      if(nwords.le.0) return
   if( nwords <= 0 ) { printf("Error in myweight_\n");}
//C--   Assign table indices
//      idPij(1,1)   =  101                               ! PQQ
//      idPij(2,1)   =  201                               ! PQG
//      idPij(3,1)   =  102                               ! PGQ
//      idPij(4,1)   =  202                               ! PGG
//      idPij(5,1)   =  101                               ! PPL
//      idPij(6,1)   =  101                               ! PMI
//      idPij(7,1)   =  101                               ! PVA

   int ic = 0; 
   int jc = 0;
   ic = 0;               // id = 1 (Pqq)
   idpij[3*jc+ic] = 101;
   ic = 1;               // id = 2 (Pqg)
   idpij[3*jc+ic] = 101;
   ic = 2;               // id = 3 (Pgq)
   idpij[3*jc+ic] = 101;
   ic = 3;               // id = 4 (Pgg)
   idpij[3*jc+ic] = 102;
   ic = 4;               // id = 5 (P+)
   idpij[3*jc+ic] = 101;
   ic = 5;               // id = 6 (P-)
   idpij[3*jc+ic] = 101;
   ic = 6;               // id = 7 (Pv)
   idpij[3*jc+ic] = 101;

   int nodelta = 0;
//C--   Fill tables      
   ic = 0;              // id = 1 (Pqq)
   makewta_(w,&idpij[3*jc+ic],&pdum_, &achi_);

   ic = 3;              // id = 7 (Pv)
//      call MakeWRS(w, idPij(1,1), PQQR, PQQS, AChi, 0)
   makewrs_(w,&idpij[3*jc+ic],&pnsg2r_, &pnsg2s_, &achi_,&nodelta);
//      call MakeWtD(w, idPij(1,1), PQQD, AChi)
   makewtd_(w,&idpij[3*jc+ic],&pnsg2d_, &achi_);
//      call MakeWtA(w, idPij(2,1), PQGA, AChi)
   makewta_(w,&idpij[3*jc+ic],&pnsg2a_, &achi_);
//      call MakeWtA(w, idPij(3,1), PGQA, AChi)
//      call MakeWtA(w, idPij(4,1), PGGA, AChi)
//      call MakeWRS(w, idPij(4,1), PGGR, PGGS, AChi, 0)
//      call MakeWtD(w, idPij(4,1), PGGD, AChi)
//C--   Done!            
//      idum  = 0
   *idum = 0;
//      call clrUmsg
   printf("ALL DONE!!!!!!!!!!!!!!!!!!!!!!!1\n");
}

double achi_( double * mu2){
   // normal mellin transform so we just return zero.
   return 1.0;
}

double pnsg2r_(double * chi, double * mu2, int * nf){
   // chi is the coefficient of the scaling variable, chi = ax. It is simply one for normal mellin transform.
   // We calculate the regular part of PNS from Equation 1.2
   // Evolution equation for the structure function g(2) (x, Q**2)
   // Vladimir M. Braun (Regensburg U.), G.P. Korchemsky (Orsay), A.N. Manashov (Regensburg U. & Barcelona U., ECM). Feb 2001. 53 pp.
   // Published in Nucl.Phys. B603 (2001) 69-124
   //printf("pnsg2r_\n");
   double Nc = 3.0;
   double cf = (Nc*Nc-1.0)/(2.0*Nc);
   double cor = 2.0;
   return( 4.0*cf/cor);
}

double pnsg2s_(double * chi, double * mu2, int * nf){
   // chi is the coefficient of the scaling variable, chi = ax. It is simply one for normal mellin transform.
   // We calculate the singular part of PNS from Equation 1.2
   // Evolution equation for the structure function g(2) (x, Q**2)
   // Vladimir M. Braun (Regensburg U.), G.P. Korchemsky (Orsay), A.N. Manashov (Regensburg U. & Barcelona U., ECM). Feb 2001. 53 pp.
   // Published in Nucl.Phys. B603 (2001) 69-124
   //printf("pnsg2s_\n");
   double x = chi[0];
   double cor = 2.0;
   return( (1.0/(1.0-x))/cor );
}

double pnsg2a_(double * chi, double * mu2, int * nf){
   // chi is the coefficient of the scaling variable, chi = ax. It is simply one for normal mellin transform.
   // Evolution equation for the structure function g(2) (x, Q**2)
   // Vladimir M. Braun (Regensburg U.), G.P. Korchemsky (Orsay), A.N. Manashov (Regensburg U. & Barcelona U., ECM). Feb 2001. 53 pp.
   // Published in Nucl.Phys. B603 (2001) 69-124
   //printf("pnsg2a_\n");
   double Nc = 3.0;
   double cf = (Nc*Nc-1.0)/(2.0*Nc);
   double cor = 2.0;
   return( -2.0*cf/cor );
}

double pnsg2d_(double * chi, double * mu2, int * nf){
   // chi is the coefficient of the scaling variable, chi = ax. It is simply one for normal mellin transform.
   // Evolution equation for the structure function g(2) (x, Q**2)
   // Vladimir M. Braun (Regensburg U.), G.P. Korchemsky (Orsay), A.N. Manashov (Regensburg U. & Barcelona U., ECM). Feb 2001. 53 pp.
   // Published in Nucl.Phys. B603 (2001) 69-124
   //printf("pnsg2d_\n");
   double Nc = 3.0;
   double cf = (Nc*Nc-1.0)/(2.0*Nc);
   double cor = 2.0;
   return( (cf + (2.0 - TMath::Pi()*TMath::Pi()/3.0)/3.0)/cor );
}

double pdum_(double * chi, double * mu2, int * nf){
   return 0.0;
}

// --------------------------------------------------
// Singlet evolution 
// Note that ptqq is the same as the non-singlet
double ptggr_(double * chi, double * mu2, int * nf){
   // eqn 5.3
   return( 4.0*3.0);
}
double ptggs_(double * chi, double * mu2, int * nf){
   // eqn 5.3
   double x = chi[0];
   return( 1.0/(1.0-x));
}
double ptgga_(double * chi, double * mu2, int * nf){
   // eqn 5.3
   double x  = chi[0];
   double t1 = 3.0*(TMath::Pi()*TMath::Pi()/3.0 -2.0);
   double t2 = 3.0*TMath::Log((1.0-x)/x)*(2.0*TMath::Pi()*TMath::Pi()/3.0-6.0);
   return(t1 + t2);
}
double ptggd_(double * chi, double * mu2, int * nf){
   return(3.0*(TMath::Pi()*TMath::Pi()/3.0 - 1.0/3.0) - 2.0*nf[0]/3.0 );
}
double ptqga_(double * chi, double * mu2, int * nf){
   // eqn 5.3
   double x  = chi[0];
   return(-4.0*nf[0]*(x-2.0*TMath::Power(1.0-x,2.0)*TMath::Log(1.0-x)));
}

void mysinglet_(double * w, int * nw, int * nwords, int * idpij, int * mxord, int * idum) {
   printf("MADE IT HERE!!!!!\n");
   *nwords = 0;
   *mxord = 1;
   int itypes[4];
   itypes[0] = 2;// Pqq and zero
   itypes[1] = 2;// Pqg and Pgg
   itypes[2] = 0;
   itypes[3] = 0;
   booktab_(w,nw,&itypes[0],nwords);
   if( nwords <= 0 ) { printf("Error in myweight_\n");}

   int ic = 0; 
   int jc = 0;
   ic = 0;               // id = 1 (Pqq)
   idpij[3*jc+ic] = 102;
   ic = 1;               // id = 2 (Pqg)
   idpij[3*jc+ic] = 201;
   ic = 2;               // id = 3 (Pgq)
   idpij[3*jc+ic] = 101;
   ic = 3;               // id = 4 (Pgg)
   idpij[3*jc+ic] = 202;
   ic = 4;               // id = 5 (P+)
   idpij[3*jc+ic] = 102;
   ic = 5;               // id = 6 (P-)
   idpij[3*jc+ic] = 102;
   ic = 6;               // id = 7 (Pv)
   idpij[3*jc+ic] = 102;

   int nodelta = 0;
//C--   Fill tables      
   // 101
   ic = 0;              // id = 7 (Pqq)
   makewrs_(w,&idpij[3*jc+ic],&pnsg2r_, &pnsg2s_, &achi_,&nodelta);
   makewtd_(w,&idpij[3*jc+ic],&pnsg2d_, &achi_);
   makewta_(w,&idpij[3*jc+ic],&pnsg2a_, &achi_);

   // 102
   ic = 2;              // id = 1 (Pgq)
   makewta_(w,&idpij[3*jc+ic],&pdum_, &achi_);

   // 201
   ic = 1;              // id = 1 (Pqg)
   makewta_(w,&idpij[3*jc+ic],&ptqga_, &achi_);

   // 202
   ic = 3;              // id = 4 (Pgg)
   makewrs_(w,&idpij[3*jc+ic],&ptggr_, &ptggs_, &achi_,&nodelta);
   makewtd_(w,&idpij[3*jc+ic],&ptggd_, &achi_);
   makewta_(w,&idpij[3*jc+ic],&ptgga_, &achi_);

   // 201
//      call MakeWtA(w, idPij(3,1), PGQA, AChi)
//      call MakeWtA(w, idPij(4,1), PGGA, AChi)
//      call MakeWRS(w, idPij(4,1), PGGR, PGGS, AChi, 0)
//      call MakeWtD(w, idPij(4,1), PGGD, AChi)

//C--   Done!            
//      idum  = 0
   *idum = 0;
//      call clrUmsg
   printf("ALL DONE!!!!!!!!!!!!!!!!!!!!!!!1\n");
}


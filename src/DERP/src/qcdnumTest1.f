C     ----------------------------------------------------------------
      program ftwo
C     ----------------------------------------------------------------
C--   LO or NLO F2 with fast convolution engine, compare with ZMSTF
C--   With the convolution engine there are various ways to arrive at
C--   the same result (depending on taste and programming style). This
C--   is illustrated here by two F2 routines: getFtwo and altFtwo.
C     ----------------------------------------------------------------
C--   Michiel Botje   h24@nikhef.nl
      implicit double precision (a-h,o-z)
      data as0/0.364/, r20/2.D0/, iord/2/, nfin/0/ !alphas, NLO, VFNS
      external func                                !input parton dists
      dimension def(-6:6,12)                       !flavor composition
      data def  /
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     + 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,   !dval
     + 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,   !uval
     + 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,   !sval
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,   !dbar
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,   !ubar
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   !sbar
     + 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,   !cval
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   !cbar
     + 52*0.    /
      dimension xmin(5),iwt(5)
      data xmin/1.D-5,0.2D0,0.4D0,0.6D0,0.75D0/       !x grid
      data iwt / 1, 2, 4, 8, 16 /                     !x grid
      data nxin/100/, iosp/3/                         !x grid, splord
      dimension qq(2),wt(2)                           !mu2 grid
      data qq/2.D0,1.D4/, wt/1.D0,1.D0/, nqin/60/     !mu2 grid
      data q2c/3.D0/, q2b/25.D0/, q0/2.0/             !thresh and mu20
      dimension xpr(5)                                !printout      
      data xpr/1.D-4,1.D-3,1.D-2,1.D-1,.5D0/          !printout
      dimension qpr(1),ft(1),f2(1),fa(1)              !printout
      dimension pro(-6:6)                             !proton
      data pro / 4., 1., 4., 1., 4., 1., 0., 1., 4., 1., 4., 1., 4. /
      character*20 fname2,fname,zname      
      data fname /'weights/unpol.wgt'/
      data fname2 /'weights/pol.wgt'/
      data zname /'weights/zmstf.wgt'/

C--   lun = 6 stdout, -6 stdout w/out banner page
      lun    = 6
      lunout = abs(lun)      
      call qcinit(lun,' ')                     !initialize
      call gxmake(xmin,iwt,5,nxin,nx,iosp)     !x-grid
      call gqmake(qq,wt,2,nqin,nq)             !mu2-grid
      call readwt(22,fname,id1,id2,nw,ierr)    !read weights
      if(ierr.ne.0) then
        call fillwt(1,id1,id2,nw)              !make weights
        call fillwt(2,id1,id2,nw)              !make weights
        call dmpwgt(1,22,fname)                !dump weights
        call dmpwgt(2,22,fname2)                !dump weights
      endif
      call setord(iord)                        !LO, NLO, NNLO
      call setalf(as0,r20)                     !input alphas
      iqc  = iqfrmq(q2c)                       !charm threshold
      iqb  = iqfrmq(q2b)                       !bottom threshold
      call setcbt(nfin,iqc,iqb,999)            !thresholds in the vfns
      iq0  = iqfrmq(q0)                        !starting scale
      call evolfg(1,func,def,iq0,eps)          !evolve all pdf's
      call evolfg(2,func,def,iq0,eps)          !evolve all pdf's
      
C--   ZMSTF weight tables
      call zmreadw(22,zname,nw,ierr)           !read weights
      if(ierr.ne.0) then
        call zmfillw(nwords)                   !make weights
        call zmdumpw(22,zname)                 !dump weights
      endif  
C--   My own weight tables
      call myF2wgt     
C--   Proton
      do i = -6,6
        pro(i) = pro(i)/9.D0
      enddo
C--   Print F2: getF2 is my own F2 calculation, altF2 is alternative
      write(lunout,'(/5X,''x'',7X,''Q2'',7X,''F2(zmstf)'',
     +                6X,''F2(getftwo)'',4X,''F2(altftwo)''/)') 
      qpr(1) = 100.D0
      do ix = 1,5
        call zmStFun(2,pro,xpr(ix),qpr,ft,1,1)
        call getFtwo(  pro,xpr(ix),qpr,f2,1,1)
        call altFtwo(  pro,xpr(ix),qpr,fa,1,1)
        write(lunout,'(1PE9.1,1PE9.1,3E15.5)') xpr(ix),qpr,ft,f2,fa
      enddo                    
      end
      
C     ----------------------------------------------------------------

C     ====================================
      subroutine getFtwo(def,x,q,f,n,ichk)
C     ====================================

C--   This routine uses the fast convolution engine with the
C--   perturbative expansion capability of fastFxK. This gives
C--   compact code with QCDNUM caring about LO or NLO and
C--   multiplication by as(mu2)

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      common /weits/f2wgts(7000),idLO,idC2G,idC2Q
      
      dimension x(*), q(*), f(*)
      dimension def(-6:6)
      dimension idw(4)
      data idw /0,0,0,0/
      
      call getord(iord)
      if(iord.eq.3) stop 'GETFTWO: cannot handle NNLO'
C--   Setup fast buffers      
      call fastIni(x,q,n,ichk)
C--   F2 quarks      
      call fastSns(1,def,7,1)        ! 1 = all quarks
      idw(1) = idLO                  ! LO
      idw(2) = idC2Q                 ! NLO 
      call fastFxK(f2wgts,idw,1,3)   ! 3 = F2(quarks) 
C--   F2 qluon
      call fastSns(1,def,0,1)        ! 1 = gluon * <e2>
      idw(1) = 0                     ! LO no gluon contribution
      idw(2) = idC2G                 ! NLO 
      call fastFxK(f2wgts,idw,1,2)   ! 2 = F2(gluon)
      call fastCpy(2,3,1)            ! 3 = F2(quarks)+F2(gluon)
C--   Interpolate      
      call fastFxq(3,f,n)  
        
      return
      end
      
C     ====================================
      subroutine altFtwo(def,x,q,f,n,ichk)
C     ====================================

C--   This routine uses the fast convolution engine but not the
C--   perturbative expansion capability of fastFxK. Here we need
C--   explicit multiplication with alfas, depending on LO or NLO

C--   Michiel Botje   h24@nikhef.nl

      implicit double precision (a-h,o-z)
      
      common /weits/f2wgts(7000),idLO,idC2G,idC2Q
      
      dimension x(*), q(*), f(*)
      dimension def(-6:6)
      dimension idw(4)
      data idw /0,0,0,0/
      external timesAs
      
      call getord(iord)
      if(iord.eq.3) stop 'ALTFTWO: cannot handle NNLO'
C--   Setup fast buffers      
      call fastIni(x,q,n,ichk)
C--   F2 LO      
      call fastSns(1,def,7,1)          ! 1 = all quarks
      call fastCpy(1,4,0)              ! 4 = F2(LO)
      if(iord.eq.2) then
C--     NLO quarks
        idw(1) = idC2Q                 ! C2Q 
        call fastFxK(f2wgts,idw,1,3)   ! 3 = quarks * C2Q 
C--     NLO qluon
        call fastSns(1,def,0,1)        ! 1 = gluon * <e2>
        idw(1) = idC2G                 ! C2G
        call fastFxK(f2wgts,idw,1,2)   ! 2 = gluon * C2G
C--     Add NLO quark and gluon and multiply by alfas        
        call fastCpy(2,3,1)            ! 3 = q*C2Q + g*C2G
        call fastKin(3,timesAs)        ! 3 = F2(NLO)
C--     Add LO and NLO        
        call fastcpy(3,4,1)            ! 4 = F2(LO) + F2(NLO)
      endif
C--   Interpolate      
      call fastFxq(4,f,n)  
        
      return
      end
      
C     ===================================================      
      double precision function timesAs(ix,iq,nf,ithresh)
C     ===================================================

C--   Returns as/2pi  (with correct renormalisation scale dependence!)

      implicit double precision (a-h,o-z)
      
      idum    = ix                     !avoid compiler warning
      idum    = nf                     !avoid compiler warning
      idum    = ithresh                !avoid compiler warning
      timesAs = getalfn(iq,-1,ierr)
      
      return
      end

C     ----------------------------------------------------------------      

C     ==================
      subroutine myF2wgt      
C     ==================

      implicit double precision (a-h,o-z)
      
      common /weits/f2wgts(7000),idLO,idC2G,idC2Q
      
      external df2Achi, df2Const1, df2C2QB, df2C2GA
      
      dimension itypes(4)
      data itypes /0,0,0,0/
      
C--   Book two type-1 tables and one type-2 table 
      itypes(1) = 2
      itypes(2) = 1     
      call booktab(f2wgts,7000,itypes,nwords)
C--   LO weights 
      idLO  = 101
      call makeWtD(f2wgts,idLO,df2Const1,df2Achi)     
C--   NLO weights
      idC2Q = 102
      call makeWtB(f2wgts,idC2Q,df2C2QB,df2Achi,0)
      idC2G = 201
      call makeWtA(f2wgts,idC2G,df2C2GA,df2Achi)
      
      return
      end 

C     =======================================
      double precision function df2Achi(qmu2)
C     =======================================

C--   Returns coefficient 'a' of scaling function chi = a*x
C--   For zero-mass coefficient functions chi = x and a = 1

      implicit double precision (a-h,o-z)

      dummy   = qmu2     !avoid compiler warning
      df2Achi = 1.D0

      return
      end
      
C     ==============================================
      double precision function df2Const1(x,qmu2,nf)
C     ==============================================

C--   Coefficient of LO delta(1-x) terms
C--   Returns the value 1.D0

      implicit double precision (a-h,o-z)

C--   Avoid compiler warnings
      dummy     = x
      dummy     = qmu2
      idumm     = nf
      df2Const1 = 1.D0

      return
      end
      
C     ============================================
      DOUBLE PRECISION FUNCTION df2C2QB(X,QMU2,NF)
C     ============================================

C--   Returns the value of C2Q(x)     singular

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      IDUM    = NF          !avoid compiler warning
      DDUM    = QMU2        !avoid compiler warning
      C1MX    = 1.-X
      df2C2QB = 3. + (5.D0/3.D0)*X + ((4.D0/3.D0)*LOG(C1MX/X)-1.) * 
     &          (1.+X**2) / C1MX
 
      RETURN
      END

C     ============================================
      DOUBLE PRECISION FUNCTION df2C2GA(X,QMU2,NF)
C     ============================================

C--   Returns the value of C2G(x,nf) regular

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      DDUM    = QMU2        !avoid compiler warning
      C1MX    = 1. - X
      C2G     = -.5 + 4.*X*C1MX + .5 * (X**2+C1MX**2) * LOG(C1MX/X)
      df2C2GA = 2.*NF*C2G
 
      RETURN
      END
      
C     ----------------------------------------------------------------

C     ======================================
      double precision function func(ipdf,x)
C     ======================================

      implicit double precision (a-h,o-z)

                     func = 0.D0
      if(ipdf.eq. 0) func = xglu(x)
      if(ipdf.eq. 1) func = xdnv(x)
      if(ipdf.eq. 2) func = xupv(x)
      if(ipdf.eq. 3) func = 0.D0
      if(ipdf.eq. 4) func = xdbar(x)
      if(ipdf.eq. 5) func = xubar(x)
      if(ipdf.eq. 6) func = xsbar(x)
      if(ipdf.eq. 7) func = 0.D0
      if(ipdf.eq. 8) func = 0.D0
      if(ipdf.eq. 9) func = 0.D0
      if(ipdf.eq.10) func = 0.D0
      if(ipdf.eq.11) func = 0.D0
      if(ipdf.eq.12) func = 0.D0

      return
      end
 
C     =================================
      double precision function xupv(x)
C     =================================

      implicit double precision (a-h,o-z)

      data au /5.107200D0/

      xupv = au * x**0.8D0 * (1.D0-x)**3.D0

      return
      end

C     =================================
      double precision function xdnv(x)
C     =================================

      implicit double precision (a-h,o-z)

      data ad /3.064320D0/

      xdnv = ad * x**0.8D0 * (1.D0-x)**4.D0

      return
      end

C     =================================
      double precision function xglu(x)
C     =================================

      implicit double precision (a-h,o-z)

      data ag    /1.7D0 /

      common /msum/ glsum, uvsum, dvsum

      xglu = ag * x**(-0.1D0) * (1.D0-x)**5.D0

      return
      end

C     ==================================
      double precision function xdbar(x)
C     ==================================

      implicit double precision (a-h,o-z)

      data adbar /0.1939875D0/

      xdbar = adbar * x**(-0.1D0) * (1.D0-x)**6.D0

      return
      end

C     ==================================
      double precision function xubar(x)
C     ==================================

      implicit double precision (a-h,o-z)

      xubar = xdbar(x) * (1.D0-x)

      return
      end

C     ==================================
      double precision function xsbar(x)
C     ==================================

      implicit double precision (a-h,o-z)

      xsbar = 0.2D0 * (xdbar(x)+xubar(x))

      return
      end








      SUBROUTINE PLAN3 (IGEOM)
C-----------------------------------------------------------------------
C
C     Compute pressure and velocity using consistent approximation spaces.     
C     Operator splitting technique.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'
      include 'GEOM' ! ROBIN_BC
      include 'MASS' ! ROBIN_BC
C
      COMMON /SCRNS/  RESV1 (LX1,LY1,LZ1,LELV)
     $ ,              RESV2 (LX1,LY1,LZ1,LELV)
     $ ,              RESV3 (LX1,LY1,LZ1,LELV)
     $ ,              DV1   (LX1,LY1,LZ1,LELV)
     $ ,              DV2   (LX1,LY1,LZ1,LELV)
     $ ,              DV3   (LX1,LY1,LZ1,LELV)
      COMMON /SCRVH/  H1    (LX1,LY1,LZ1,LELV)
     $ ,              H2    (LX1,LY1,LZ1,LELV)
      COMMON /MYROB/  H2_ROB(LX1,LY1,LZ1,LELV ) ! ROBIN_BC
      n = nx1*ny1*nz1*nelv ! ROBIN_BC
C
      IF (IGEOM.EQ.1) THEN
C
C        Old geometry
C
         CALL MAKEF
C
      ELSE
C
C        New geometry, new b.c.
C
         intype = -1
         call sethlm  (h1,h2,intype)
         ! this is moved to ophinv!!!
C------------------------------------------------------------------!
C--------Create H2 that takes into account for Robin BC
         call rzero(H2_ROB,N) ! ROBIN_BC
         call bcneusc_mod(H2_ROB) ! ROBIN_BC
C--------Add contribution of Robin BC to standard Neumann condition!
         call add2 (H2,H2_ROB,N) ! ROBIN_BC
C------------------------------------------------------------------!
         call cresvif (resv1,resv2,resv3,h1,h2)

         call ophinv  (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
         call opadd2  (vx,vy,vz,dv1,dv2,dv3)
c
         call incomprn(vx,vy,vz,pr)
C
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE LAGPRES 
C--------------------------------------------------------------------
C
C     Keep old pressure values
C
C--------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'

      common /cgeom/ igeom

      IF (NBDINP.EQ.3.and.igeom.le.2) THEN
         NTOT2 = lx2*ly2*lz2*NELV
         CALL COPY (PRLAG,PR,NTOT2)
      ENDIF
      RETURN
      END
C
      subroutine cresvif (resv1,resv2,resv3,h1,h2)
C---------------------------------------------------------------------
C
C     Compute startresidual/right-hand-side in the velocity solver
C
C---------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      REAL           RESV1 (LX1,LY1,LZ1,1)
      REAL           RESV2 (LX1,LY1,LZ1,1)
      REAL           RESV3 (LX1,LY1,LZ1,1)
      REAL           H1    (LX1,LY1,LZ1,1)
      REAL           H2    (LX1,LY1,LZ1,1)
      COMMON /SCRUZ/ W1    (LX1,LY1,LZ1,LELV)
     $ ,             W2    (LX1,LY1,LZ1,LELV)
     $ ,             W3    (LX1,LY1,LZ1,LELV)

      common /cgeom/ igeom
      integer in_slipper

      NTOT1 = lx1*ly1*lz1*NELV
      NTOT2 = lx2*ly2*lz2*NELV
      if (igeom.eq.2) CALL LAGVEL 
      CALL BCDIRVC (VX,VY,VZ,v1mask,v2mask,v3mask)
      CALL BCNEUTR
C
      call extrapp (pr,prlag)
      call opgradt (resv1,resv2,resv3,pr)
      CALL OPADD2  (RESV1,RESV2,RESV3,BFX,BFY,BFZ)
      CALL OPHX    (W1,W2,W3,VX,VY,VZ,H1,H2)
      CALL OPSUB2  (RESV1,RESV2,RESV3,W1,W2,W3)
C
      RETURN
      END
c-----------------------------------------------------------------------
C ********************************************************************
      SUBROUTINE BCNEUSC_MOD(S,slip_dir) ! PICELLA
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'NEKUSE'
      include 'SHS'
      include 'ADJOINT'
      DIMENSION S(LX1,LY1,LZ1,LELV)
      CHARACTER CB*3
      integer i2dx,i2dy,mem
      real NU
      integer slip_dir
      real v_addon
      NU = abs(param(02))
      !if (nio.eq.0) write (6,*) 'SLIP_LENGTH' , SLIP_LENGTH , 'NU' , NU
      NFACES =2*NDIM
      NXYZ = NX1*NY1*NZ1
      IFIELD =1 ! flow variables
      NEL = NELFLD(IFIELD)
      NTOT = NXYZ*NEL
      CALL RZERO(S,NTOT)
c     if(nio.eq.0)
c     $write (6 ,*) ' SUBROUTINE BCNEUSC_MOD '
c    $ , NFACES , NEL , NTOT , LELV
      DO 1000 IE =1 , NEL
      DO 1000 IFACE =1 , NFACES
      ieg = lglel (ie)
      in_slipper = object_inverse(ieg,iobj_store)
      CB = CBC (IFACE,IE,IFIELD)
      IF (CB.EQ .'SYM'.and.in_slipper.ne.0) THEN ! Apply the condition to ’SYM ’ faces only ...
      IA =0
C
C IA is areal counter, assumes advancing fastest index first. (IX ...IY...IZ)
C
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFACE )
      DO 100 IZ = KZ1 , KZ2
      DO 100 IY = KY1 , KY2
      DO 100 IX = KX1 , KX2
      IA = IA + 1
!     S (IX,IY,IZ,IE)=1/ slip_length
      call map_3Dtoslipper(I2dx,I2dy,mem,IX,IY,IZ,ie,IFACE,1) ! on obj1
      !slip_length = slipper(I2dx,I2dy,1,mem,1) +1e-10 
      slip_length = slipper(I2dx,I2dy,1,mem,1) ! note! this is inverse
      ! slip length now!!! so we multiply not times! 
      S(IX,IY,IZ,IE) = S(IX,IY,IZ,IE) +
     $ NU*Slip_length*AREA(IA,1,IFACE,IE)/BM1(IX,IY,IZ,IE)
      
C     if(nio.eq.0) write(6,*) IE, S(IX,IY,IZ,IE)
 100  CONTINUE
      endif
 1000 CONTINUE
      RETURN
      END
C

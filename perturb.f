c-----------------------------------------------------------------------
      subroutine fluidp (igeom)
c
c     Driver for perturbation velocity
c
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'

      do jp=1,npert

         if (nio.eq.0.and.igeom.eq.2) write(6,1) istep,time,jp
   1     format(i9,1pe14.7,' Perturbation Solve:',i5)

         call perturbv (igeom)

      enddo

      jp=0   ! set jp to zero, for baseline flow

      return
      end
c-----------------------------------------------------------------------
      subroutine perturbv (igeom)
c
c
c     Solve the convection-diffusion equation for the perturbation field, 
c     with projection onto a div-free space.
c
c
      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'
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
      real L2, tmpH1, SEMI, LINF
      n = nx1*ny1*nz1*nelv ! ROBIN_BC

c
      ifield = 1
c
      if (igeom.eq.1) then
c
c        Old geometry, old velocity
c
         call makefp
         call lagfieldp
c
      else
c
c        New geometry, new velocity
c
         intype = -1
         call sethlm   (h1,h2,intype)
C------------------------------------------------------------------!
C--------Create H2 that takes into account for Robin BC
         call rzero(H2_ROB,N) ! ROBIN_BC
         call bcneusc_mod(H2_ROB) ! ROBIN_BC
C--------Add contribution of Robin BC to standard Neumann condition!
         call add2 (H2,H2_ROB,N) ! ROBIN_BC
C------------------------------------------------------------------!
         call cresvipp (resv1,resv2,resv3,h1,h2)
         call ophinv   (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
         call opadd2   (vxp(1,jp),vyp(1,jp),vzp(1,jp),dv1,dv2,dv3)
         call incomprp (vxp(1,jp),vyp(1,jp),vzp(1,jp),prp(1,jp))
c
        !! HARRY
        ! zero mass flow rate vol_flow!
        call vol_flowp
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine lagfieldp
c
c     Keep old Vp-field(s)
c
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
c
      do ilag=nbdinp-1,2,-1
         call opcopy
     $     (vxlagp(1,ilag  ,jp),vylagp(1,ilag  ,jp),vzlagp(1,ilag  ,jp)
     $     ,vxlagp(1,ilag-1,jp),vylagp(1,ilag-1,jp),vzlagp(1,ilag-1,jp))
      enddo
      call opcopy(vxlagp(1,1,jp),vylagp(1,1,jp),vzlagp(1,1,jp)
     $           ,vxp   (1,jp)  ,vyp   (1,jp)  ,vzp   (1,jp) )
c
      return
      end
c-----------------------------------------------------------------------
      subroutine makefp
c
c     Make rhs for velocity perturbation equation
c
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'
      include 'ADJOINT'
 
                                              call makeufp
      if (ifnav.and.(.not.ifchar).and.(.not.ifadj))call advabp
      if (ifnav.and.(.not.ifchar).and.(     ifadj))call advabp_adjoint
      if (iftran)                              call makextp
                                               call makebdfp
c
      return
      end
c--------------------------------------------------------------------
      subroutine makeufp
c
c     Compute and add: (1) user specified forcing function (FX,FY,FZ)
c
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
C
      time = time-dt
      call nekuf   (bfxp(1,jp),bfyp(1,jp),bfzp(1,jp))
      call opcolv2 (bfxp(1,jp),bfyp(1,jp),bfzp(1,jp)
     $                              ,vtrans(1,1,1,1,ifield),bm1)
      time = time+dt
c
      return
      end
c--------------------------------------------------------------------
      subroutine advabp
C
C     Eulerian scheme, add convection term to forcing function
C     at current time step.
C
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
C
      COMMON /SCRNS/ TA1 (LX1*LY1*LZ1*LELV)
     $ ,             TA2 (LX1*LY1*LZ1*LELV)
     $ ,             TA3 (LX1*LY1*LZ1*LELV)
     $ ,             TB1 (LX1*LY1*LZ1*LELV)
     $ ,             TB2 (LX1*LY1*LZ1*LELV)
     $ ,             TB3 (LX1*LY1*LZ1*LELV)
C
      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv
c
      if (if3d) then
         call opcopy  (tb1,tb2,tb3,vx,vy,vz)                   ! Save velocity
         call opcopy  (vx,vy,vz,vxp(1,jp),vyp(1,jp),vzp(1,jp)) ! U <-- dU
         call convop  (ta1,tb1)                                ! du.grad U
         call convop  (ta2,tb2)
         call convop  (ta3,tb3)
         call opcopy  (vx,vy,vz,tb1,tb2,tb3)  ! Restore velocity
c
         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
            bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
         enddo
c
         call convop  (ta1,vxp(1,jp))       !  U.grad dU
         call convop  (ta2,vyp(1,jp))
         call convop  (ta3,vzp(1,jp))
c
         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
            bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
         enddo
c
      else  ! 2D
c
         call opcopy  (tb1,tb2,tb3,vx,vy,vz)                   ! Save velocity
         call opcopy  (vx,vy,vz,vxp(1,jp),vyp(1,jp),vzp(1,jp)) ! U <-- dU
         call convop  (ta1,tb1)                                ! du.grad U
         call convop  (ta2,tb2)
         call opcopy  (vx,vy,vz,tb1,tb2,tb3)  ! Restore velocity
c
         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
         enddo
c
         call convop  (ta1,vxp(1,jp))       !  U.grad dU
         call convop  (ta2,vyp(1,jp))
c
         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
         enddo
c
      endif
c
      return
      end
c--------------------------------------------------------------------
      subroutine advabp_adjoint
C
C     Eulerian scheme, add convection term to forcing function
C     at current time step for backward part of adjoint: 
C     Convective term is now (U.Grad)u - (Grad U)^T .u
C     instead of (u.Grad)U + (U.Grad)u in above subroutine advabp
C
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'GEOM'
      include 'ADJOINT'
C
      COMMON /SCRNS/ TA1 (LX1*LY1*LZ1*LELV)
     $ ,             TA2 (LX1*LY1*LZ1*LELV)
     $ ,             TA3 (LX1*LY1*LZ1*LELV)
     $ ,             TB1 (LX1*LY1*LZ1*LELV)
     $ ,             TB2 (LX1*LY1*LZ1*LELV)
     $ ,             TB3 (LX1*LY1*LZ1*LELV)


      real fact,factx,facty
C
      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv   !dimensionn arrays 
      NTOT      = lx1*ly1*lz1*NELT

c
      if (if3d) then
         call opcopy  (tb1,tb2,tb3,vx,vy,vz)                   ! Save velocity
         call opcopy  (vx,vy,vz,vxp(1,jp),vyp(1,jp),vzp(1,jp)) ! U <-- u
c     
         call convop_adj  (ta1,ta2,ta3,tb1,tb2,tb3,vx,vy,vz) ! u.grad U^T

         call opcopy  (vx,vy,vz,tb1,tb2,tb3)  ! Restore velocity
c
         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
            bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
         enddo
c
c               
c
         call convop  (ta1,vxp(1,jp))       !  U.grad u
         call convop  (ta2,vyp(1,jp))
         call convop  (ta3,vzp(1,jp))
c
         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)+tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)+tmp*ta2(i)
            bfzp(i,jp) = bfzp(i,jp)+tmp*ta3(i)
         enddo
c
         if (ifheat) then       ! dt.(grad T)
c                               ! term coming from the temperature convection
            call opcolv3 (ta1,ta2,ta3,dTdx,dTdy,dTdz,tp)
c
            do i=1,ntot1
               tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,2)
               bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
               bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
               bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
            enddo
         endif
c
      else  ! 2D

         call opcopy  (tb1,tb2,tb3,vx,vy,vz)                   ! Save velocity
         call opcopy  (vx,vy,vz,vxp(1,jp),vyp(1,jp),vzp(1,jp)) 

         call convop_adj  (ta1,ta2,ta3,tb1,tb2,tb3,vx,vy,vz) ! u.((grad U)^T)

         call opcopy  (vx,vy,vz,tb1,tb2,tb3) ! Restore velocity

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
         enddo

         call convop  (ta1,vxp(1,jp))       !  U.grad u
         call convop  (ta2,vyp(1,jp))
c
         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)+tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)+tmp*ta2(i)
         enddo

         if (ifheat) then       ! dt.(grad T)^T
                                ! term coming from the temperature convection
            call opcolv3 (ta1,ta2,ta3,dTdx,dTdy,dTdz,tp)
c
            do i=1,ntot1
               tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,2)
               bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
               bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
            enddo
         endif
         
      endif
c
      return
      end
c--------------------------------------------------------------------
      subroutine convop_adj  (bdux,bduy,bduz,udx,udy,udz,cx,cy,cz)

      include 'SIZE'
      include 'TOTAL'

      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $     , uf1(ltd),uf2(ltd),uf3(ltd),uf4(ltd),uf5(ltd),uf6(ltd)
      real urx(ltd),usx(ltd),utx(ltd)
      real ury(ltd),usy(ltd),uty(ltd)
      real urz(ltd),usz(ltd),utz(ltd)
      real bdux(1),bduy(1),bduz(1),u(1),cx(1),cy(1),cz(1)
      real udx(1),udy(1),udz(1)
      logical ifuf,ifcf            ! u and/or c already on fine mesh?
      integer e
      real bdrx(1), bdry(1),bdrz (1)

      call set_dealias_rx
      nxyz1 = lx1*ly1*lz1
c     AM DING DING 
      nxyzd = lxd*lyd*lzd
      nxyzu = nxyz1
      nxyzc = nxyz1
      ntot1=lx1*ly1*lz1*nelv
      ic = 1                    ! pointer to vector field C
      iu = 1                    ! pointer to scalar field u
      ib = 1                    ! pointer to scalar field Bdu
      if(if3d)then
         do e=1,nelv
                                ! map coarse velocity to fine mesh (C-->F)
            call intp_rstd(fx,cx(ic),lx1,lxd,if3d,0) ! 0 --> forward
            call intp_rstd(fy,cy(ic),lx1,lxd,if3d,0) 
            call intp_rstd(fz,cz(ic),lx1,lxd,if3d,0) 
               
            call intp_rstd(uf1,udx(iu),lx1,lxd,if3d,0) ! 0 --> forward
            call grad_rst(urx,usx,utx,uf1,lxd,if3d)
            
            call intp_rstd(uf2,udy(iu),lx1,lxd,if3d,0) 
            call grad_rst(ury,usy,uty,uf2,lxd,if3d)
            
            call intp_rstd(uf3,udz(iu),lx1,lxd,if3d,0) 
            call grad_rst(urz,usz,utz,uf3,lxd,if3d)
            
            do i=1,nxyzd        ! mass matrix included, per DFM (4.8.5)
               uf4(i)=fx(i)*(rx(i,1,e)*urx(i)+rx(i,4,e)*usx(i)
     $              +rx(i,7,e)*utx(i))+
     $              fy(i)*(rx(i,1,e)*ury(i)+rx(i,4,e)*usy(i)
     $              +rx(i,7,e)*uty(i))+
     $              fz(i)*(rx(i,1,e)*urz(i)+rx(i,4,e)*usz(i)
     $              +rx(i,7,e)*utz(i))
               uf5(i)=fx(i)*(rx(i,2,e)*urx(i)+rx(i,5,e)*usx(i)
     $              +rx(i,8,e)*utx(i))+
     $              fy(i)*(rx(i,2,e)*ury(i)+rx(i,5,e)*usy(i)
     $              +rx(i,8,e)*uty(i))+
     $              fz(i)*(rx(i,2,e)*urz(i)+rx(i,5,e)*usz(i)
     $              +rx(i,8,e)*utz(i))
               uf6(i)=fx(i)*(rx(i,3,e)*urx(i)+rx(i,6,e)*usx(i)
     $              +rx(i,9,e)*utx(i))+
     $              fy(i)*(rx(i,3,e)*ury(i)+rx(i,6,e)*usy(i)
     $              +rx(i,9,e)*uty(i))+
     $              fz(i)*(rx(i,3,e)*urz(i)+rx(i,6,e)*usz(i)
     $              +rx(i,9,e)*utz(i))
            enddo

            call intp_rstd(bdux(ib),uf4,lx1,lxd,if3d,1) ! Project back to coarse
            call intp_rstd(bduy(ib),uf5,lx1,lxd,if3d,1)
            call intp_rstd(bduz(ib),uf6,lx1,lxd,if3d,1)

            ic = ic + nxyzc
            iu = iu + nxyzu
            ib = ib + nxyz1
         enddo
         call invcol2     (bdux,bm1,ntot1) ! local mass inverse
         call invcol2     (bduy,bm1,ntot1) ! local mass inverse
         call invcol2     (bduz,bm1,ntot1) ! local mass inverse
      else
         do e=1,nelv

                               ! map coarse velocity to fine mesh (C-->F)
            call intp_rstd(fx,cx(ic),lx1,lxd,if3d,0) ! 0 --> forward
            call intp_rstd(fy,cy(ic),lx1,lxd,if3d,0) 

            call intp_rstd(uf1,udx(iu),lx1,lxd,if3d,0) 
            call grad_rst(urx,usx,utx,uf1,lxd,if3d)

            call intp_rstd(uf2,udy(iu),lx1,lxd,if3d,0) 
            call grad_rst(ury,usy,uty,uf2,lxd,if3d)

            do i=1,nxyzd       
               uf4(i) = fx(i)*(rx(i,1,e)*urx(i)+rx(i,3,e)*usx(i))+
     $              fy(i)*(rx(i,1,e)*ury(i)+rx(i,3,e)*usy(i))
               uf5(i) = fx(i)*(rx(i,2,e)*urx(i)+rx(i,4,e)*usx(i))+
     $              fy(i)*(rx(i,2,e)*ury(i)+rx(i,4,e)*usy(i))
            enddo

            call intp_rstd(bdux(ib),uf4,lx1,lxd,if3d,1)
            call intp_rstd(bduy(ib),uf5,lx1,lxd,if3d,1)

            ic = ic + nxyzc
            iu = iu + nxyzu
            ib = ib + nxyz1
         end do
         call invcol2     (bdux,bm1,ntot1) ! local mass inverse
         call invcol2     (bduy,bm1,ntot1) ! local mass inverse
      end if


      return
      end
c--------------------------------------------------------------------
      subroutine makextp
c
c     Add extrapolation terms to perturbation source terms
c
c     (nek5 equivalent for velocity is "makeabf")
c
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
C
      common /scrns/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)
c
      ntot1 = lx1*ly1*lz1*nelv
c
      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)
      call add3s2 (ta1,exx1p(1,jp),exx2p(1,jp),ab1,ab2,ntot1)
      call add3s2 (ta2,exy1p(1,jp),exy2p(1,jp),ab1,ab2,ntot1)
      call copy   (exx2p(1,jp),exx1p(1,jp),ntot1)
      call copy   (exy2p(1,jp),exy1p(1,jp),ntot1)
      call copy   (exx1p(1,jp),bfxp (1,jp),ntot1)
      call copy   (exy1p(1,jp),bfyp (1,jp),ntot1)
      call add2s1 (bfxp(1,jp),ta1,ab0,ntot1)
      call add2s1 (bfyp(1,jp),ta2,ab0,ntot1)
      if (if3d) then
         call add3s2 (ta3,exz1p(1,jp),exz2p(1,jp),ab1,ab2,ntot1)
         call copy   (exz2p(1,jp),exz1p(1,jp),ntot1)
         call copy   (exz1p(1,jp),bfzp (1,jp),ntot1)
         call add2s1 (bfzp(1,jp),ta3,ab0,ntot1)
      endif
c
      return
      end
c--------------------------------------------------------------------
      subroutine makebdfp
C
C     Add contributions to perturbation source from lagged BD terms.
C
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
C
      COMMON /SCRNS/ TA1(LX1,LY1,LZ1,LELV)
     $ ,             TA2(LX1,LY1,LZ1,LELV)
     $ ,             TA3(LX1,LY1,LZ1,LELV)
     $ ,             TB1(LX1,LY1,LZ1,LELV)
     $ ,             TB2(LX1,LY1,LZ1,LELV)
     $ ,             TB3(LX1,LY1,LZ1,LELV)
     $ ,             H2 (LX1,LY1,LZ1,LELV)
C
      ntot1 = lx1*ly1*lz1*nelv
      const = 1./dt
      call cmult2(h2,vtrans(1,1,1,1,ifield),const,ntot1)
      call opcolv3c (tb1,tb2,tb3
     $              ,vxp(1,jp),vyp(1,jp),vzp(1,jp),bm1,bd(2))
C
      do ilag=2,nbd
         if (ifgeom) then
            call opcolv3c(ta1,ta2,ta3,vxlagp(1,ilag-1,jp),
     $                                vylagp(1,ilag-1,jp),
     $                                vzlagp(1,ilag-1,jp),
     $                                bm1lag(1,1,1,1,ilag-1),bd(ilag+1))
         else
            call opcolv3c(ta1,ta2,ta3,vxlagp(1,ilag-1,jp),
     $                                vylagp(1,ilag-1,jp),
     $                                vzlagp(1,ilag-1,jp),
     $                                bm1                   ,bd(ilag+1))
         endif
         call opadd2  (tb1,tb2,tb3,ta1,ta2,ta3)
      enddo
      call opadd2col (bfxp(1,jp),bfyp(1,jp),bfzp(1,jp),tb1,tb2,tb3,h2)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine cresvipp(resv1,resv2,resv3,h1,h2)
c
c     Account for inhomogeneous Dirichlet boundary contributions 
c     in rhs of perturbation eqn.
c                                               n
c     Also, subtract off best estimate of grad p
c
      include 'SIZE'
      include 'TOTAL'
      real           resv1 (lx1,ly1,lz1,1)
      real           resv2 (lx1,ly1,lz1,1)
      real           resv3 (lx1,ly1,lz1,1)
      real           h1    (lx1,ly1,lz1,1)
      real           h2    (lx1,ly1,lz1,1)
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             prextr(lx2,ly2,lz2,lelv)
      real tmph1,semi,l2,linf
c
      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv
c
c-----------------------------------------------------------------------
      ! HARRY 
      ! This was stolen from Valerio/Guilluam
c     add a diagonal term to the operator for the ON/on/O/o boundary conditions if the adjoint equations are solved
      call adjonbc(H2)
c-----------------------------------------------------------------------
      call bcdirvc (vxp(1,jp),vyp(1,jp),vzp(1,jp)
     $             ,v1mask,v2mask,v3mask)
      call rzero(prextr,lx1*ly1*lz1*lelv)
      call extrapprp(prextr)
      call opgradt(resv1,resv2,resv3,prextr)
      call opadd2(resv1,resv2,resv3,bfxp(1,jp),bfyp(1,jp),bfzp(1,jp))
      call ophx  (w1,w2,w3,vxp(1,jp),vyp(1,jp),vzp(1,jp),h1,h2)
      call opsub2(resv1,resv2,resv3,w1,w2,w3)
c
      return
      end
c--------------------------------------------------------------------
      subroutine heatp (igeom)
C
C     Driver for temperature or passive scalar.
C
C     Current version:
C     (1) Varaiable properties.
C     (2) Implicit time stepping.
C     (3) User specified tolerance for the Helmholtz solver
C         (not based on eigenvalues).
C     (4) A passive scalar can be defined on either the 
C         temperatur or the velocity mesh.
C     (5) A passive scalar has its own multiplicity (B.C.).  
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      INCLUDE 'SOLN'
C
      do jp=1,npert
         do ifield=2,nfield
            INTYPE        = -1
            IF (.NOT.IFTMSH(IFIELD)) IMESH = 1
            IF (     IFTMSH(IFIELD)) IMESH = 2
            CALL UNORM
            CALL SETTOLT
            CALL CDSCALP (IGEOM)
         enddo
      enddo
c
      jp=0   ! set jp to zero, for baseline flow
c
      return
      end
C
C-----------------------------------------------------------------------
      subroutine cdscalp (igeom)
C-----------------------------------------------------------------------
C
C     Solve the convection-diffusion equation for passive scalar IPSCAL
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'MVGEOM'
      INCLUDE 'SOLN'
      INCLUDE 'MASS'
      INCLUDE 'TSTEP'
      COMMON  /CPRINT/ IFPRINT
      LOGICAL          IFPRINT
      LOGICAL          IFCONV
C
      COMMON /SCRNS/ TA(LX1,LY1,LZ1,LELT)
     $              ,TB(LX1,LY1,LZ1,LELT)
      COMMON /SCRVH/ H1(LX1,LY1,LZ1,LELT)
     $              ,H2(LX1,LY1,LZ1,LELT)
c
      include 'ORTHOT'
      ifld1 = ifield-1
      napproxt(1,ifld1) = laxtt
C
      IF (IGEOM.EQ.1) THEN
C
C        Old geometry
C
         CALL MAKEQP
         CALL LAGSCALP
C
      ELSE
C
         IF (IFPRINT) THEN
            IF (IFIELD.EQ.2.AND.NIO.EQ.0) THEN
               WRITE (6,*) ' Temperature/Passive scalar solution'
            ENDIF
         ENDIF
         if1=ifield-1
         write(name4t,1) if1
    1    format('TEM',i1)
C
C        New geometry
C
         NEL    = NELFLD(IFIELD)
         NTOT   = lx1*ly1*lz1*NEL
C
         INTYPE = 0
         IF (IFTRAN) INTYPE = -1
         CALL SETHLM  (H1,H2,INTYPE)
         CALL BCNEUSC (TA,-1)
         CALL ADD2    (H2,TA,NTOT)
         CALL BCDIRSC (TP(1,IFIELD-1,jp))
         CALL AXHELM  (TA,TP (1,IFIELD-1,jp),H1,H2,IMESH,1)
         CALL SUB3    (TB,BQP(1,IFIELD-1,jp),TA,NTOT)
         CALL BCNEUSC (TA,1)
         CALL ADD2    (TB,TA,NTOT)
c
         CALL HMHOLTZ (name4t,TA,TB,H1,H2
     $                 ,TMASK(1,1,1,1,IFIELD-1)
     $                 ,TMULT(1,1,1,1,IFIELD-1)
     $                 ,IMESH,TOLHT(IFIELD),NMXT(ifield-1),1)
c
         CALL ADD2    (TP(1,IFIELD-1,jp),TA,NTOT)
C
         CALL BCNEUSC (TA,1)
         CALL ADD2 (BQP(1,IFIELD-1,jp),TA,NTOT)
C
      ENDIF
C
      return
      end
C
      subroutine makeqp
C----------------------------------------------------------------------
C
C     Generate forcing function for the solution of a passive scalar.
C     !! NOTE: Do not change the content of the array BQ until the current
C              time step is completed 
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      include 'SOLN'
                                                       CALL MAKEUQP
      IF (IFADVC(IFIELD).AND.(.NOT.IFCHAR))            CALL CONVABP
      IF (IFTRAN)                                      CALL MAKEABQP
      IF ((IFTRAN.AND..NOT.IFCHAR).OR.
     $    (IFTRAN.AND..NOT.IFADVC(IFIELD).AND.IFCHAR)) CALL MAKEBDQP
c     IF (IFADVC(IFIELD).AND.IFCHAR)                   CALL CONVCHP
      IF (IFADVC(IFIELD).AND.IFCHAR) write(6,*) 'no convchp'
      IF (IFADVC(IFIELD).AND.IFCHAR) call exitt
c
      return
      end
C
      subroutine makeuqp
C---------------------------------------------------------------------
C
C     Fill up user defined forcing function and collocate will the
C     mass matrix on the Gauss-Lobatto mesh.
C
C---------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'MASS'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
C
      NTOT = lx1*ly1*lz1*NELFLD(IFIELD)
C
      time = time-dt                           ! time is tn
c
      call rzero   ( bqp(1,ifield-1,jp) ,    ntot)
      call setqvol ( bqp(1,ifield-1,jp)          )
      call col2    ( bqp(1,ifield-1,jp) ,bm1,ntot)
c
      time = time+dt                           ! restore time
C
      return
      end
C
      subroutine convabp
C---------------------------------------------------------------
C
C     Eulerian scheme, add convection term to forcing function 
C     at current time step.
C
C---------------------------------------------------------------
      include 'SIZE'
      include 'ADJOINT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
c
      common /scruz/ ta (lx1,ly1,lz1,lelt)
     $             , ua (lx1,ly1,lz1,lelt)
     $             , ub (lx1,ly1,lz1,lelt)
     $             , uc (lx1,ly1,lz1,lelt)
      real coeff
c
      nel = nelfld(ifield)
      ntot1 = lx1*ly1*lz1*nel
c
      if (.not.ifadj) then
         call opcopy(ua,ub,uc,vx,vy,vz)
         call opcopy(vx,vy,vz,vxp(1,jp),vyp(1,jp),vzp(1,jp))
         call convop(ta,t(1,1,1,1,ifield-1)) ! dU.grad T
         call opcopy(vx,vy,vz,ua,ub,uc)
      endif
c
      call convop  (ua,tp(1,ifield-1,jp)) ! U.grad dT
c
      if (.not.ifadj) then
         call add2 (ta,ua,ntot1)          ! U.grad dT + dU.grad T
      else
         call copy (ta,ua,ntot1)
         coeff=-1.0
         call cmult(ta,coeff,ntot1)       ! -U.grad dT
                                          ! the second term depends on the buoyancy
      endif
c
      call col2    (ta,vtrans(1,1,1,1,ifield),ntot1) !vtrans (U.grad dT + dU.grad T)
      call subcol3 (bqp(1,ifield-1,jp),bm1,ta,ntot1)
c
      return
      end
C
      subroutine makeabqp
C-----------------------------------------------------------------------
C
C     Sum up contributions to 3rd order Adams-Bashforth scheme.
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
C
      COMMON /SCRUZ/ TA (LX1,LY1,LZ1,LELT)
C
      AB0   = AB(1)
      AB1   = AB(2)
      AB2   = AB(3)
      NEL   = NELFLD(IFIELD)
      NTOT1 = lx1*ly1*lz1*NEL
C
      CALL ADD3S2 (TA,VGRADT1P(1,IFIELD-1,jp),
     $                VGRADT2P(1,IFIELD-1,jp),AB1,AB2,NTOT1)
      CALL COPY   (   VGRADT2P(1,IFIELD-1,jp),
     $                VGRADT1P(1,IFIELD-1,jp),NTOT1)
      CALL COPY   (   VGRADT1P(1,IFIELD-1,jp),
     $                     BQP(1,IFIELD-1,jp),NTOT1)
      CALL ADD2S1 (BQP(1,IFIELD-1,jp),TA,AB0,NTOT1)
C
      return
      end
C
      subroutine makebdqp
C-----------------------------------------------------------------------
C
C     Add contributions to Q from lagged BD terms.
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'MASS'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
C
      COMMON /SCRNS/ TA (LX1,LY1,LZ1,LELT)
     $ ,             TB (LX1,LY1,LZ1,LELT)
     $ ,             H2 (LX1,LY1,LZ1,LELT)
C
      NEL   = NELFLD(IFIELD)
      NTOT1 = lx1*ly1*lz1*NEL
      CONST = 1./DT
      CALL COPY  (H2,VTRANS(1,1,1,1,IFIELD),NTOT1)
      CALL CMULT (H2,CONST,NTOT1)
C
      CALL COL3  (TB,BM1,TP(1,IFIELD-1,jp),NTOT1)
      CALL CMULT (TB,BD(2),NTOT1)
C
      DO 100 ILAG=2,NBD
         IF (IFGEOM) THEN
            CALL COL3 (TA,BM1LAG(1,1,1,1,ILAG-1),
     $                    TLAGP (1,ILAG-1,IFIELD-1,jp),NTOT1)
         ELSE
            CALL COL3 (TA,BM1,
     $                    TLAGP (1,ILAG-1,IFIELD-1,jp),NTOT1)
         ENDIF
         CALL ADD2S2(TB,TA,BD(ilag+1),NTOT1)
 100  CONTINUE
C
      CALL COL2 (TB,H2,NTOT1)
      CALL ADD2 (BQP(1,IFIELD-1,jp),TB,NTOT1)
C
      return
      end
C
      subroutine lagscalp
C-----------------------------------------------------------------------
C
C     Keep old passive scalar field(s) 
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
C
      NTOT1 = lx1*ly1*lz1*NELFLD(IFIELD)
C
      DO 100 ILAG=NBDINP-1,2,-1
         CALL COPY (TLAGP(1,ILAG  ,IFIELD-1,jp),
     $              TLAGP(1,ILAG-1,IFIELD-1,jp),NTOT1)
 100  CONTINUE
C
      CALL COPY (TLAGP(1,1,IFIELD-1,jp),TP(1,IFIELD-1,jp),NTOT1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine incomprp (ux,uy,uz,up)
c
c     Project U onto the closest incompressible field
c
c     Input:  U     := (ux,uy,uz)
c
c     Output: updated values of U, iproj, proj; and
c             up    := pressure correction req'd to impose div U = 0
c
c
c     Dependencies: ifield ==> which "density" (vtrans) is used.
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
c
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)
      COMMON /SCRCH/ PREXTR(LX2,LY2,LZ2,LELV)
      logical ifprjp

c
      if (icalld.eq.0) tpres=0.0
      icalld=icalld+1
      npres=icalld
c
      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      intype = 1
      dtbd   = bd(1)/dt

      call rzero   (h1,ntot1)
c     call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call cmult2  (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)
      call invers2 (h2inv,h2,ntot1)

      call opdiv   (dp,ux,uy,uz)
      call chsign  (dp,ntot2)
      call ortho   (dp)


C******************************************************************


      ifprjp=.false.    ! project out previous pressure solutions?
      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0) ifprjp=.true.

      ! Most likely, the following can be commented out. (pff, 1/6/2010)
c     if (npert.gt.1.or.ifbase)            ifprjp=.false.
cpff  if (ifprjp)   call setrhs  (dp,h1,h2,h2inv)

                    call esolver (dp,h1,h2,h2inv,intype)

cpff  if (ifprjp)   call gensoln (dp,h1,h2,h2inv)

cNOTE:  The "cpff" comments added 11/24/17 to avoid old-style projection,
cNOTE:  which should be replaced with something more updated.

C******************************************************************

      call opgradt (w1 ,w2 ,w3 ,dp)
      call opbinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
      call opadd2  (ux ,uy ,uz ,dv1,dv2,dv3)

      call extrapprp(prextr)
      call lagpresp
      call add3(up,prextr,dp,ntot2)

      return
      end
c------------------------------------------------------------------------
      subroutine extrapprp (prextr)
C
C     Pressure extrapolation
C
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      COMMON /CTMP0/ DPR (LX2,LY2,LZ2,LELV)
      REAL        PREXTR (LX2,LY2,LZ2,LELV)
C
      ntot2 = lx2*ly2*lz2*nelv
      if (nbd.eq.1.or.nbd.eq.2) then
         call copy (prextr,prp(1,JP),ntot2)
      elseif (nbd.eq.3) then
         const = dtlag(1)/dtlag(2)
         call sub3 (dpr,prp(1,JP),prlagp(1,1,JP),ntot2)
         call cmult(dpr,const,ntot2)
         call add3 (prextr,prp(1,JP),dpr,ntot2)
      elseif (nbd.gt.3) then
         write (6,*) 'Pressure extrapolation cannot be completed'
         write (6,*) 'Try a lower-order temporal scheme'
         call exitt
      endif
      return
      end
C-------------------------------------------------------------------
      subroutine lagpresp
C
C     save old pressure values
C
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
C
      if (nbdinp.eq.3) then
         ntot2 = lx2*ly2*lz2*nelv
         call copy (prlagp(1,1,JP),prp(1,JP),ntot2)
      endif
      return
      end
C-------------------------------------------------------------------
      subroutine lyap_scale ! Rescale / orthogonalize perturbation fields
c
c
      include 'SIZE'
      include 'TOTAL'
c
      real sigma(0:lpert)
c
      ntotv = lx1*ly1*lz1*nelv
      ntotp = lx2*ly2*lz2*nelv
      ntott = lx1*ly1*lz1*nelt
c
      do j=1,npert
         call normvc(h1,semi,pl2,plinf,vxp(1,j),vyp(1,j),vzp(1,j))
         sigma(j) = pl2
         if (pl2.gt.0) then
            write(6,*) 'this is pl2:',pl2
            scale = 1./pl2
            call opcmult(vxp(1,j),vyp(1,j),vzp(1,j),scale)
            call   cmult(tp(1,1,j),scale,ntott)
            call   cmult(prp(1,j) ,scale,ntotp)
         endif
c
c        Have to do lag terms as well, etc
c        
c        Also, must orthogonalize
      enddo
c
      if (nio.eq.0) then
         if (npert.eq.1) write(6,1) istep,time,(sigma(j),j=1,npert)
         if (npert.eq.2) write(6,2) istep,time,(sigma(j),j=1,npert)
         if (npert.eq.3) write(6,3) istep,time,(sigma(j),j=1,npert)
         if (npert.eq.4) write(6,4) istep,time,(sigma(j),j=1,npert)
         if (npert.eq.5) write(6,5) istep,time,(sigma(j),j=1,npert)
         if (npert.eq.6) write(6,6) istep,time,(sigma(j),j=1,npert)
         if (npert.eq.7) write(6,7) istep,time,(sigma(j),j=1,npert)
         if (npert.eq.8) write(6,8) istep,time,(sigma(j),j=1,npert)
         if (npert.eq.9) write(6,9) istep,time,(sigma(j),j=1,npert)
      endif
    1 format(i8,1p2e12.4,' pgrow')
    2 format(i8,1p3e12.4,' pgrow')
    3 format(i8,1p4e12.4,' pgrow')
    4 format(i8,1p5e12.4,' pgrow')
    5 format(i8,1p6e12.4,' pgrow')
    6 format(i8,1p7e12.4,' pgrow')
    7 format(i8,1p8e12.4,' pgrow')
    8 format(i8,1p9e12.4,' pgrow')
    9 format(i8,1p10e12.4,' pgrow')
c
      return
      end
c-----------------------------------------------------------------------
      subroutine out_pert  ! dump perturbation .fld files
c
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
c
      character*3 s3
c
      do jpp=1,npert
         write(s3,3) jpp
 3       format('p',i2.2)
         call outpost2
     $   (vxp(1,jpp),vyp(1,jpp),vzp(1,jpp),prp(1,jpp),tp(1,1,jpp),1,s3)
c        call writehist
c    $     (vxp(1,jpp),vyp(1,jpp),vzp(1,jpp),tp(1,1,jpp),jpp)
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine pert_add2s2(i,j,scale)   ! xi = xi + scale * xj
      include 'SIZE'
      include 'TOTAL'

      ntotp = lx2*ly2*lz2*nelv
      ntotv = lx1*ly1*lz1*nelv
      ntott = lx1*ly1*lz1*nelt

      call add2s2(vxp(1,i),vxp(1,j),scale,ntotv)
      call add2s2(vyp(1,i),vyp(1,j),scale,ntotv)
      if (if3d) call add2s2(vzp(1,i),vzp(1,j),scale,ntotv)
      call add2s2(prp(1,i),prp(1,j),scale,ntotp)

      do l=1,lorder-1
         call add2s2(vxlagp(1,l,i),vxlagp(1,l,j),scale,ntotv)
         call add2s2(vylagp(1,l,i),vylagp(1,l,j),scale,ntotv)
         if (if3d) call add2s2(vzlagp(1,l,i),vzlagp(1,l,j),scale,ntotv)
         if (l.le.lorder-2) 
     $      call add2s2(prlagp(1,l,i),prlagp(1,l,j),scale,ntotp)
      enddo

      call add2s2(exx1p(1,i),exx1p(1,j),scale,ntotv)
      call add2s2(exy1p(1,i),exy1p(1,j),scale,ntotv)
      if (if3d) call add2s2(exz1p(1,i),exz1p(1,j),scale,ntotv)

      call add2s2(exx2p(1,i),exx2p(1,j),scale,ntotv)
      call add2s2(exy2p(1,i),exy2p(1,j),scale,ntotv)
      if (if3d) call add2s2(exz2p(1,i),exz2p(1,j),scale,ntotv)

      if (ifheat) then
        do k=0,npscal
          k1=k+1
          ntotk = lx1*ly1*lz1*nelfld(k+2)
          call add2s2(tp(1,k1,i),tp(1,k1,j),scale,ntotk)
          do l=1,lorder-1
            call add2s2(tlagp(1,l,k1,i),tlagp(1,l,k1,j),scale,ntotk)
          enddo
          call add2s2(vgradt1p(1,k1,i),vgradt1p(1,k1,j),scale,ntotk)
          call add2s2(vgradt2p(1,k1,i),vgradt2p(1,k1,j),scale,ntotk)
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      function pert_inner_prod(i,j) ! weighted inner product vi^T vj
      include 'SIZE'
      include 'TOTAL'

      common/normset/pran, ray, rayc

      ntotv=lx1*ly1*lz1*nelv
      ntott=lx1*ly1*lz1*nelt

      s1 = rayc*glsc3(vxp(1,i),bm1,vxp(1,j),ntotv)
      s2 = rayc*glsc3(vyp(1,i),bm1,vyp(1,j),ntotv)
      s3 = 0
      if (if3d) s3 = rayc*glsc3(vzp(1,i),bm1,vzp(1,j),ntotv)

      t1 = 0
      if (ifheat) t1=pran*ray*ray*glsc3(tp(1,1,i),bm1,tp(1,1,j),ntott)

      pert_inner_prod = (s1+s2+s3+t1)/volvm1

      return
      end
c-----------------------------------------------------------------------
      subroutine pert_ortho_norm ! orthogonalize and rescale pert. arrays
      include 'SIZE'
      include 'TOTAL'

      do k=1,npert
         call pert_ortho_norm1(k)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pert_ortho_norm1 (k) ! orthogonalize k against 1...k-1
      include 'SIZE'
      include 'TOTAL'

      do j=1,k-1
         scale = -pert_inner_prod(j,k)
         call pert_add2s2(k,j,scale)   ! xi = xi + scale * xj
      enddo
      scale = pert_inner_prod(k,k)
      if (scale.gt.0) scale = 1./scale
      call rescalepert(pertnorm,scale,k)

      return
      end
c-----------------------------------------------------------------------
      function opnorm2(v1,v2,v3)
      include 'SIZE'
      include 'TOTAL'
c
      real v1(1) , v2(1), v3(1)
      real normsq1,normsq2,normsq3,opnorm
c
      ntotv=lx1*ly1*lz1*nelv
      normsq1=glsc3(v1,bm1,v1,ntotv)
      normsq2=glsc3(v2,bm1,v2,ntotv)
      if(if3d) then
         normsq3=glsc3(v3,bm1,v3,ntotv)
      else
         normsq3=0
      endif

      opnorm2=normsq1+normsq2+normsq3
      if (opnorm2.gt.0) opnorm2=sqrt(opnorm2/volvm1)
      return
      end
c-----------------------------------------------------------------------

      function Tnorm(temp)
      include 'SIZE'
      include 'TOTAL'

      real temp(*)
c
      ntotv = lx1*ly1*lz1*nelv
      Tnorm = sqrt( glsc3(temp,BM1,temp,ntotv) /voltm1)
c
      return
      end
c--------------------------------------------
      function dmnorm1(v1,v2,v3,temp)
c     Norm weighted by mass matrix
      include 'SIZE'
      include 'TOTAL'

      real v1(1),v2(1),v3(1),temp(1)
      real normsq1,normsq2,normsq3,normsqT,dMnorm
      common/normset/pran, ray, rayc

      ntotv=lx1*ly1*lz1*nelv
      normsq1=(rayc)*glsc3(v1,BM1,v1,ntotv)
      normsq2=(rayc)*glsc3(v2,BM1,v2,ntotv)
      if(if3d) then
         normsq3=(rayc)*glsc3(v3,BM1,v3,ntotv)
      else
         normsq3=0
      endif

      if(ifheat) then
         normsqT = (pran*ray*ray)*glsc3(temp,BM1,temp,ntotv)
      else
         normsqT = 0
      endif

      dmnorm1=sqrt((normsq1+normsq2+normsq3+normsqT)/volvm1)

      return
      end

c---------------------------------------------------------------
      subroutine opscale(v1,v2,v3,temp,alpha)
c     v =  alpha*v
      include 'SIZE'
      include 'INPUT'

      real alpha
      real v1(1),v2(1),v3(1),temp(1)

      ltotv=lx1*ly1*lz1*lelv
      ltott=lx1*ly1*lz1*lelt

      call cmult(v1,alpha,ltotv)
      call cmult(v2,alpha,ltotv)
      if (if3d)   call cmult(v3,alpha,ltotv)
      if (ifheat) call cmult(temp,alpha,ltott*ldimt)

      return
      end

c---------------------------------------------------------------
      subroutine opscaleV(v1,v2,v3,alpha)
c     v =  alpha*v
      include 'SIZE'
      include 'INPUT'
      real alpha
      real v1(*),v2(*),v3(*)

      ntotv=lx1*ly1*lz1*nelv

      call cmult(v1,alpha,ntotv)
      call cmult(v2,alpha,ntotv)

      if (if3d)   call cmult(v3,alpha,ntotv)
c
      return
      end

c-----------------------------------------------------------------------
      subroutine computelyap
      include 'SIZE'
      include 'TOTAL'

      do jpp=1,npert         ! Loop through each Lyapunov eigenvalue
         call computelyap1
     $     (vxp(1,jpp),vyp(1,jpp),vzp(1,jpp),tp(1,1,jpp),jpp)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine computelyap1(vxq,vyq,vzq,tq,jpp)
      include 'SIZE'
      include 'TOTAL'

      real vxq(1),vyq(1),vzq(1),tq(1)
      real lyapsum,twt
      common /pertsave/ timestart,tinitial

      integer icount
      save    icount
      data    icount /0/

      logical       if_restart,if_ortho_lyap
      common/restar/if_restart,if_ortho_lyap

      character*132 lyprestart
      common/restflename/lyprestart  !file for restart data

      twt = param(126) !time to wait to start computing exponents

      if (nio.eq.0) 
     $  write(6,8) istep,icount,time,twt,(lyap(k,jpp),k=1,3),jpp
    8 format(i9,i4,2f8.4,1p3e12.4,i3,'clyap')

      if(time.lt.twt) then
c
c        For a  fresh simulation, then we wait 5 vertical diffusion 
c        times before we start measuring, so in this case just rescale
c
         pertnorm = dmnorm1(vxq,vyq,vzq,tq)
         pertinvnorm = 1.0/pertnorm
         call rescalepert(pertnorm,pertinvnorm,jpp)
         lyap(3,jpp) = pertnorm !store latest norm
         timestart   = time
         tinitial    = time
         icount      = 0
         return
      else
         if (jpp.eq.1) icount = icount + 1
      endif

      irescale = param(128)
      if (icount.eq.irescale) then

         lyapsum     = lyap(2,jpp)
         oldpertnorm = lyap(3,jpp)
         pertnorm=dmnorm1(vxq,vyq,vzq,tq)
c
         lyap(1,jpp) = log(pertnorm/oldpertnorm)/(time-timestart)
         lyapsum     = lyapSum + log(pertnorm/oldpertnorm)
         lyap(2,jpp) = lyapSum

         if(nid.eq.0) then        ! write out results to the .lyp file
 
           write(6 ,1) istep,time,lyap(1,jpp),lyapsum,pertnorm,jpp
           write(79,2) time,lyap(1,jpp),lyapsum,pertnorm,oldpertnorm,jpp
 1          format(i9,1p4e17.8,i4,'lyap')
 2          format(1p5e17.8,i4,'lyap')
c           call flushbuffer(79)
 
            if (jpp.eq.1) open(unit=96,file=lyprestart)
            write(96,9) lyapsum,timestart,jpp
 9          format(1p2e19.11,i9)
            if (jpp.eq.npert) close(96)
         endif

         pertinvnorm = 1.0/pertnorm
         call rescalepert(pertnorm,pertinvnorm,jpp)
         lyap(3,jpp) = pertnorm  !save current pertnorm as old pertnorm

         if (jpp.eq.npert) then
            icount    = 0
            timestart = time
         endif

         if_ortho_lyap = .false.
         if (param(125).gt.0) if_ortho_lyap = .true.
         if (jpp.eq.npert .and. if_ortho_lyap) call pert_ortho_norm

      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine rescalepert(pertnorm,pertinvnorm,jpp)
      include 'SIZE'
      include 'TOTAL'

      ntotp = lx2*ly2*lz2*nelv

      call opscale                     !normalize vectors to unit norm
     $      (vxp(1,jpp),vyp(1,jpp),vzp(1,jpp),tp(1,1,jpp),pertinvnorm)
      call cmult(prp(1,jpp),pertinvnorm,ntotp)

      call opscale(exx1p(1,jpp),exy1p(1,jpp),exz1p(1,jpp)
     $                           ,vgradt1p(1,1,jpp),pertinvnorm)
      call opscale(exx2p(1,jpp),exy2p(1,jpp),exz2p(1,jpp)
     $                           ,vgradt2p(1,1,jpp),pertinvnorm)

      ltotv = lx1*ly1*lz1*lelv
      ltotp = lx2*ly2*lz2*lelv

      call cmult( tlagp(1,1,1,jpp),pertinvnorm,ltotv*(lorder-1)*ldimt)
      call cmult(vxlagp(1,1,jpp),pertinvnorm,ltotv*(lorder-1))
      call cmult(vylagp(1,1,jpp),pertinvnorm,ltotv*(lorder-1))
      call cmult(vzlagp(1,1,jpp),pertinvnorm,ltotv*(lorder-1))
      call cmult(prlagp(1,1,jpp),pertinvnorm,ltotp*(Lorder-2))

      if (nio.eq.0) write(6,1) istep,pertnorm,pertinvnorm,jpp,'PNORM'
  1   format(i4,1p2e12.4,i4,a5)
      pertnorm = pertnorm*pertinvnorm

      return
      end
c-----------------------------------------------------------------------
      ! HARRY 
      ! Also stolen from Valerio/Guillaume
c-----------------------------------------------------------------------

      subroutine adjonbc(H2)
      implicit none

      include 'SIZE'
      include 'MASS'
      include 'PARALLEL'
      include 'SOLN'
      include 'ADJOINT'
      include 'TSTEP'
      include 'INPUT'
      include 'GEOM'

      real H2(lx1,ly1,lz1,1)
      character cb*3
      integer nfaces,nxyz,nel,ntot,iel,iface,ieg
      integer kx1,kx2,ky1,ky2,kz1,kz2,ix,iy,iz,ia
      logical ifalgn,ifnorx,ifnory,ifnorz

      nfaces=2*ndim
      nxyz  =nx1*ny1*nz1
      nel   =nelfld(ifield)
      ntot  =nxyz*nel
c
c     Add diagonal terms to the matrix for adjoint O/o and ON/on
c     boundary conditions. Does nothing for direct equations.
c
      if (ifadj) then
c
c     check which faces have O/o/ON/on conditions and modify the H2 matrix
c     matrix accordingly. The equations must include a Robin term related
c     to the base flow velocity normal to the wall.
c
         do iel=1,nel
            do iface=1,nfaces
               ieg=lglel(iel)
               cb = cbc(iface,iel,ifield)
               if (cb.eq.'O  '.or.cb.eq.'o  '
     &              .or.cb.eq.'ON '.or.cb.eq.'on ') then
                  ia = 0
c
c     Loop on the points of the face
                  call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
                  do iz=kz1,kz2
                     do iy=ky1,ky2
                        do ix=kx1,kx2
                           ia = ia + 1
c     add an U_n coefficient to H2
                           H2(ix,iy,iz,iel)=H2(ix,iy,iz,iel)+
     &                          area(ia,1,iface,iel)/bm1(ix,iy,iz,iel)*
     &                          (unx(ia,1,iface,iel)*vx(ix,iy,iz,iel)
     &                          +uny(ia,1,iface,iel)*vy(ix,iy,iz,iel)
     &                          +unz(ia,1,iface,iel)*vz(ix,iy,iz,iel))
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end if

      end subroutine


c-----------------------------------------------------------------------
      subroutine vol_flowp
              ! HARRY
              ! you've just ripped off vol_flow but for the pertubation
              ! solver now
              ! very untest.... probably full of bugs
c
c
c     Adust flow volume at end of time step to keep flow rate fixed by
c     adding an appropriate multiple of the linear solution to the Stokes
c     problem arising from a unit forcing in the X-direction.  This assumes
c     that the flow rate in the X-direction is to be fixed (as opposed to Y-
c     or Z-) *and* that the periodic boundary conditions in the X-direction
c     occur at the extreme left and right ends of the mesh.
c
c     pff 6/28/98
c
      include 'SIZE'
      include 'TOTAL'
c
c     Swap the comments on these two lines if you don't want to fix the
c     flow rate for periodic-in-X (or Z) flow problems.
c
      parameter (kx1=lx1,ky1=ly1,kz1=lz1,kx2=lx2,ky2=ly2,kz2=lz2)
c
      common /cvflow_ap/ vxc(kx1,ky1,kz1,lelv)
     $                , vyc(kx1,ky1,kz1,lelv)
     $                , vzc(kx1,ky1,kz1,lelv)
     $                , prc(kx2,ky2,kz2,lelv)
     $                , vdc(kx1*ky1*kz1*lelv,2)
      common /cvflow_rp/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)
      common /cvflow_ip/ icvflow,iavflow
      common /cvflow_cp/ chv(3)
      character*1 chv

      real bd_vflow,dt_vflow
      save bd_vflow
      data bd_vflow,dt_vflow /-99.,-99./
      real legit_volflow

      logical ifcomp
      logical flow_angle      ! if we want an angle
      real x_vec(3)           ! a UNIT vector pointing in that direction
      common /harry_volflow_l/ flow_angle
      real tmp_scale
      real x_len, y_len, z_len
      common /harry_volflow_r/ x_vec, x_len, y_len, z_len

c     Check list:

c     param (55) -- volume flow rate, if nonzero
c     forcing in X? or in Z?




      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

      ! HARRY we WANT zero mass flow rate!
      !if (param(55).eq.0.) return
      if (param(55).eq.0.0.and.param(54).eq.0.0) return

      ! HARRY make sure it's a unit vector
      if(flow_angle) then
      if(abs(x_vec(1)**2 + x_vec(2)**2 + x_vec(3)**2 - 1.0).gt.1e-9)
     $      then
      if(nio.eq.0) write(6,*) 'VOLFLOW: rescaling flow vector'
      if(nio.eq.0) write(6,*) 'Before: ', x_vec(1),x_vec(2),x_vec(3)

      tmp_scale = sqrt(x_vec(1)**2 + x_vec(2)**2 + x_vec(3)**2) 
      call cmult(x_vec,1.0/tmp_scale,3)
      if(nio.eq.0) write(6,*) 'After:  ', x_vec(1),x_vec(2),x_vec(3)
      endif
      endif
      if (kx1.eq.1) then
         write(6,*) 'ABORT. Recompile vol_flow with kx1=lx1, etc.'
         call exitt
      endif

      icvflow   = 1                                  ! Default flow dir. = X
      if (param(54).ne.0) icvflow = abs(param(54))
      iavflow   = 0                                  ! Determine flow rate
      if (param(54).lt.0) iavflow = 1                ! from mean velocity
      flow_rate = param(55)

      chv(1) = 'X'
      chv(2) = 'Y'
      chv(3) = 'Z'

c     If either dt or the backwards difference coefficient change,
c     then recompute base flow solution corresponding to unit forcing:

      ifcomp = .false.
      if (dt.ne.dt_vflow.or.bd(1).ne.bd_vflow.or.ifmvbd) ifcomp=.true.
      if (.not.ifcomp) then
         ifcomp=.true.
         do i=1,ntot1
            if (vdiff (i,1,1,1,1).ne.vdc(i,1)) goto 20
            if (vtrans(i,1,1,1,1).ne.vdc(i,2)) goto 20
         enddo
         ifcomp=.false.  ! If here, then vdiff/vtrans unchanged.
   20    continue
      endif
      call gllog(ifcomp,.true.)
      
      call copy(vdc(1,1),vdiff (1,1,1,1,1),ntot1)
      call copy(vdc(1,2),vtrans(1,1,1,1,1),ntot1)
      dt_vflow = dt
      bd_vflow = bd(1)

      if (ifcomp) call compute_vol_solnp(vxc,vyc,vzc,prc)

      if(flow_angle) then
      ! HARRY
      ! here is if we're at an angle
      ! we want int_\Omega U . x_vc
      ! NOTE, we will NOT divide by the length now!! we do the whole
      ! volume integral
      current_flow = glsc2(vxp,bm1,ntot1)*x_vec(1)  
     $             + glsc2(vyp,bm1,ntot1)*x_vec(2)  
     $             + glsc2(vzp,bm1,ntot1)*x_vec(3)  
      if (iavflow.eq.1) then
         ! again, whole volume. this should be C_flux
         flow_rate = param(55)*volvm1
      endif
      else
      if (icvflow.eq.1) current_flow=glsc2(vxp,bm1,ntot1)/domain_length 
      if (icvflow.eq.2) current_flow=glsc2(vyp,bm1,ntot1)/domain_length 
      if (icvflow.eq.3) current_flow=glsc2(vzp,bm1,ntot1)/domain_length 

      if (iavflow.eq.1) then
         xsec = volvm1 / domain_length
         flow_rate = param(55)*xsec
      endif
      endif

      delta_flow = flow_rate-current_flow

c     Note, this scale factor corresponds to FFX, provided FFX has
c     not also been specified in userf.   If ffx is also specified
c     in userf then the true FFX is given by ffx_userf + scale.

      ! HARRY
      ! hopefully base_flow is ALSO modified to be a VOLUME if we're
      ! dealing wth flow_angle
      scale = delta_flow/base_flow
      saved_scale = scale
      scale_vf(icvflow) = scale
      if (nio.eq.0) write(6,1) istep,chv(icvflow)
     $   ,time,scale,delta_flow,current_flow,flow_rate
    1    format(i11,'  Volflow ',a1,11x,1p5e13.4)

      call add2s2(vxp,vxc,scale,ntot1)
      call add2s2(vyp,vyc,scale,ntot1)
      call add2s2(vzp,vzc,scale,ntot1)
      call add2s2(prp,prc,scale,ntot2)

      ! HARRY
      ! I'm suspicious that with inhomegenious BC's we're not going to
      ! hit the mass_flowRate on the dot
      if (icvflow.eq.1) legit_volflow=glsc2(vxp,bm1,ntot1)/volvm1 
      if (icvflow.eq.2) legit_volflow=glsc2(vyp,bm1,ntot1)/volvm1 
      if (icvflow.eq.3) legit_volflow=glsc2(vzp,bm1,ntot1)/volvm1 
      if(nio.eq.0) write(6,*) 'HARRY legit volflow = ',legit_volflow
      if(nio.eq.0) write(6,*) 'HARRY scale  = ',scale

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_vol_solnp(vxc,vyc,vzc,prc)
c
c     Compute the solution to the time-dependent Stokes problem
c     with unit forcing, and find associated flow rate.
c
c     pff 2/28/98
c
      include 'SIZE'
      include 'TOTAL'
c
      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)
c
      common /cvflow_rp/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)
      common /cvflow_ip/ icvflow,iavflow
      common /cvflow_cp/ chv(3)
      character*1 chv
c
      integer icalld
      save    icalld
      data    icalld/0/
      logical flow_angle      ! if we want an angle
      real x_vec(3)           ! a vector pointing in that direction
      real x_len, y_len, z_len
      common /harry_volflow_l/ flow_angle
      common /harry_volflow_r/ x_vec, x_len, y_len, z_len
c
c
      ntot1 = lx1*ly1*lz1*nelv
      if (icalld.eq.0) then
         icalld=icalld+1
         xlmin = glmin(xm1,ntot1)
         xlmax = glmax(xm1,ntot1)
         ylmin = glmin(ym1,ntot1)          !  for Y!
         ylmax = glmax(ym1,ntot1)
         zlmin = glmin(zm1,ntot1)          !  for Z!
         zlmax = glmax(zm1,ntot1)
c
         ! HARRY
         ! so I can have these later
         x_len = xlmax - xlmin
         y_len = ylmax - ylmin
         z_len = zlmax - zlmin
         if (icvflow.eq.1) domain_length = xlmax - xlmin
         if (icvflow.eq.2) domain_length = ylmax - ylmin
         if (icvflow.eq.3) domain_length = zlmax - zlmin
c
      endif
c
      if (ifsplit) then
c        call plan2_vol(vxc,vyc,vzc,prc)
         call plan4_volp(vxc,vyc,vzc,prc)
      else
         call plan3_volp(vxc,vyc,vzc,prc)
      endif
c
c     Compute base flow rate
c 
      if(flow_angle) then
      ! AGAIN, now we have volume integrals
      ! this is int_\Omega U' . x_vec
      base_flow = glsc2(vxc,bm1,ntot1)*x_vec(1)  
     $             + glsc2(vyc,bm1,ntot1)*x_vec(2)  
     $             + glsc2(vzc,bm1,ntot1)*x_vec(3)  
      else
      if (icvflow.eq.1) base_flow = glsc2(vxc,bm1,ntot1)/domain_length
      if (icvflow.eq.2) base_flow = glsc2(vyc,bm1,ntot1)/domain_length
      if (icvflow.eq.3) base_flow = glsc2(vzc,bm1,ntot1)/domain_length
      endif
c
      if (nio.eq.0 .and. loglevel.gt.2) write(6,1) 
     $   istep,chv(icvflow),base_flow,domain_length,flow_rate
    1    format(i11,'  basflowp ',a1,11x,1p3e13.4)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine plan2_volp(vxc,vyc,vzc,prc)
c
c     Compute pressure and velocity using fractional step method.
c     (classical splitting scheme).
c
c
      ! HARRY
      ! I'm not using this one
      return
      end
c-----------------------------------------------------------------------
      subroutine plan3_volp(vxc,vyc,vzc,prc)
c
c     Compute pressure and velocity using fractional step method.
c     (PLAN3).
c
c
      include 'SIZE'
      include 'TOTAL'

c
      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)
C
      logical flow_angle      ! if we want an angle
      real x_vec(3)           ! a vector pointing in that direction
      common /harry_volflow_l/ flow_angle
      common /harry_volflow_r/ x_vec

      COMMON /SCRNS/ rw1   (LX1,LY1,LZ1,LELV)
     $ ,             rw2   (LX1,LY1,LZ1,LELV)
     $ ,             rw3   (LX1,LY1,LZ1,LELV)
     $ ,             dv1   (LX1,LY1,LZ1,LELV)
     $ ,             dv2   (LX1,LY1,LZ1,LELV)
     $ ,             dv3   (LX1,LY1,LZ1,LELV)
     $ ,             RESPR (LX2,LY2,LZ2,LELV)
      COMMON /SCRVH/ H1    (LX1,LY1,LZ1,LELV)
     $ ,             H2    (LX1,LY1,LZ1,LELV)
      COMMON /SCRHI/ H2INV (LX1,LY1,LZ1,LELV)
      common /cvflow_ip/ icvflow,iavflow

c
c
c     Compute velocity, 1st part 
c
      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      ifield = 1
c
      if(flow_angle) then
      ! HARRY
      ! this is if we have x_vec defining a flow vector
      ! I guess this assumes x_vec is UNIT vector
         call copy     (rw1,bm1,ntot1)
         call copy     (rw2,bm1,ntot1)
         call copy     (rw3,bm1,ntot1)
         call cmult    (rw1,x_vec(1),ntot1)
         call cmult    (rw2,x_vec(2),ntot1)
         call cmult    (rw3,x_vec(3),ntot1)
      else
      if (icvflow.eq.1) then
         call copy     (rw1,bm1,ntot1)
         call rzero    (rw2,ntot1)
         call rzero    (rw3,ntot1)
      elseif (icvflow.eq.2) then
         call rzero    (rw1,ntot1)
         call copy     (rw2,bm1,ntot1)
         call rzero    (rw3,ntot1)
      else
         call rzero    (rw1,ntot1)        ! Z-flow!
         call rzero    (rw2,ntot1)        ! Z-flow!
         call copy     (rw3,bm1,ntot1)    ! Z-flow!
      endif
      endif
      intype = -1
      ! HARRY
      ! this is the sketchiest part.
      
      ! This was fluidp
      !call cresvipp (resv1,resv2,resv3,h1,h2)
      !call ophinv   (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
      !call opadd2   (vxp(1,jp),vyp(1,jp),vzp(1,jp),dv1,dv2,dv3)
      !call incomprp (vxp(1,jp),vyp(1,jp),vzp(1,jp),prp(1,jp))

      ! I'm trying to base this off fluidp and incomprnp
      call sethlm   (h1,h2,intype) ! I think this should be the same
      ! also this
      ! I'm getting rid of this... because it's for inhomogenious BCs
      !call cresvipp (resv1,resv2,resv3,h1,h2)
      call ophinv   (vxc,vyc,vzc,rw1,rw2,rw3,h1,h2,tolhv,nmxv)
      ! now fluidp calls opadd2(vxp,vyp...
      ! Following vol_flow here (wtf does this call even do??)
      call ssnormd  (vxc,vyc,vzc)
c
c     Compute pressure  (from "incompr")
c
      intype = 1
      dtinv  = 1./dt
c
      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call cmult   (h2,dtinv,ntot1)
      call invers2 (h2inv,h2,ntot1)
      call opdiv   (respr,vxc,vyc,vzc)
      call chsign  (respr,ntot2)
      call ortho   (respr)
      ! so far this is the same as vol_flow and incomprnp
c
c
c     Set istep=0 so that h1/h2 will be re-initialized in eprec
      i_tmp = istep
      istep = 0
      call esolver (respr,h1,h2,h2inv,intype)
      istep = i_tmp
c
      call opgradt (rw1,rw2,rw3,respr)
      call opbinv  (dv1,dv2,dv3,rw1,rw2,rw3,h2inv)
      call opadd2  (vxc,vyc,vzc,dv1,dv2,dv3)
      ! here is where incomprnp continues with lag terms etc
      ! they call extrapprp(prextr) first, then they add
      ! I'm worried this is not enough...
      ! it's too similar to the original volflow
c
      call cmult2  (prc,respr,bd(1),ntot2)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine plan4_volp(vxc,vyc,vzc,prc)

c     Compute pressure and velocity using fractional step method.
c     (Tombo splitting scheme).
      ! HARRY
      ! or this one

      end
c-----------------------------------------------------------------------


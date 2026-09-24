!
!     F-barES-FEM-T4: internal forces for U3, the F-barES-FEM-T4 volumetric domain.
!
!     f_(i,c) = (K V_h) * tbar(c,i) * thetabar,   thetabar = sbar . u
!
!     i.e. the rank-one operator e_c3d_u3 assembles, contracted the only way
!     round that is consistent with it: the smoothed divergence thetabar comes
!     from the TRIAL row sbar (eq. 6-11, the c-time cyclically smoothed J), the
!     nodal force is spread by the TEST row tbar (eq. 1, unsmoothed).  Swapping
!     them recovers the transpose -- a different element, and the sign that the
!     Petrov-Galerkin structure of eq. (17) has been dropped.
!
!     In the symmetric (Galerkin) pairing, marked by '_' in lakon(3:3), the
!     stiffness is xkv*sbar^T*sbar, so the consistent internal force spreads
!     the same thetabar by sbar, not tbar.  THE TWO MUST USE THE SAME PAIRING:
!     the first results call assembles f = K*u_prescribed into the right hand
!     side, and a tbar spread there was silently dropping the symmetric
!     coupling -- measured as H1 1.0286 against 0.479 at n=4, K/G=100, with
!     the assembled matrix itself exact.
!
!     Same u3vol call as the stiffness, so the two cannot drift apart; see
!     resultsmech_u2 for why that matters.  No stress is written.
!
      subroutine resultsmech_u3(co,kon,ipkon,lakon,ne,v,elcon,nelcon,
     &     ielmat,mi,ncmat_,ntmat_,fn,calcul_fn,nal,nelem,iperturb)
!
      implicit none
!
      character*8 lakon(*)
      integer mi(*)
      integer kon(*),ipkon(*),ielmat(mi(3),*),ncmat_,ntmat_,
     &     nelcon(2,*),nelem,ne,calcul_fn,nal,i,c,nope,indexe,
     &     konl(255),imat,ncyc,ilen,isym,iperturb(*)
      real*8 co(3,*),v(0:mi(2),*),elcon(0:ncmat_,ntmat_,*),
     &     fn(0:mi(2),*),vh,xkv,sbar(3,255),tbar(3,255),theta,fac
      character*16 cval
!
      nope=ichar(lakon(nelem)(8:8))
      indexe=ipkon(nelem)
      do i=1,nope
        konl(i)=kon(indexe+i)
      enddo
!
!     NLGEOM: the finite-strain force of eqs. (1)-(17) through the
!     shared kernel, which builds the same E/S chain u3vol does.
!
      if(iperturb(1).gt.1) then
        if(calcul_fn.eq.1) then
          imat=ielmat(1,nelem)
          ncyc=1
          call getenv('CCX_FBAR_C',cval)
          ilen=len_trim(cval)
          if(ilen.gt.0) read(cval(1:ilen),*) ncyc
          call u3nlforce(co,kon,ipkon,lakon,ne,konl,nope,v,elcon,
     &         nelcon,ielmat,mi,ncmat_,ntmat_,imat,ncyc,fn,nelem)
        endif
        nal=nal+3*nope
        return
      endif
!
!     the '_' marker selects the symmetric (Galerkin) pairing, exactly as
!     e_c3d_u3 and mafillsmas.f read it.
!
      isym=0
      if(lakon(nelem)(3:3).eq.'_') isym=1
!
!     the same CCX_FBAR_C the stiffness read: a mismatch here is precisely the
!     stiffness/force drift this file exists to prevent.
!
      ncyc=1
      call getenv('CCX_FBAR_C',cval)
      ilen=len_trim(cval)
      if(ilen.gt.0) read(cval(1:ilen),*) ncyc
      if((ncyc.lt.0).or.(ncyc.gt.3)) then
        write(*,*) '*ERROR in resultsmech_u3: CCX_FBAR_C =',ncyc
        write(*,*) '       must be 0..3'
        call exit(201)
      endif
!
      imat=ielmat(1,nelem)
      call u3vol(co,kon,ipkon,lakon,ne,konl,nope,vh,sbar,tbar,xkv,
     &     nelem,ielmat,elcon,nelcon,mi,ncmat_,ntmat_,imat,ncyc)
      if(vh.le.0.d0) return
!
      theta=0.d0
      do i=1,nope
        do c=1,3
          theta=theta+sbar(c,i)*v(c,konl(i))
        enddo
      enddo
      fac=xkv*theta
!
      if(calcul_fn.eq.1) then
        if(isym.eq.0) then
          do i=1,nope
            do c=1,3
              fn(c,konl(i))=fn(c,konl(i))+fac*tbar(c,i)
            enddo
          enddo
        else
          do i=1,nope
            do c=1,3
              fn(c,konl(i))=fn(c,konl(i))+fac*sbar(c,i)
            enddo
          enddo
        endif
      endif
!
      nal=nal+3*nope
!
      return
      end

!
!     F-barES-FEM-T4: internal forces for U2, the ES-FEM deviatoric smoothing domain.
!
!     f_(i,c) = V_h * sum_a gt(a,i) * sig_(c,a),  sig = 2 mu dev(sym(grad u))
!
!     which is exactly sum_(j,d) s(ii,jj) u_(j,d) for the s that e_c3d_u2
!     assembles -- expanded once on paper so the recovery is one pass over the
!     ring instead of a (3*nope)^2 product.  Both routines get their geometry
!     from the SAME u2edge call, which is the point: patch 0002 records what
!     happens when the stiffness is patched and the forces are not, and the
!     reaction cross-check -- the only reference-free convergence evidence this
!     repository has -- stops meaning anything.
!
!     No stress is written.  The smoothing domain is not an element with a
!     shape function, so stx is left alone and U2 must stay out of every
!     *EL PRINT set.  The base tets (U4) are null too, so an F-barES-FEM-T4
!     run carries no element stress at all: read the answer from
!     displacements and reactions.
!
      subroutine resultsmech_u2(co,kon,ipkon,lakon,ne,v,elcon,nelcon,
     &     ielmat,mi,ncmat_,ntmat_,fn,calcul_fn,nal,nelem,iperturb)
!
      implicit none
!
      character*8 lakon(*)
      integer mi(*)
      integer kon(*),ipkon(*),ielmat(mi(3),*),ncmat_,ntmat_,
     &     nelcon(2,*),nelem,ne,calcul_fn,nal,i,a,c,nope,indexe,
     &     konl(255),imat,nfound,iperturb(*)
      real*8 co(3,*),v(0:mi(2),*),elcon(0:ncmat_,ntmat_,*),
     &     fn(0:mi(2),*),vh,gt(3,255),e,un,um,g(3,3),sig(3,3),tr
      character*16 cval
!
      nope=ichar(lakon(nelem)(8:8))
      indexe=ipkon(nelem)
      do i=1,nope
        konl(i)=kon(indexe+i)
      enddo
!
!     NLGEOM: U2 contributes NO separate force.  The finite-strain force
!     of eq. (17) evaluates the full Hencky stress of Fbar and pairs it
!     with the ES-FEM motion, so its linearisation IS U2 + U3; adding the
!     small-strain deviatoric force on top of it would count the
!     deviatoric half twice.  The stiffness split is unaffected: U2 and U3
!     remain the two halves of df/du at the reference configuration.
!
      if(iperturb(1).gt.1) then
        nal=nal+3*nope
        return
      endif
!
      imat=ielmat(1,nelem)
      e=elcon(1,1,imat)
      un=elcon(2,1,imat)
      um=e/(2.d0*(1.d0+un))
!
      call u2edge(co,kon,ipkon,lakon,ne,konl,nope,vh,gt,nelem,
     &     ielmat,mi,imat,nfound)
      if(vh.le.0.d0) return
!
!     g(a,c) = du_a/dx_c, smoothed over the edge domain
!
      do a=1,3
        do c=1,3
          g(a,c)=0.d0
        enddo
      enddo
      do i=1,nope
        do a=1,3
          do c=1,3
            g(a,c)=g(a,c)+gt(c,i)*v(a,konl(i))
          enddo
        enddo
      enddo
      tr=g(1,1)+g(2,2)+g(3,3)
!
      do c=1,3
        do a=1,3
          sig(c,a)=um*(g(c,a)+g(a,c))
        enddo
        sig(c,c)=sig(c,c)-2.d0*um/3.d0*tr
      enddo
!
      if(calcul_fn.eq.1) then
        do i=1,nope
          do c=1,3
            do a=1,3
              fn(c,konl(i))=fn(c,konl(i))+vh*gt(a,i)*sig(c,a)
            enddo
          enddo
        enddo
      endif
!
      nal=nal+3*nope
!
      return
      end

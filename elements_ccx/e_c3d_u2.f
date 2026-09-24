!
!     F-barES-FEM-T4: U2 -- edge-based smoothed DEVIATORIC stiffness (ES-FEM-T4).
!
!     The deviatoric half of F-barES-FEM-T4 [Onishi, Iida & Amaya, Int. J.
!     Comput. Methods 15(7) 1845003 (2018)].  One element per edge smoothing
!     domain, per material.  Its stiffness is eqs. (1), (4) and (13):
!
!         K_dev^h = V_h * Bt^T D_dev Bt
!
!     with Bt the edge-smoothed strain-displacement operator of eq. (1) and
!     D_dev = 2 mu (T - m m^T / 3).  Deriving from eqs. (4) and (11): F^iso
!     depends only on Ftilde, and Fbar differs from Ftilde^iso by a spherical
!     factor, so dev(H) never sees the cyclic smoothing of J -- the deviatoric
!     response is pure ES-FEM and this element is the whole of it.
!
!     Together with U3 (the volumetric half, eqs. 6-11) this replaces the base
!     tets' stiffness entirely, so the base tets must contribute NOTHING: the
!     generator retypes them to U4, the null base tetrahedron, which keeps
!     them in the model as the geometry U2 and U3 read.
!
!         *USER ELEMENT,TYPE=U2,NODES=<ring size>,INTEGRATIONPOINTS=1,MAXDOF=3
!
!     konl(1) and konl(2) ARE THE TWO EDGE NODES.  The rest is the other nodes
!     of the tets sharing that edge; u2edge finds the tets themselves.  On the
!     campaign's LMESH_m0p0240 cell this ring is 6.3 nodes on average and 14 at
!     worst, so unlike U3 it sits far inside the element-matrix capacity.
!
!     Keep U2 out of any *EL PRINT set: it has no material volume of its own,
!     and printoutelem.f would pick a shape function from its node count.
!
      subroutine e_c3d_u2(co,kon,lakonl,s,sm,ff,nelem,elcon,nelcon,
     &     ielmat,mi,ncmat_,ntmat_,ipkon,lakon,ne,stiffness)
!
      implicit none
!
      character*8 lakonl,lakon(*)
      integer mi(*)
      integer kon(*),ipkon(*),ielmat(mi(3),*),ncmat_,ntmat_,
     &     nelcon(2,*),nelem,ne,stiffness,i,j,c,d,ii,jj,nope,indexe,
     &     imat,konl(255),nfound
      real*8 co(3,*),s(765,765),sm(765,765),ff(765),
     &     elcon(0:ncmat_,ntmat_,*),
     &     e,un,um,vh,gt(3,255),fac
!
      nope=ichar(lakonl(8:8))
!
!     Hard capacity guard, as in e_c3d_u6: nope follows mesh valence, and
!     overrunning s used to corrupt the assembly silently rather than stop.
!     255 nodes (765 DOF) is the encoding limit -- lakon(8:8) is one byte and
!     userelements.f:83 already refuses more -- so a ring that trips this
!     cannot be a ccx user element at all.
!
      if(nope.gt.255) then
        write(*,*) '*ERROR in e_c3d_u2: edge element ',nelem
        write(*,*) '       spans ',nope,' nodes (',3*nope,' DOF).'
        write(*,*) '       The limit is 255 nodes (765 DOF) -- the most'
        write(*,*) '       that lakon(8:8) can encode.'
        call exit(201)
      endif
      indexe=ipkon(nelem)
      do i=1,nope
        konl(i)=kon(indexe+i)
      enddo
!
!
!     Zero the ACTUAL extent, and zero sm too: the caller allocates s/sm/ff
!     once and never re-zeroes them between elements, and U2 is static-only
!     and never writes sm (patch 0008).
!
      do i=1,3*nope
        ff(i)=0.d0
        do j=1,3*nope
          s(i,j)=0.d0
          sm(i,j)=0.d0
        enddo
      enddo
      if(stiffness.eq.0) return
!
      imat=ielmat(1,nelem)
      if(nelcon(1,imat).ne.2) then
        write(*,*) '*ERROR in e_c3d_u2: element',nelem,' needs an'
        write(*,*) '       isotropic *ELASTIC card (2 constants)'
        call exit(201)
      endif
      e=elcon(1,1,imat)
      un=elcon(2,1,imat)
      um=e/(2.d0*(1.d0+un))
!
      call u2edge(co,kon,ipkon,lakon,ne,konl,nope,vh,gt,nelem,
     &     ielmat,mi,imat,nfound)
      if(vh.le.0.d0) return
      if(nfound.eq.0) then
        write(*,*) '*ERROR in e_c3d_u2: element',nelem,' found no tet'
        write(*,*) '       containing its edge ',konl(1),'-',konl(2)
        call exit(201)
      endif
!
!     2 mu dev(eps):dev(eps) over the smoothing domain.  Identical in form to
!     e_c3d_u5 with vol -> vh and g -> gt:
!         eps:eps = (delta_cd (g_i.g_j) + g_i[d] g_j[c])/2
!         div div = g_i[c] g_j[d]
!
      do i=1,nope
        do c=1,3
          ii=3*(i-1)+c
          do j=1,nope
            do d=1,3
              jj=3*(j-1)+d
              fac=0.d0
              if(c.eq.d) fac=gt(1,i)*gt(1,j)+gt(2,i)*gt(2,j)
     &             +gt(3,i)*gt(3,j)
              s(ii,jj)=vh*(um*(fac+gt(d,i)*gt(c,j))
     &             -2.d0*um/3.d0*gt(c,i)*gt(d,j))
            enddo
          enddo
        enddo
      enddo
!
      return
      end

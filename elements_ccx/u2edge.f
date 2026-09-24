!
!     F-barES-FEM-T4: U2 -- the edge-based (ES-FEM) DEVIATORIC smoothing domain.
!
!     U2 is the deviatoric and U3 the volumetric smoothing domain of the
!     F-barES-FEM-T4 element; U4 is the null base tetrahedron whose only
!     role is to carry the geometry.
!
!     Geometry shared by e_c3d_u2 (stiffness) and resultsmech_u2 (internal
!     forces), so the two cannot drift apart -- the failure that leaves a patch
!     element reporting an equilibrium_gap of 1.0.
!
!     This is eq. (1) of Onishi, Iida & Amaya, IJCM 15(7) 1845003 (2018),
!     restricted to the deviatoric part.  For edge h:
!
!         V_h        = sum over base tets containing BOTH edge nodes of V_e/6
!         gt(c,i)    = (1/V_h) sum over those tets of (V_e/6) grad(L_i)|_e
!
!     and the stiffness is V_h * Bt^T D_dev Bt, which in index form is exactly
!     the C3D4 deviatoric expression with vol -> V_h and g -> gt.  The /6 is the edge
!     count of a T4, so each tet hands its whole volume to its six edges and
!     sum_h V_h = sum_e V_e exactly.
!
!     THE EDGE IS konl(1)-konl(2).  The rest of the connectivity is the other
!     nodes of the tets sharing it, in any order.  The geometry is not
!     recoverable from the node list alone -- which subsets of the ring form
!     tets is not determined by it -- so the element finds its own tets through
!     a node->element map over the base 'U4' elements, built once and cached
!     under a lock.
!
!     SMOOTHING NEVER CROSSES A MATERIAL INTERFACE: only tets of the edge's
!     own material contribute.  Averaging the deviatoric
!     response across a 1000x modulus contrast is not the method, and an
!     interface edge legitimately gets one smoothing domain per phase, each
!     with its own V_h.  Summed over both, sum_h V_h = sum_e V_e still holds.
!
      subroutine u2edge(co,kon,ipkon,lakon,ne,konl,nope,vh,gt,
     &     nelem,ielmat,mi,imatf,nfound)
!
      implicit none
!
      character*8 lakon(*)
      integer mi(*)
      integer kon(*),ipkon(*),ne,konl(*),nope,nelem,nfound,
     &     i,j,k,c,ie,ipos,n4(4),ielmat(mi(3),*),imatf,na,nb,
     &     ihit1,ihit2
      real*8 co(3,*),vh,gt(3,*),xl(3,4),shp(4,4),xsj,vol,w
!
      integer mapdone,maxn,nlen
      integer, allocatable, save :: nstart(:),nlist(:)
      save mapdone,maxn
      data mapdone /0/
!
      call fbarlock()
      if(mapdone.eq.0) then
        maxn=0
        nlen=0
        do i=1,ne
          if(ipkon(i).lt.0) cycle
          if(lakon(i)(1:2).ne.'U4') cycle
          do j=1,4
            k=kon(ipkon(i)+j)
            if(k.gt.maxn) maxn=k
            nlen=nlen+1
          enddo
        enddo
        allocate(nstart(maxn+2))
        allocate(nlist(max(nlen,1)))
        do i=1,maxn+2
          nstart(i)=0
        enddo
        do i=1,ne
          if(ipkon(i).lt.0) cycle
          if(lakon(i)(1:2).ne.'U4') cycle
          do j=1,4
            nstart(kon(ipkon(i)+j)+1)=nstart(kon(ipkon(i)+j)+1)+1
          enddo
        enddo
        do i=2,maxn+2
          nstart(i)=nstart(i)+nstart(i-1)
        enddo
        do i=1,ne
          if(ipkon(i).lt.0) cycle
          if(lakon(i)(1:2).ne.'U4') cycle
          do j=1,4
            k=kon(ipkon(i)+j)
            nstart(k)=nstart(k)+1
            nlist(nstart(k))=i
          enddo
        enddo
        do i=maxn+1,2,-1
          nstart(i)=nstart(i-1)
        enddo
        nstart(1)=0
        mapdone=1
      endif
      call fbarunlock()
!
      vh=0.d0
      nfound=0
      do i=1,nope
        do c=1,3
          gt(c,i)=0.d0
        enddo
      enddo
!
      na=konl(1)
      nb=konl(2)
      if((na.gt.maxn).or.(nb.gt.maxn)) then
        write(*,*) '*ERROR in u2edge: element',nelem,' edge node'
        write(*,*) '       ',na,' or ',nb,' is in no U4 element'
        call exit(201)
      endif
!
!     Walk the tets at the FIRST edge node and keep those that also contain
!     the second.  Scanning one node's list is enough -- a tet containing the
!     edge contains both its nodes -- and it is the shorter loop.
!
      do ie=nstart(na)+1,nstart(na+1)
        i=nlist(ie)
        if(ielmat(1,i).ne.imatf) cycle
        ihit1=0
        ihit2=0
        do j=1,4
          n4(j)=kon(ipkon(i)+j)
          if(n4(j).eq.na) ihit1=1
          if(n4(j).eq.nb) ihit2=1
        enddo
        if((ihit1.eq.0).or.(ihit2.eq.0)) cycle
        nfound=nfound+1
        do j=1,4
          do k=1,3
            xl(k,j)=co(k,n4(j))
          enddo
        enddo
        call shape4tet(0.25d0,0.25d0,0.25d0,xl,xsj,shp,3)
        vol=dabs(xsj)/6.d0
        w=vol/6.d0
        vh=vh+w
        do j=1,4
          ipos=0
          do k=1,nope
            if(konl(k).eq.n4(j)) then
              ipos=k
              exit
            endif
          enddo
          if(ipos.eq.0) then
            write(*,*) '*ERROR in u2edge: element',nelem,' node',
     &           n4(j),' not in its own connectivity'
            call exit(201)
          endif
          do c=1,3
            gt(c,ipos)=gt(c,ipos)+w*shp(c,j)
          enddo
        enddo
      enddo
!
      if(vh.gt.0.d0) then
        do i=1,nope
          do c=1,3
            gt(c,i)=gt(c,i)/vh
          enddo
        enddo
      endif
!
      return
      end

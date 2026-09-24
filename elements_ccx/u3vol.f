!
!     F-barES-FEM-T4: U3 -- the F-barES-FEM-T4 VOLUMETRIC smoothing chain, eqs. (6)-(11).
!
!     U3 is the volumetric smoothing domain of the F-barES-FEM-T4 element
!     (U2 the deviatoric, U4 the null base tetrahedron).
!
!     Geometry and operator shared by e_c3d_u3 (stiffness) and resultsmech_u3
!     (internal forces), so the two cannot drift apart.
!
!     For edge h this builds the row of  S = E A^c  over elements, where
!
!         (Q x)_n = (1/V_n) sum_{e at n} x_e V_e/4          eq. (6)
!         (P y)_e = (1/4) sum_{n in e} y_n                  eq. (7)
!         A       = P Q                                     eq. (7) closes a cycle
!         (E z)_h = (1/V_h) sum_{e at h} z_e V_e/6          eq. (8)
!
!     and then contracts it with the element divergence to give the nodal
!     coefficients of theta_bar_h:
!
!         theta_bar_h = sum_a sbar(:,a) . u_a ,
!         sbar(c,a)   = sum_e (E A^c)_he grad(L_a)|_e
!
!     It also returns the UNSMOOTHED edge row, tbar = (E D_div)_h, because
!     eq. (17) pairs the stress with the stretching of Ftilde, not of Fbar:
!     the volumetric operator is tbar^T (K V_h) sbar and is NOT symmetric.
!     See sections 2f and 4 of the companion paper.
!
!     Eq. (10) of the paper states outright that the whole chain is a weighted
!     mean of the RAW element J, which is what makes it this single linear
!     operator.  The input to eq. (6) carries no tilde -- read that off the
!     typeset equations, not a text extraction; pdftotext drops every tilde and
!     bar in this paper and the diacritics carry the meaning.
!
!     STENCIL WIDTH.  Measured on the campaign's LMESH_m0p0240 soft phase
!     (117437 tets, 36323 nodes, 184572 edges):
!
!         c = 1   mean  33.7 nodes/edge   p99  78   max 173  ( 519 DOF)
!         c = 2   mean  93.6 nodes/edge   p99 259   max 494  (1482 DOF)
!
!     c = 1 is deliverable as an element and c = 2 IS NOT: userelements.f:83
!     rejects NODES > 255, and mastruct.c reads the count from a single
!     character of the label.  c = 2 needs a direct global assembly pass
!     instead of the element abstraction.
!
      subroutine u3vol(co,kon,ipkon,lakon,ne,konl,nope,vh,sbar,tbar,
     &     xkv,nelem,ielmat,elcon,nelcon,mi,ncmat_,ntmat_,imatf,ncyc)
!
      implicit none
!
      character*8 lakon(*)
      integer mi(*)
      integer kon(*),ipkon(*),ne,konl(*),nope,nelem,ncyc,
     &     i,j,k,c,ie,ipos,n4(4),ielmat(mi(3),*),nelcon(2,*),
     &     ncmat_,ntmat_,imatf,na,nb,ihit1,ihit2,ic,nl,ll,je,
     &     imat
      real*8 vn
      real*8 co(3,*),vh,sbar(3,*),tbar(3,*),xkv,xl(3,4),shp(4,4),
     &     xsj,vol,w,elcon(0:ncmat_,ntmat_,*),ek,eun,gsave(3,4),
     &     wgt
!
!     element weight lists for the chain: wcur is (E A^k)_h, wnxt its image
!     after one more cycle.  Sized for the measured worst case with headroom.
!
      integer maxel
      parameter(maxel=4096)
      integer elcur(maxel),elnxt(maxel),ncur,nnxt
      real*8 wcur(maxel),wnxt(maxel)
!
      integer mapdone,maxn,nlen
      integer, allocatable, save :: nstart(:),nlist(:)
      real*8, allocatable, save :: volel(:)
      save mapdone,maxn
      data mapdone /0/
!
!     The map, the nodal volumes V_n of eq. (6) and the element volumes are
!     built once and read by every thread.  The test MUST stay inside the
!     critical region: a double-checked flag outside it would let a
!     thread read a half-built map and silently corrupted the assembly.
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
        allocate(volel(ne))
        do i=1,maxn+2
          nstart(i)=0
        enddo
        do i=1,ne
          volel(i)=0.d0
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
!
!       element volumes.  V_n of eq. (6) is NOT cached: it is per material,
!       and a node on the brine/ice interface has a different V_n in each
!       phase.  Caching the plain sum over all elements at the node and then
!       skipping foreign ones in the walk would divide by an inflated V_n,
!       break the unit row sum of Q, and stop the chain preserving a constant
!       J -- i.e. fail the patch test, quietly and only at the interface.
!       It is a short loop over the node's own element list, so it is
!       recomputed where it is used.
!
        do i=1,ne
          if(ipkon(i).lt.0) cycle
          if(lakon(i)(1:2).ne.'U4') cycle
          do j=1,4
            do k=1,3
              xl(k,j)=co(k,kon(ipkon(i)+j))
            enddo
          enddo
          call shape4tet(0.25d0,0.25d0,0.25d0,xl,xsj,shp,3)
          volel(i)=dabs(xsj)/6.d0
        enddo
        mapdone=1
      endif
      call fbarunlock()
!
      vh=0.d0
      xkv=0.d0
      do i=1,nope
        do c=1,3
          sbar(c,i)=0.d0
          tbar(c,i)=0.d0
        enddo
      enddo
!
      na=konl(1)
      nb=konl(2)
      if((na.gt.maxn).or.(nb.gt.maxn)) then
        write(*,*) '*ERROR in u3vol: element',nelem,' edge node'
        write(*,*) '       ',na,' or ',nb,' is in no U4 element'
        call exit(201)
      endif
!
!     ---- eq. (8): the elements at edge h, weights V_e/6 / V_h ------------
!
      ncur=0
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
        ncur=ncur+1
        if(ncur.gt.maxel) then
          write(*,*) '*ERROR in u3vol: element',nelem,' edge ring'
          write(*,*) '       exceeds maxel =',maxel
          call exit(201)
        endif
        elcur(ncur)=i
        wcur(ncur)=volel(i)/6.d0
        vh=vh+volel(i)/6.d0
      enddo
      if(vh.le.0.d0) return
      do i=1,ncur
        wcur(i)=wcur(i)/vh
      enddo
!
!     K-weighted volume.  An edge domain never spans two materials here (the
!     smoothing is split by material), so K is that material's, but it is read
!     per contributing element: it makes a
!     mixed domain impossible to get silently wrong.
!
      do i=1,ncur
        imat=ielmat(1,elcur(i))
        ek=elcon(1,1,imat)
        eun=elcon(2,1,imat)
        xkv=xkv+ek/(3.d0*(1.d0-2.d0*eun))*wcur(i)*vh
      enddo
!
!     ---- tbar: the UNSMOOTHED edge divergence row, eq. (17)'s test space --
!
      do i=1,ncur
        call u3grad(co,kon,ipkon,elcur(i),gsave,n4)
        do j=1,4
          ipos=0
          do k=1,nope
            if(konl(k).eq.n4(j)) then
              ipos=k
              exit
            endif
          enddo
          if(ipos.eq.0) then
            write(*,*) '*ERROR in u3vol: element',nelem,' node',
     &           n4(j),' not in its own connectivity (tbar)'
            call exit(201)
          endif
          do c=1,3
            tbar(c,ipos)=tbar(c,ipos)+wcur(i)*gsave(c,j)
          enddo
        enddo
      enddo
!
!     ---- eqs. (6)-(7): c cycles of A = P Q, on the element weight list ----
!
      do ic=1,ncyc
        nnxt=0
        do i=1,ncur
          je=elcur(i)
          wgt=wcur(i)
          do j=1,4
            k=kon(ipkon(je)+j)
!           (P y)_e picks up 1/4 of each of its nodes; (Q x)_n spreads
!           x_e V_e/4 / V_n.  Composed and walked backwards from the edge,
!           element e' at node k receives wgt * (1/4) * V_e'/4 / V_n.
!           V_n of eq. (6), restricted to this material
            vn=0.d0
            do ll=nstart(k)+1,nstart(k+1)
              if(ielmat(1,nlist(ll)).ne.imatf) cycle
              vn=vn+volel(nlist(ll))/4.d0
            enddo
            if(vn.le.0.d0) cycle
            do ll=nstart(k)+1,nstart(k+1)
              nl=nlist(ll)
              if(ielmat(1,nl).ne.imatf) cycle
              ipos=0
              do ie=1,nnxt
                if(elnxt(ie).eq.nl) then
                  ipos=ie
                  exit
                endif
              enddo
              if(ipos.eq.0) then
                nnxt=nnxt+1
                if(nnxt.gt.maxel) then
                  write(*,*) '*ERROR in u3vol: element',nelem,
     &                 ' chain exceeds maxel =',maxel
                  call exit(201)
                endif
                elnxt(nnxt)=nl
                wnxt(nnxt)=0.d0
                ipos=nnxt
              endif
              wnxt(ipos)=wnxt(ipos)
     &             +wgt*0.25d0*volel(nl)/4.d0/vn
            enddo
          enddo
        enddo
        ncur=nnxt
        do i=1,ncur
          elcur(i)=elnxt(i)
          wcur(i)=wnxt(i)
        enddo
      enddo
!
!     ---- sbar: contract the chain with the element divergence -------------
!
      do i=1,ncur
        call u3grad(co,kon,ipkon,elcur(i),gsave,n4)
        do j=1,4
          ipos=0
          do k=1,nope
            if(konl(k).eq.n4(j)) then
              ipos=k
              exit
            endif
          enddo
          if(ipos.eq.0) then
            write(*,*) '*ERROR in u3vol: element',nelem,' node',
     &           n4(j),' not in its own connectivity (sbar).'
            write(*,*) '       The generator and the element disagree'
            write(*,*) '       about the c =',ncyc,' stencil.'
            call exit(201)
          endif
          do c=1,3
            sbar(c,ipos)=sbar(c,ipos)+wcur(i)*gsave(c,j)
          enddo
        enddo
      enddo
!
      return
      end
!
!
      subroutine u3grad(co,kon,ipkon,ie,g,n4)
!
!     Shape-function gradients and node list of one straight tet.
!
      implicit none
      integer kon(*),ipkon(*),ie,n4(4),j,k
      real*8 co(3,*),g(3,4),xl(3,4),shp(4,4),xsj
      do j=1,4
        n4(j)=kon(ipkon(ie)+j)
        do k=1,3
          xl(k,j)=co(k,n4(j))
        enddo
      enddo
      call shape4tet(0.25d0,0.25d0,0.25d0,xl,xsj,shp,3)
      do j=1,4
        do k=1,3
          g(k,j)=shp(k,j)
        enddo
      enddo
      return
      end

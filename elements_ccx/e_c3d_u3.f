!
!     F-barES-FEM-T4: U3 -- edge-based smoothed VOLUMETRIC stiffness, F-barES-FEM-T4.
!
!     The volumetric half of the method [Onishi, Iida & Amaya, Int. J. Comput.
!     Methods 15(7) 1845003 (2018)], eqs. (6)-(11) and (17).  One element per
!     edge smoothing domain, per material:
!
!         K_vol^h = (K V_h) * tbar^T sbar
!
!         sbar = (E A^c D_div)_h   trial: the c-time cyclically smoothed J
!         tbar = (E D_div)_h       test:  the UNSMOOTHED edge divergence
!
!     THIS MATRIX IS NOT SYMMETRIC, and that is the definition of F-bar, not
!     an approximation: the stress comes from the modified gradient while the
!     virtual work is paired with the unmodified one (eq. 17).  The Galerkin
!     form sbar^T sbar is a different element, within 0.02% of this one on
!     C1111.  The element therefore needs the asymmetric assembly path
!     (nasym=1, mafillsmas.f, PARDISO mtype=11).
!
!         *USER ELEMENT,TYPE=U3,NODES=<stencil>,INTEGRATIONPOINTS=1,MAXDOF=3
!
!     konl(1) and konl(2) ARE THE TWO EDGE NODES; the rest is the full E A^c
!     support, which the generator must compute with exactly the same walk
!     u3vol does -- u3vol stops with a message naming the element if it finds
!     a node the connectivity does not carry.
!
!     CAPACITY.  c = 1 reaches 173 nodes and fits the 255-node encoding
!     limit; c = 2 reaches 494 and cannot be a user element at all.
!
      subroutine e_c3d_u3(co,kon,lakonl,s,sm,ff,nelem,elcon,nelcon,
     &     ielmat,mi,ncmat_,ntmat_,ipkon,lakon,ne,stiffness,nasym,
     &     iperturb,vold)
!
      implicit none
!
      character*8 lakonl,lakon(*)
      integer mi(*)
      integer kon(*),ipkon(*),ielmat(mi(3),*),ncmat_,ntmat_,
     &     nelcon(2,*),nelem,ne,stiffness,i,j,c,d,ii,jj,nope,nring,
     &     indexe,imat,konl(255),ncyc,nasym,isym,iperturb(*)
      real*8 vold(0:mi(2),*)
      real*8 co(3,*),s(765,765),sm(765,765),ff(765),
     &     elcon(0:ncmat_,ntmat_,*),
     &     vh,xkv,sbar(3,255),tbar(3,255),tmax
!
      character*16 cval
      integer ilen
!
      nope=ichar(lakonl(8:8))
!
!     '_' is the Galerkin marker: fbares.py --symmetric writes it as the
!     ring-size character.  It decodes past any node count, so mastruct and
!     mafillsmas allocate the full support, and the element assembles
!     K_vol = K V_h sbar^T sbar below -- a Gram matrix, positive semidefinite
!     at every K.  The pairing is carried by the deck itself, with no
!     environment switch.  nasym still comes from ncyc below, so the element
!     rides the asymmetric storage path, which stores both triangles and
!     factors any general matrix; for a symmetric one the two are equal.
!     The default (a letter) is Onishi's Petrov-Galerkin tbar^T sbar.
!
      isym=0
      if(lakonl(3:3).eq.'_') isym=1
!
!     The edge ring size travels in lakon(3:3) as LETTERS[nring-1], so that
!     mastruct.c and mafillsmas.f can skip the identically-zero outer block of
!     K_vol = tbar^T sbar -- tbar is the UNSMOOTHED edge divergence and is
!     nonzero only on the ring.  The generator (fbares.py) writes the ring
!     first in the connectivity, so rows 1..nring are the ring.
!
      nring=ichar(lakonl(3:3))-ichar('A')+1
      if(nring.lt.1.or.nring.gt.nope) nring=nope
      if(isym.eq.1) nring=nope
!
      if(nope.gt.255) then
        write(*,*) '*ERROR in e_c3d_u3: edge element ',nelem
        write(*,*) '       spans ',nope,' nodes (',3*nope,' DOF).'
        write(*,*) '       The limit is 255 nodes (765 DOF) -- the most'
        write(*,*) '       that lakon(8:8) can encode.  The c=1 stencil'
        write(*,*) '       reaches 173 nodes and fits;'
        write(*,*) '       c=2 reaches 494 and cannot be an element.'
        call exit(201)
      endif
      indexe=ipkon(nelem)
      do i=1,nope
        konl(i)=kon(indexe+i)
      enddo
!
!     Zero the ACTUAL extent, and zero sm too -- see e_c3d_u2 / patch 0008.
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
!     number of cyclic smoothings, eq. (6)-(7).  The paper recommends 1 or 2
!     for nu <= 0.49 and states that the optimum varies with Poisson's ratio.
!
      ncyc=1
      call getenv('CCX_FBAR_C',cval)
      ilen=len_trim(cval)
      if(ilen.gt.0) read(cval(1:ilen),*) ncyc
      if((ncyc.lt.0).or.(ncyc.gt.3)) then
        write(*,*) '*ERROR in e_c3d_u3: CCX_FBAR_C =',ncyc
        write(*,*) '       must be 0..3'
        call exit(201)
      endif
!
      imat=ielmat(1,nelem)
      if(nelcon(1,imat).ne.2) then
        write(*,*) '*ERROR in e_c3d_u3: element',nelem,' needs an'
        write(*,*) '       isotropic *ELASTIC card (2 constants)'
        call exit(201)
      endif
!
!     NLGEOM: the finite-strain consistent tangent replaces the
!     small-strain operator.  Its rows are the ring and its columns the
!     support -- the same ring x support structure the deck encodes -- so
!     the assembly is unchanged; nasym stays 1.
!
      if(iperturb(1).gt.1) then
        call u3nltan(co,kon,ipkon,lakon,ne,konl,nope,vold,elcon,
     &       nelcon,ielmat,mi,ncmat_,ntmat_,imat,ncyc,s,nring,nelem)
        nasym=1
        return
      endif
!
      call u3vol(co,kon,ipkon,lakon,ne,konl,nope,vh,sbar,tbar,xkv,
     &     nelem,ielmat,elcon,nelcon,mi,ncmat_,ntmat_,imat,ncyc)
      if(vh.le.0.d0) return
!     THE RING ORDERING IS NOT TAKEN ON TRUST.
!
!     The row loop below stops at nring, so if the generator did not put the
!     edge ring first in the connectivity -- or if lakon(3:3) disagrees with
!     what it wrote -- real tbar rows would be dropped SILENTLY and the
!     volumetric operator would be quietly wrong, in a way no patch test can
!     see (a uniform field has nothing for the missing rows to carry).
!
!     u3vol builds tbar only from the tets that contain the edge, so it MUST
!     vanish past nring.  Checking that here costs one pass over the stencil
!     and turns the whole class of ordering bugs into a named element number.
!
      tmax=0.d0
      do i=1,nring
        do c=1,3
          tmax=max(tmax,dabs(tbar(c,i)))
        enddo
      enddo
      do i=nring+1,nope
        do c=1,3
          if(dabs(tbar(c,i)).gt.1.d-12*tmax) then
            write(*,*) '*ERROR in e_c3d_u3: element',nelem
            write(*,*) '       lakon(3:3) says the ring is the'
            write(*,*) '       first',nring,' nodes, but tbar is'
            write(*,*) '       nonzero at position',i,' node',konl(i)
            write(*,*) '       Deck and element disagree on the ring:'
            write(*,*) '       regenerate with a matching fbares.py'
            write(*,*) '       (ring first, nring in the label).'
            call exit(201)
          endif
        enddo
      enddo
!
!     K_vol = xkv * tbar^T sbar.  Row index carries the TEST space, column the
!     TRIAL space; swapping them is the transpose and is wrong.  The row loop
!     runs only over the ring (tbar vanishes past nring); the column loop spans
!     the whole support, which sbar carries.  In the symmetric pairing the
!     stress is paired with the smoothed divergence itself, sbar^T sbar, dense
!     over the support and positive semidefinite.
!
      if(isym.eq.0) then
        do i=1,nring
          do c=1,3
            ii=3*(i-1)+c
            do j=1,nope
              do d=1,3
                jj=3*(j-1)+d
                s(ii,jj)=xkv*tbar(c,i)*sbar(d,j)
              enddo
            enddo
          enddo
        enddo
      else
        do i=1,nope
          do c=1,3
            ii=3*(i-1)+c
            do j=1,nope
              do d=1,3
                jj=3*(j-1)+d
                s(ii,jj)=xkv*sbar(c,i)*sbar(d,j)
              enddo
            enddo
          enddo
        enddo
      endif
!
!     Tell the caller the global matrix is asymmetric.  At ncyc = 0 the chain
!     is the identity, sbar = tbar, and the element is symmetric -- that case
!     is plain selective ES-FEM-T4 and needs no asymmetric path.  The
!     symmetric pairing is symmetric by construction at every ncyc.
!
      if(ncyc.gt.0) nasym=1
!
      return
      end

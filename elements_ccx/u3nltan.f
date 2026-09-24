!
!     F-barES-FEM-T4: U3 CONSISTENT TANGENT, finite strain (NLGEOM).
!
!     The exact derivative of the u3nlforce internal force with respect to
!     the element's nodal displacements, through the SAME operator chain
!     the force uses:
!
!         du -> dF -> dF~ -> (dJ~, dJbar) -> dFbar -> dT -> dA -> df
!
!     with the ring weights (E row) at both ends, the cyclic chain (S row)
!     inside dJbar, and the Hencky tangent taken in the eigenbasis of
!     B = Fbar Fbar^T.  Both terms of the current-volume factor count:
!
!         df = sum_ring ww V_h [ dJ~ (g A) + J~ (g dA) ].
!
!     ROWS AND COLUMNS.  The force reaches the ring nodes only (the spread
!     carries the E row) but depends on every node of the chain through
!     Jbar.  The tangent therefore has RING rows and SUPPORT columns --
!     exactly the ring x support structure mastruct.c allocates from
!     lakon(3:3) -- and this routine fills rows 1..3*nring of s, so deck
!     and element matrix stay in agreement.
!
!     Per support column the chain is walked in reverse mode: every
!     derivative is accumulated where it is used, and the result lands in
!     the ring rows.  The Hencky map dT/dFbar is applied through the
!     eigenbasis, with the degenerate divided difference guarded to its
!     limit 1/(2 w_i).
!
!     The Python prototype this was derived from is `FbarNL.tangent` in
!     tests/smoothing_proto.py, verified against central differences and
!     against the small-strain operator (tests/verify_tangent.py).
!
      subroutine u3nltan(co,kon,ipkon,lakon,ne,konl,nope,v,elcon,nelcon,
     &     ielmat,mi,ncmat_,ntmat_,imatf,ncyc,s,nring,nelem)
!
      implicit none
!
      character*8 lakon(*)
      integer mi(*)
      integer kon(*),ipkon(*),ne,konl(*),nope,nring,nelem,ncyc,
     &     i,j,k,c,ie,ipos,n4(4),ielmat(mi(3),*),nelcon(2,*),
     &     ncmat_,ntmat_,imatf,na,nb,ihit1,ihit2,ic,nl,ll,je,
     &     imat,a,q,lb,apos,m,ir,icx,ic2,la,itch
      real*8 v(0:mi(2),*),elcon(0:ncmat_,ntmat_,*),s(765,765)
      real*8 co(3,*),xl(3,4),shp(4,4),xsj
      real*8 fe(3,3),ft(3,3),fi(3,3),fb(3,3),bb(3,3),hh(3,3),aa(3,3),
     &     tt(3,3),g(3,4),ev(3),vv(3,3),ffun(3,3),
     &     dft(3,3),dfb(3,3),dt(3,3),da(3,3),dbe(3,3),dhe(3,3),
     &     db(3,3),dhh(3,3),tmp(3,3),
     &     dw,c1,djt,detj,detb,trh,trd,third,
     &     ek,eg,e,un,vh,scr,scr2,alpha,beta,onethird
!
      integer maxel
      parameter(maxel=4096)
      integer elcur(maxel),elnxt(maxel),elring(maxel),ncur,nnxt,nring2
      real*8 wcur(maxel),wnxt(maxel),wring(maxel)
      integer iring(maxel)
      integer, allocatable :: loc4a(:,:),cstart(:),clist(:),
     &     rstart(:),rlist(:)
      real*8, allocatable :: gcha(:,:,:),finva(:,:,:),feca(:,:,:),
     &     deta(:)
!
      integer mapdone,maxn,nlen
      integer, allocatable, save :: nstart(:),nlist(:)
      real*8, allocatable, save :: volel(:)
      save mapdone,maxn
      data mapdone /0/
!
      onethird=1.d0/3.d0
!
!     Node -> U4 element map, built once (same as u3nlforce; saved locals
!     are per subroutine, so each routine keeps its own copy).
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
!     ---- eq. (8): ring of the edge ---------------------------------------
!
      na=konl(1)
      nb=konl(2)
      if((na.gt.maxn).or.(nb.gt.maxn)) then
        write(*,*) '*ERROR in u3nltan: element',nelem,' edge node'
        write(*,*) '       ',na,' or ',nb,' is in no U4 element'
        call exit(201)
      endif
      nring2=0
      vh=0.d0
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
        nring2=nring2+1
        if(nring2.gt.maxel) then
          write(*,*) '*ERROR in u3nltan: element',nelem,' ring'
          call exit(201)
        endif
        elring(nring2)=i
        wring(nring2)=volel(i)/6.d0
        vh=vh+volel(i)/6.d0
      enddo
      if(vh.le.0.d0) return
      do i=1,nring2
        wring(i)=wring(i)/vh
        elcur(i)=elring(i)
        wcur(i)=wring(i)
      enddo
      ncur=nring2
!
!     ---- eqs. (6)-(7): the same cyclic walk u3nlforce runs ----------------
!
      do ic=1,ncyc
        nnxt=0
        do i=1,ncur
          je=elcur(i)
          do j=1,4
            k=kon(ipkon(je)+j)
            scr=0.d0
            do ll=nstart(k)+1,nstart(k+1)
              if(ielmat(1,nlist(ll)).ne.imatf) cycle
              scr=scr+volel(nlist(ll))/4.d0
            enddo
            if(scr.le.0.d0) cycle
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
                  write(*,*) '*ERROR in u3nltan: element',nelem,
     &                 ' chain exceeds maxel'
                  call exit(201)
                endif
                elnxt(nnxt)=nl
                wnxt(nnxt)=0.d0
                ipos=nnxt
              endif
              wnxt(ipos)=wnxt(ipos)
     &             +wcur(i)*0.25d0*volel(nl)/4.d0/scr
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
!     ---- chain kinematics, cached once for the column loop ---------------
!
      allocate(gcha(3,4,max(ncur,1)))
      allocate(finva(3,3,max(ncur,1)))
      allocate(feca(3,3,max(ncur,1)))
      allocate(deta(max(ncur,1)))
      allocate(loc4a(4,max(ncur,1)))
      allocate(cstart(nope+1))
      allocate(clist(4*max(ncur,1)))
      allocate(rstart(nope+1))
      allocate(rlist(4*max(nring2,1)))
      do i=1,ncur
        je=elcur(i)
        call u3grad(co,kon,ipkon,je,g,n4)
        do a=1,4
          do j=1,3
            gcha(j,a,i)=g(j,a)
          enddo
        enddo
        do j=1,4
          loc4a(j,i)=0
          do k=1,nope
            if(konl(k).eq.n4(j)) then
              loc4a(j,i)=k
              exit
            endif
          enddo
          if(loc4a(j,i).eq.0) then
            write(*,*) '*ERROR in u3nltan: element',nelem,' node',n4(j)
            write(*,*) '       of the chain is not in its connectivity'
            call exit(201)
          endif
        enddo
        do c=1,3
          do j=1,3
            fe(c,j)=0.d0
          enddo
          fe(c,c)=1.d0
        enddo
        do c=1,3
          do j=1,4
            do k=1,3
              fe(k,c)=fe(k,c)+v(k,n4(j))*g(c,j)
            enddo
          enddo
        enddo
        call u3det3(fe,deta(i))
        if(deta(i).le.0.d0) then
!         A trial state can invert an element with the full Newton
!         step; ccx's own cutback handles that, so return a null block
!         for this trial instead of aborting the run (e_c3d flags the
!         same state with nmethod=0 and carries on).
          return
        endif
        call u3inv3(fe,finva(1,1,i))
        do c=1,3
          do j=1,3
            feca(c,j,i)=fe(c,j)
          enddo
        enddo
      enddo
!
!     for each local node, the chain entries containing it (CSR)
!
      do lb=1,nope+1
        cstart(lb)=0
      enddo
      do i=1,ncur
        do j=1,4
          cstart(loc4a(j,i)+1)=cstart(loc4a(j,i)+1)+1
        enddo
      enddo
      do lb=2,nope+1
        cstart(lb)=cstart(lb)+cstart(lb-1)
      enddo
      do i=1,ncur
        do j=1,4
          cstart(loc4a(j,i))=cstart(loc4a(j,i))+1
          clist(cstart(loc4a(j,i)))=i
        enddo
      enddo
      do lb=nope,1,-1
        cstart(lb+1)=cstart(lb)
      enddo
      cstart(1)=0
!
!     ring entries must be a subset of the chain; index them
!
      do ir=1,nring2
        itch=0
        do i=1,ncur
          if(elcur(i).eq.elring(ir)) then
            itch=i
            exit
          endif
        enddo
        if(itch.eq.0) then
          write(*,*) '*ERROR in u3nltan: element',nelem,
     &         ' ring element is not in its chain'
          call exit(201)
        endif
        iring(ir)=itch
      enddo
!
!     for each local node, the ring entries containing it (CSR)
!
      do lb=1,nope+1
        rstart(lb)=0
      enddo
      do ir=1,nring2
        do j=1,4
          la=loc4a(j,iring(ir))
          if(la.gt.nring) then
            write(*,*) '*ERROR in u3nltan: element',nelem,
     &           ' a ring node is beyond nring in the connectivity'
            call exit(201)
          endif
          rstart(la+1)=rstart(la+1)+1
        enddo
      enddo
      do lb=2,nope+1
        rstart(lb)=rstart(lb)+rstart(lb-1)
      enddo
      do ir=1,nring2
        do j=1,4
          la=loc4a(j,iring(ir))
          rstart(la)=rstart(la)+1
          rlist(rstart(la))=ir
        enddo
      enddo
      do lb=nope,1,-1
        rstart(lb+1)=rstart(lb)
      enddo
      rstart(1)=0
!
!     ---- ring kinematics: F~, J~, Fbar, Hencky T --------------------------
!
      imat=ielmat(1,elring(1))
      e=elcon(1,1,imat)
      un=elcon(2,1,imat)
      ek=e/(3.d0*(1.d0-2.d0*un))
      eg=e/(2.d0*(1.d0+un))
!
      do k=1,3
        do j=1,3
          ft(k,j)=0.d0
        enddo
      enddo
      do ir=1,nring2
        icx=iring(ir)
        do c=1,3
          do j=1,3
            ft(c,j)=ft(c,j)+wring(ir)*feca(c,j,icx)
          enddo
        enddo
      enddo
      call u3det3(ft,detj)
      if(detj.le.0.d0) return
      call u3inv3(ft,fi)
      detb=0.d0
      do i=1,ncur
        detb=detb+wcur(i)*deta(i)
      enddo
      if(detb.le.0.d0) return
      alpha=(detb/detj)**onethird
      beta=1.d0/(3.d0*detb)
      do c=1,3
        do j=1,3
          fb(c,j)=alpha*ft(c,j)
        enddo
      enddo
      do c=1,3
        do j=1,3
          bb(c,j)=0.d0
        enddo
      enddo
      do c=1,3
        do j=1,3
          do k=1,3
            bb(c,j)=bb(c,j)+fb(c,k)*fb(j,k)
          enddo
        enddo
      enddo
      call u3eig3(bb,ev,vv)
      trh=0.d0
      do j=1,3
        if(ev(j).le.0.d0) return
        trh=trh+0.5d0*dlog(ev(j))
      enddo
      do c=1,3
        do j=1,3
          hh(c,j)=0.d0
          do k=1,3
            hh(c,j)=hh(c,j)+vv(c,k)*0.5d0*dlog(ev(k))*vv(j,k)
          enddo
        enddo
      enddo
      do c=1,3
        do j=1,3
          scr=0.d0
          if(c.eq.j) scr=1.d0
          tt(c,j)=ek*trh*scr
     &         +2.d0*eg*(hh(c,j)-trh*onethird*scr)
        enddo
      enddo
      do c=1,3
        do j=1,3
          aa(c,j)=0.d0
          do k=1,3
            aa(c,j)=aa(c,j)+fi(c,k)*tt(k,j)
          enddo
        enddo
      enddo
!
!     Hencky divided differences (h_i - h_j)/(w_i - w_j), limit 1/(2w)
!
      do i=1,3
        do j=1,3
          if(i.eq.j) then
            ffun(i,j)=0.5d0/ev(i)
          else
            dw=ev(i)-ev(j)
            if(dabs(dw).gt.1.d-10*dabs(ev(i))) then
              ffun(i,j)=(0.5d0*dlog(ev(i))
     &             -0.5d0*dlog(ev(j)))/dw
            else
              ffun(i,j)=0.5d0/ev(i)
            endif
          endif
        enddo
      enddo
!
!     ---- the column loop: one support component at a time ----------------
!
      do lb=1,nope
        do q=1,3
          do c=1,3
            do j=1,3
              dft(c,j)=0.d0
            enddo
          enddo
!
!     dJbar: every chain element on this node contributes
!
          scr2=0.d0
          do ic2=cstart(lb)+1,cstart(lb+1)
            icx=clist(ic2)
            apos=1
            do j=1,4
              if(loc4a(j,icx).eq.lb) then
                apos=j
                exit
              endif
            enddo
            scr=0.d0
            do m=1,3
              scr=scr+finva(m,q,icx)*gcha(m,apos,icx)
            enddo
            scr2=scr2+wcur(icx)*deta(icx)*scr
          enddo
          detb=scr2
!
!     dF~: ring weights only
!
          do ic2=rstart(lb)+1,rstart(lb+1)
            ir=rlist(ic2)
            icx=iring(ir)
            apos=1
            do j=1,4
              if(loc4a(j,icx).eq.lb) then
                apos=j
                exit
              endif
            enddo
            do j=1,3
              dft(q,j)=dft(q,j)+wring(ir)*gcha(j,apos,icx)
            enddo
          enddo
!
!     dFbar, dT, dA
!
          c1=0.d0
          do m=1,3
            do j=1,3
              c1=c1+fi(m,j)*dft(j,m)
            enddo
          enddo
          djt=detj*c1
          do c=1,3
            do j=1,3
              dfb(c,j)=alpha*(dft(c,j)-onethird*c1*ft(c,j))
     &             +beta*fb(c,j)*detb
            enddo
          enddo
          do c=1,3
            do j=1,3
              db(c,j)=0.d0
              do k=1,3
                db(c,j)=db(c,j)+dfb(c,k)*fb(j,k)+fb(c,k)*dfb(j,k)
              enddo
            enddo
          enddo
          do c=1,3
            do j=1,3
              dbe(c,j)=0.d0
              do k=1,3
                do m=1,3
                  dbe(c,j)=dbe(c,j)+vv(k,c)*db(k,m)*vv(m,j)
                enddo
              enddo
            enddo
          enddo
          do c=1,3
            do j=1,3
              dhe(c,j)=ffun(c,j)*dbe(c,j)
            enddo
          enddo
          do c=1,3
            do j=1,3
              dhh(c,j)=0.d0
              do k=1,3
                do m=1,3
                  dhh(c,j)=dhh(c,j)+vv(c,k)*dhe(k,m)*vv(j,m)
                enddo
              enddo
            enddo
          enddo
          trd=dhh(1,1)+dhh(2,2)+dhh(3,3)
          do c=1,3
            do j=1,3
              scr=0.d0
              if(c.eq.j) scr=1.d0
              dt(c,j)=ek*trd*scr
     &             +2.d0*eg*(dhh(c,j)-trd*onethird*scr)
            enddo
          enddo
!     tmp = dF~ A must be COMPLETE before dA is formed: fusing the two
!     loops reads tmp(k,j) before row k has been computed whenever the
!     column's dF~ row is not row 1.
!
          do c=1,3
            do j=1,3
              tmp(c,j)=0.d0
              do k=1,3
                tmp(c,j)=tmp(c,j)+dft(c,k)*aa(k,j)
              enddo
            enddo
          enddo
          do c=1,3
            do j=1,3
              da(c,j)=0.d0
              do k=1,3
                da(c,j)=da(c,j)-fi(c,k)*tmp(k,j)
     &               +fi(c,k)*dt(k,j)
              enddo
            enddo
          enddo
!
!     spread over the ring: V_h [ dJ~ (g A) + J~ (g dA) ]
!
!     The spread reaches EVERY ring node: the dependence on the column
!     is carried by dJ~ and dA, not by which nodes the ring elements hold.
!
          do ir=1,nring2
            icx=iring(ir)
            do a=1,4
              la=loc4a(a,icx)
              do c=1,3
                scr=0.d0
                do k=1,3
                  scr=scr+gcha(k,a,icx)*aa(k,c)
                enddo
                scr2=0.d0
                do k=1,3
                  scr2=scr2+gcha(k,a,icx)*da(k,c)
                enddo
                s(3*(la-1)+c,3*(lb-1)+q)=s(3*(la-1)+c,3*(lb-1)+q)
     &               +wring(ir)*vh*(djt*scr+detj*scr2)
              enddo
            enddo
          enddo
        enddo
      enddo
!
      deallocate(gcha,finva,feca,deta,loc4a)
      deallocate(cstart,clist,rstart,rlist)
!
      return
      end

!
!     F-barES-FEM-T4: U3 FINITE-STRAIN internal force.
!
!     The Hencky-elastic force of eqs. (1)-(17) at a given displacement
!     field, for one U3 edge domain:
!
!         f_(a,c) = sum_h ww_he (J~ V_h) g_e[a,k] (F~^-1 T)_k,c
!
!     with the edge-smoothed F~ = E F (eq. 1), the cyclic chain
!     Jbar = S J(F) (eqs. 6-8), the recombination
!     Fbar = (Jbar/J~)^(1/3) F~ (eq. 11), and the Hencky stress of Fbar
!     (eqs. 12-13).  It is the exact finite-strain counterpart of the
!     small-strain path in u3vol.f: as u -> 0, T -> K tr(eps~) I + 2G dev,
!     F~ -> I and J~ -> 1, so its linearisation is the element matrix
!     e_c3d_u3 assembles.
!
!     The geometry and the operator chain are the SAME ones u3vol.f
!     builds: u3grad, the V_e/6 weighting of eq. (8), the V_n of eq. (6)
!     restricted to one material, and the cyclic walk of eq. (7).  The
!     weight list is copied before the cycles (the ring, = E row) and then
!     cycled in place (the chain, = S row), so the force uses exactly the
!     operator the stiffness uses.
!
!     This file is called from resultsmech_u3.f when iperturb(1) = 1
!     (NLGEOM).  No stress is written, as everywhere in an F-bar deck.
!
      subroutine u3nlforce(co,kon,ipkon,lakon,ne,konl,nope,v,elcon,
     &     nelcon,ielmat,mi,ncmat_,ntmat_,imatf,ncyc,fn,nelem)
!
      implicit none
!
      character*8 lakon(*)
      integer mi(*)
      integer kon(*),ipkon(*),ne,konl(*),nope,nelem,ncyc,
     &     i,j,k,c,ie,ipos,n4(4),ielmat(mi(3),*),nelcon(2,*),
     &     ncmat_,ntmat_,imatf,na,nb,ihit1,ihit2,ic,nl,ll,je,
     &     imat,a
      real*8 v(0:mi(2),*),elcon(0:ncmat_,ntmat_,*),fn(0:mi(2),*)
      real*8 co(3,*),xl(3,4),shp(4,4),xsj
      real*8 fe(3,3),ft(3,3),fb(3,3),bb(3,3),hh(3,3),aa(3,3),
     &     tt(3,3),fi(3,3),g(3,4),ev(3),vv(3,3),
     &     vh,scr,detj,detb,detf,trh,ek,eg,e,un,third
!
      integer maxel
      parameter(maxel=4096)
      integer elcur(maxel),elnxt(maxel),elring(maxel),ncur,nnxt,nring
      real*8 wcur(maxel),wnxt(maxel),wring(maxel)
!
      integer mapdone,maxn,nlen
      integer, allocatable, save :: nstart(:),nlist(:)
      real*8, allocatable, save :: volel(:)
      save mapdone,maxn
      data mapdone /0/
!
!     Node -> U4 element map, built once.  Kept separate from the copy in
!     u3vol.f on purpose: saved locals are per subroutine, and sharing
!     them would need a module or a common block in a tree that uses
!     neither.  A few MB on the production mesh.
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
!     ---- eq. (8): the elements at edge h, weights V_e/6 / V_h ------------
!
      na=konl(1)
      nb=konl(2)
      if((na.gt.maxn).or.(nb.gt.maxn)) then
        write(*,*) '*ERROR in u3nlforce: element',nelem,' edge node'
        write(*,*) '       ',na,' or ',nb,' is in no U4 element'
        call exit(201)
      endif
      nring=0
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
        nring=nring+1
        if(nring.gt.maxel) then
          write(*,*) '*ERROR in u3nlforce: element',nelem,' ring'
          write(*,*) '       exceeds maxel =',maxel
          call exit(201)
        endif
        elring(nring)=i
        wring(nring)=volel(i)/6.d0
        vh=vh+volel(i)/6.d0
      enddo
      if(vh.le.0.d0) return
      do i=1,nring
        wring(i)=wring(i)/vh
        elcur(i)=elring(i)
        wcur(i)=wring(i)
      enddo
      ncur=nring
!
!     ---- material of the edge domain and ring kinematics: F~ = E F ------
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
      do i=1,nring
        je=elring(i)
        call u3grad(co,kon,ipkon,je,g,n4)
        do k=1,3
          do j=1,3
            fe(k,j)=0.d0
          enddo
          fe(k,k)=1.d0
        enddo
        do c=1,3
          do j=1,4
            do k=1,3
              fe(k,c)=fe(k,c)+v(k,n4(j))*g(c,j)
            enddo
          enddo
        enddo
        do c=1,3
          do k=1,3
            ft(k,c)=ft(k,c)+wring(i)*fe(k,c)
          enddo
        enddo
      enddo
      call u3det3(ft,detj)
      if(detj.le.0.d0) return
!
!     ---- eqs. (6)-(7): c cycles of A = P Q, on the weight list -----------
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
                  write(*,*) '*ERROR in u3nlforce: element',nelem,
     &                 ' chain exceeds maxel =',maxel
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
!     ---- Jbar = S J(F) ---------------------------------------------------
!
      detb=0.d0
      do i=1,ncur
        je=elcur(i)
        call u3grad(co,kon,ipkon,je,g,n4)
        do k=1,3
          do j=1,3
            fe(k,j)=0.d0
          enddo
          fe(k,k)=1.d0
        enddo
        do c=1,3
          do j=1,4
            do k=1,3
              fe(k,c)=fe(k,c)+v(k,n4(j))*g(c,j)
            enddo
          enddo
        enddo
        call u3det3(fe,detf)
        detb=detb+wcur(i)*detf
      enddo
      if(detb.le.0.d0) return
!
!     ---- eq. (11): Fbar = (Jbar/J~)^(1/3) F~ ------------------------------
!
      scr=(detb/detj)**(1.d0/3.d0)
      do c=1,3
        do k=1,3
          fb(k,c)=scr*ft(k,c)
        enddo
      enddo
!
!     ---- eqs. (12)-(13): Hencky stress of Fbar ----------------------------
!
      do c=1,3
        do k=1,3
          bb(k,c)=0.d0
        enddo
      enddo
      do c=1,3
        do k=1,3
          do j=1,3
            bb(k,c)=bb(k,c)+fb(k,j)*fb(c,j)
          enddo
        enddo
      enddo
      call u3eig3(bb,ev,vv)
      trh=0.d0
      do j=1,3
        if(ev(j).le.0.d0) return
        trh=trh+0.5d0*dlog(ev(j))
      enddo
      third=trh/3.d0
      do c=1,3
        do k=1,3
          hh(k,c)=0.d0
          do j=1,3
            hh(k,c)=hh(k,c)+vv(k,j)*0.5d0*dlog(ev(j))*vv(c,j)
          enddo
        enddo
      enddo
      do c=1,3
        do k=1,3
          scr=0.d0
          if(k.eq.c) scr=1.d0
          tt(k,c)=ek*trh*scr+2.d0*eg*(hh(k,c)-third*scr)
        enddo
      enddo
!
!     ---- eq. (17) push-forward and spread ---------------------------------
!
      call u3inv3(ft,fi)
      do c=1,3
        do k=1,3
          aa(k,c)=0.d0
          do j=1,3
            aa(k,c)=aa(k,c)+fi(k,j)*tt(j,c)
          enddo
        enddo
      enddo
      do i=1,nring
        je=elring(i)
        call u3grad(co,kon,ipkon,je,g,n4)
        do a=1,4
          do c=1,3
            scr=0.d0
            do k=1,3
              scr=scr+g(k,a)*aa(k,c)
            enddo
            fn(c,n4(a))=fn(c,n4(a))
     &           +wring(i)*detj*vh*scr
          enddo
        enddo
      enddo
!
      return
      end
!
!
      subroutine u3det3(a,d)
!
!     Determinant of a 3x3 matrix by the rule of Sarrus.
!
      implicit none
      real*8 a(3,3),d
      d=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
     &     -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
     &     +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      return
      end
!
!
      subroutine u3inv3(a,ai)
!
!     Inverse of a 3x3 matrix via the adjugate; assumes det > 0.
!
      implicit none
      real*8 a(3,3),ai(3,3),d
      call u3det3(a,d)
      ai(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))/d
      ai(1,2)=(a(1,3)*a(3,2)-a(1,2)*a(3,3))/d
      ai(1,3)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))/d
      ai(2,1)=(a(2,3)*a(3,1)-a(2,1)*a(3,3))/d
      ai(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))/d
      ai(2,3)=(a(1,3)*a(2,1)-a(1,1)*a(2,3))/d
      ai(3,1)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))/d
      ai(3,2)=(a(1,2)*a(3,1)-a(1,1)*a(3,2))/d
      ai(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))/d
      return
      end
!
!
      subroutine u3eig3(a,w,v)
!
!     Cyclic Jacobi eigendecomposition of a symmetric 3x3 matrix:
!     a = v diag(w) v^T with the eigenvectors in the COLUMNS of v.
!     Thirty sweeps is far past convergence for 3x3; the loop exits on
!     the first sweep whose off-diagonal weight is negligible.
!
      implicit none
      integer i,j,k,ip,iq,isw
      real*8 a(3,3),w(3),v(3,3),b(3,3),off,tau,t,c,s,apq,
     &     bpk,bqk,bkp,bkq,vkp,vkq
      do i=1,3
        do j=1,3
          b(i,j)=a(i,j)
          v(i,j)=0.d0
        enddo
        v(i,i)=1.d0
      enddo
      do isw=1,30
        off=dabs(b(1,2))+dabs(b(1,3))+dabs(b(2,3))
        if(off.le.1.d-15*(dabs(b(1,1))+dabs(b(2,2))+dabs(b(3,3))
     &       +1.d-30)) exit
        do ip=1,2
          do iq=ip+1,3
            apq=b(ip,iq)
            if(apq.eq.0.d0) cycle
            tau=(b(iq,iq)-b(ip,ip))/(2.d0*apq)
            if(tau.ge.0.d0) then
              t=1.d0/(tau+dsqrt(1.d0+tau*tau))
            else
              t=-1.d0/(-tau+dsqrt(1.d0+tau*tau))
            endif
            c=1.d0/dsqrt(1.d0+t*t)
            s=t*c
            do k=1,3
              bkp=b(k,ip)
              bkq=b(k,iq)
              b(k,ip)=c*bkp-s*bkq
              b(k,iq)=s*bkp+c*bkq
            enddo
            do k=1,3
              bpk=b(ip,k)
              bqk=b(iq,k)
              b(ip,k)=c*bpk-s*bqk
              b(iq,k)=s*bpk+c*bqk
            enddo
            do k=1,3
              vkp=v(k,ip)
              vkq=v(k,iq)
              v(k,ip)=c*vkp-s*vkq
              v(k,iq)=s*vkp+c*vkq
            enddo
          enddo
        enddo
      enddo
      do i=1,3
        w(i)=b(i,i)
      enddo
      return
      end

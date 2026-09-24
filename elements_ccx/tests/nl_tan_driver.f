!
!     Test driver for u3nltan: reads a mesh, a displacement field and a
!     list of U3 edge elements (nope, nring, konl), calls the tangent for
!     each and writes the global (row, column, value) triplets of the
!     assembled matrix.  Used by tests/verify_nltan.py, which compares the
!     result with FbarNL.tangent.
!
      program tldrv
      implicit none
      integer nn,ne,nedge,i,j,k,e,ipkon(200000),kon(800000),n4(4)
      integer ielmat(3,200000),mi(10),ncmat_,ntmat_,nelcon(2,10)
      real*8 co(3,200000),elcon(0:10,10,10),v(0:3,200000)
      real*8 s(765,765)
      character*8 lakon(200000)
      integer konl(255),nope,nring,imatf,ncyc,nelem,gr,gc,jj,ll
      read(5,*) nn,ne
      do i=1,nn
        read(5,*) co(1,i),co(2,i),co(3,i)
      enddo
      do e=1,ne
        read(5,*) (kon(4*(e-1)+j),j=1,4)
        ipkon(e)=4*(e-1)
        lakon(e)='U4'
        do j=1,3
          ielmat(j,e)=1
        enddo
      enddo
      read(5,*) ncmat_,ntmat_
      read(5,*) nelcon(1,1),nelcon(2,1)
      read(5,*) elcon(1,1,1),elcon(2,1,1)
      do i=1,nn
        read(5,*) (v(k,i),k=1,3)
      enddo
      mi(1)=1
      mi(2)=3
      imatf=1
      ncyc=1
      read(5,*) nedge
      do k=1,nedge
        read(5,*) nope,nring
        read(5,*) (konl(i),i=1,nope)
        do jj=1,765
          do ll=1,765
            s(jj,ll)=0.d0
          enddo
        enddo
        nelem=k
        call u3nltan(co,kon,ipkon,lakon,ne,konl,nope,v,elcon,nelcon,
     &       ielmat,mi,ncmat_,ntmat_,imatf,ncyc,s,nring,nelem)
        do jj=1,3*nring
          gr=3*(konl((jj-1)/3+1)-1)+mod(jj-1,3)+1
          do ll=1,3*nope
            if(s(jj,ll).ne.0.d0) then
              gc=3*(konl((ll-1)/3+1)-1)+mod(ll-1,3)+1
              write(6,'(i8,1x,i8,1x,e24.16)') gr,gc,s(jj,ll)
            endif
          enddo
        enddo
      enddo
      end

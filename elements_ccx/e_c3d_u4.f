!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2025 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine e_c3d_u4(s,sm,ff,stiffness,mass,buckling,rhsi,coriolis)
!
!     F-barES-FEM-T4: U4 -- the null base tetrahedron of the F-barES-FEM-T4 element.
!
!     The base tets of an F-bar deck are retyped to U4 so that the U2
!     (deviatoric) and U3 (volumetric) smoothing domains can read their
!     geometry through the user-element node->element map.  U4 itself
!     carries no stiffness, so the 4-node (12-DOF) element matrices are
!     zeroed and nothing else is done.
!
      implicit none
!
      integer stiffness,mass(*),buckling,rhsi,coriolis,i,j
!
      real*8 s(765,765),sm(765,765),ff(765)
!
      if(stiffness.eq.1) then
         do i=1,12
            do j=1,12
               s(i,j)=0.d0
            enddo
         enddo
      endif
!
      if((mass(1).eq.1).or.(buckling.eq.1)) then
         do i=1,12
            do j=1,12
               sm(i,j)=0.d0
            enddo
         enddo
      endif
!
      if(rhsi.eq.1) then
         do i=1,12
            ff(i)=0.d0
         enddo
      endif
!
      return
      end

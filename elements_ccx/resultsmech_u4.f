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
      subroutine resultsmech_u4(fn,calcul_fn,calcul_qa,mi,nelem)
!
!     F-barES-FEM-T4: U4 -- the null base tetrahedron of the F-barES-FEM-T4 element.
!
!     U4 carries no stiffness (e_c3d_u4 zeroes its 12-DOF block) and so
!     contributes no internal force.  fn is already zeroed by resultsini
!     before the stress pass, so there is nothing to do here: the previous
!     version wrote fn(i,nelem)=0 with nelem the ELEMENT index into fn's
!     NODE-indexed second dimension, which overflows fn whenever the element
!     count exceeds the node count and segfaults in the stress pass.
!
      implicit none
!
      integer calcul_fn,calcul_qa,mi(*),nelem
!
      real*8 fn(0:mi(2),*)
!
      return
      end

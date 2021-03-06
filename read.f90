SUBROUTINE read_dim(numele,numnods,nedge,nelem_ghost, &
              nnode_ghost,ngauss,nnodBC,ndofBC,ndofOP,nvdw, &
              ngauss_vdw,ng_tot,nneigh,ninrange,indentation)
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)

read(90,*) 
read(90,*) 
read(90,*) 
read(90,*) numele
read(90,*) 
read(90,*) numnods
read(90,*) 
read(90,*) nedge
read(90,*) 
read(90,*) nelem_ghost
read(90,*) 
read(90,*) nnode_ghost
read(90,*) 
read(90,*) ngauss
read(90,*) 
read(90,*) nnodBC
read(90,*) 
read(90,*) ndofBC
read(90,*) 
read(90,*) ndofOP
read(90,*) 
read(90,*) nvdw
if (nvdw.eq.1) then
  read(90,*) 
  read(90,*) ngauss_vdw
  read(90,*) 
  read(90,*) ng_tot
  read(90,*)
  read(90,*) nneigh
  read(90,*) 
  read(90,*) ninrange
endif
END SUBROUTINE read_dim

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE read_general(ylength,mat1,nW_hat,crit,imperfect,fact_imp)
USE data_mat
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(material) :: mat1
DIMENSION :: crit(2)

read(91,*) 
read(91,*) 
read(91,*) 
read(91,'(d30.17)') ylength
read(91,*) 
read(91,'(d30.17)') mat1%A0
read(91,*) 
read(91,*) mat1%nCode_Pot
read(91,*) 
if (mat1%nCode_Pot.eq.1) then
  read(91,'(d30.17)') mat1%Vs(1)
  read(91,'(d30.17)') mat1%Vs(2)
  read(91,'(d30.17)') mat1%Va(1)
  read(91,'(d30.17)') mat1%Va(2)
else if (mat1%nCode_Pot.eq.2) then
  read(91,'(d30.17)')  mat1%A1
  read(91,'(d30.17)')  mat1%Vs(1)
  read(91,'(d30.17)')  mat1%Vs(2)
  read(91,'(d30.17)')  mat1%Vs(3)
  read(91,'(d30.17)')  mat1%Va(1)
  read(91,'(d30.17)')  mat1%Va(2)
  read(91,'(d30.17)')  mat1%Va(3)
else if (mat1%nCode_Pot.eq.22) then
  continue
else if (mat1%nCode_Pot.eq.3) then
  read(91,'(d30.17)') mat1%Vs(1)
  read(91,'(d30.17)') mat1%Va(1)
endif
read(91,*) 
read(91,'(2d30.17)') mat1%E(1,:)
read(91,'(2d30.17)') mat1%E(2,:)
read(91,'(2d30.17)') mat1%E(3,:)
read(91,*) 
read(91,'(d30.17)') mat1%s0
read(91,*) 
read(91,*) nW_hat
read(91,*) 
read(91,'(d30.17)') crit(1)
read(91,'(d30.17)') crit(2)
read(91,*) 
read(91,*) imperfect
read(91,*) 
read(91,'(d30.17)') fact_imp

END SUBROUTINE read_general
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE read_mesh(mesh0)
USE data_mesh
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(mesh) :: mesh0

read(94,*)
read(94,*)
read(94,*)
do i=1,mesh0%numele
  read(94,*)
  read(94,*) ielem, mesh0%connect(ielem)%vertices(1:3)
  read(94,*) mesh0%connect(ielem)%num_neigh_elem
  read(94,*) mesh0%connect(ielem)%num_neigh_vert
  do jj=1,12
    read(94,*) mesh0%connect(ielem)%neigh_elem(jj), mesh0%connect(ielem)%neigh_vert(jj)
  enddo
  read(94,*) mesh0%connect(ielem)%code_bc(1:3)
enddo
!read(94,*)
!do inod=1,mesh0%numnods
!  read(94,'(13i6)') mesh0%ntable(inod,:)
!enddo
read(94,*)
do inod=1,mesh0%nedge
  read(94,'(3i6)') mesh0%nghost_tab(inod,:)
enddo
read(94,*) 
do inod=1,mesh0%nelem_ghost
  read(94,'(3i6)') mesh0%elem_ghost(inod)
enddo
read(94,*)
do inod=1,mesh0%nnode_ghost
  read(94,'(3i6)') mesh0%node_ghost(inod)
enddo

END SUBROUTINE read_mesh
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE read_zero(numele,F0,J0)
USE data_tensor22
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(tensor22) :: F0(numele)
REAL(8) :: J0(numele)

read(95,*) 
read(95,*) 
do ielem=1,numele
  read(95,'(d30.17)') J0(ielem)
  read(95,'(2d30.17)') F0(ielem)%val(1,1),F0(ielem)%val(1,2)
  read(95,'(2d30.17)') F0(ielem)%val(2,1),F0(ielem)%val(2,2)
enddo


END SUBROUTINE read_zero
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE read_BCs(BCs)
USE data_BC
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(BC_data):: BCs
PARAMETER (pi0=3.141592653589793238d0)

read(92,*) 
read(92,*) 
read(92,*) 
read(92,*) BCs%nloadstep
read(92,*)
read(92,*) BCs%nCodeLoad
read(92,*) 
do i=1,BCs%ndofBC
  read(92,*) BCs%mdofBC(i)
enddo
read(92,*)
do i=1,BCs%ndofOP
  read(92,*) BCs%mdofOP(i)
enddo
read(92,*) 
do i=1,BCs%nnodBC
  read(92,*) BCs%mnodBC(i,1:2)
enddo
read(92,*) 
read(92,'(3d30.17)') BCs%rotation(1,1:3)
read(92,'(3d30.17)') BCs%rotation(2,1:3)
read(92,'(3d30.17)') BCs%rotation(3,1:3)
read(92,*) 
read(92,'(3d30.17)') BCs%xc(1:3)
read(92,*) 
read(92,'(d30.17)') BCs%value

!!$  BCs%value=5.d0*pi0/180.d0
!!$  BCs%rotation(1,1:3)=[ dcos(BCs%value/BCs%nloadstep),dsin(BCs%value/BCs%nloadstep),0.d0]
!!$  BCs%rotation(2,1:3)=[-dsin(BCs%value/BCs%nloadstep),dcos(BCs%value/BCs%nloadstep),0.d0]
!!$  BCs%rotation(3,1:3)=[0.d0,0.d0,1.d0]
!!$  BCs%xc(1:3)=[6.25d0,0.d0,0.d0]


END SUBROUTINE read_BCs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE read_vdw(vdw1,nneigh)
USE data_vdw
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(vdw_data):: vdw1


if (vdw1%nvdw.eq.1) then
   vdw1%neval=0
   open(unit=93,file='nano_vdw.dat',status='unknown')
   read(93,*) 
   read(93,*) vdw1%meval
   read(93,*) 
   read(93,'(d30.17)') vdw1%r_cut
   read(93,*)
   read(93,'(d30.17)') vdw1%r_bond
   read(93,*) 
   read(93,'(d30.17)') vdw1%sig
   read(93,*)
   read(93,'(d30.17)') vdw1%a
   read(93,*) 
   read(93,'(d30.17)') vdw1%y0
   read(93,*) 
   read(93,'(2d30.17)') vdw1%Vcut(1:2)
   read(93,*) 
   do i=1,12
      read(93,'(10d30.17)') vdw1%shapef(1:vdw1%ngauss_vdw,i)
   enddo
   read(93,*) 
   read(93,'(10d30.17)') vdw1%weight(1:vdw1%ngauss_vdw)
   read(93,*) 
   if (nneigh.gt.0) then
      do i=1,vdw1%ng_tot
         read(93,*) vdw1%near(i,0:nneigh)
      enddo
   else
      vdw1%ntub_pos(0)=0
      do i=1,vdw1%ntubes
         read(56,*) vdw1%ntub_pos(i)
      enddo
   endif
   close(93)
endif

!!$!**************************************************************************************
!!$!**************************************************************************************
!!$
!!$open(unit=32,file='bondadd.dat',status='unknown')
!!$read(32,*) 
!!$read(32,*) vdw1%ngaussfind
!!$read(32,*) vdw1%bonddensity
!!$
!!$ALLOCATE(vdw1%indexelem(vdw1%ngaussfind),vdw1%shapefff(12,vdw1%ngaussfind) ,STAT=istat)
!!$ALLOCATE(vdw1%shapeff(2,12) , vdw1%weightt(1:2),STAT=istat)
!!$ALLOCATE(vdw1%x1(vdw1%ngaussfind,3),vdw1%x2(vdw1%ngaussfind,3) ,STAT=istat)
!!$
!!$
!!$do i=1,vdw1%ngaussfind
!!$   read(32,*) vdw1%indexelem(i)
!!$enddo
!!$
!!$
!!$do i=1,vdw1%ngaussfind
!!$   read(32,'(10d30.17)') vdw1%shapefff(1:12,i)
!!$enddo
!!$
!!$do i=1,12
!!$   read(32,'(10d30.17)') vdw1%shapeff(1:2,i)
!!$enddo
!!$
!!$read(32,'(10d30.17)') vdw1%weightt(1:2)
!!$
!!$close(32)



!**************************************************************************************
!*************************************************************************************

END SUBROUTINE read_vdw
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE read_indent(vdw1)
USE data_vdw
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(vdw_data):: vdw1

if (vdw1%indentation.eq.1) then
  open(unit=98,file='nano_indent.dat',status='unknown')
  read(98,*)
  read(98,'(d30.17)') vdw1%denter_a
  read(98,*)
  read(98,'(d30.17)') vdw1%denter_r
  read(98,*)
  read(98,*) vdw1%mdent_eval
  read(98,*)
  read(98,*)vdw1%denter_x
  close(98)
endif

write(*,*) vdw1%denter_a, vdw1%denter_r, vdw1%mdent_eval
write(*,*) vdw1%denter_x

END SUBROUTINE read_indent

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
SUBROUTINE read_config(numnods,numele,ngauss,x0,eta)
USE data_vector2
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
REAL(8) :: x0(3*(numnods))
TYPE(vector2) :: eta(ngauss,(numele))

read(96,*)
read(96,*) 
read(96,*) 
do i=1,numnods
  read(96,'(3d30.17)') x0(3*i-2:3*i)
enddo
read(96,*) 
do ielem=1,numele
  do igauss=1,ngauss
    read(96,'(2d30.17)') eta(igauss,ielem)%val(1:2)
  enddo
enddo

END SUBROUTINE read_config
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE write_config(numnods,numele,ngauss,x0,eta)
USE data_vector2
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
REAL(8) :: x0(3*(numnods))
TYPE(vector2) :: eta(ngauss,(numele))

write(97,*) 'Config Data'
write(97,*) '-----------'
write(97,*) 'Nodal positions'
do i=1,numnods
  write(97,'(3d30.17)') x0(3*i-2:3*i)
enddo
write(97,*) 'Inner displacements'
do ielem=1,numele
  do igauss=1,ngauss
    write(97,'(2d30.17)') eta(igauss,ielem)%val(1:2)
  enddo
enddo

END SUBROUTINE write_config
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


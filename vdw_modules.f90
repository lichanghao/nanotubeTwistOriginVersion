
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE Vvdw(a,A_p,sig,y0,Vs)
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
REAL(8), INTENT(IN) :: a,A_p,sig,y0
REAL(8), INTENT(OUT) :: Vs(2)

a1=sig/a
a6=(a1)**6
a12=a6*a6
y06=y0**6


Vs(1)=0.5d0*y06*a12 - a6
Vs(2)=a1*6.d0/sig*(-y06*a12 + a6)

Vs=Vs*A_p/(sig**6)

END SUBROUTINE Vvdw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE Vvdw_cut(a,A_p,sig,y0,r_cut,Vaux,Vs)
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
REAL(8), INTENT(IN) :: a,A_p,sig,y0,r_cut,Vaux(2)
REAL(8), INTENT(OUT) :: Vs(2) 

call Vvdw(a,A_p,sig,y0,Vs)

if (a.lt.r_cut) then
  Vs(1)=Vs(1)-Vaux(2)*(a-r_cut)-Vaux(1)
  Vs(2)=Vs(2)-Vaux(2)
else
  Vs=0.d0
endif
!write(*,*)"111111111",a,vs(1),vs(2)
END SUBROUTINE Vvdw_cut

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!****************************************************************************************************
!****************************************************************************************************
!****************************************************************************************************
!****************************************************************************************************


SUBROUTINE Vvdw_cutbond(a,Vs)
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
REAL(8), INTENT(IN) :: a
REAL(8), INTENT(OUT) :: Vs(2) 
REAL(8):: AA
!call Vvdw(a,A_p,sig,y0,Vs)

!!$if (a.lt.r_cut) then
!!$  Vs(1)=Vs(1)-Vaux(2)*(a-r_cut)-Vaux(1)
!!$  Vs(2)=Vs(2)-Vaux(2)
!!$else
!!$  Vs=0.d0
!!$endif
AA=25.00
if(a.ge.0.34.and.a.lt.0.64)then
Vs(1)=AA/2.0*(a-0.34)**2
Vs(2)=AA*(a-0.34)
else
Vs(1)=0.0
Vs(2)=0.0
endif
!write(*,*)"22222222",a
!if(a.gt.0.4)then
!vs(1)=0.0
!vs(2)=0.0
!endif

END SUBROUTINE Vvdw_cutbond



!****************************************************************************************************
!****************************************************************************************************
!****************************************************************************************************
!****************************************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE xlocate_gauss(vdw1,mesh0,x0)
  USE data_mesh
  USE data_vdw
  implicit REAL(8) (a-h,o-z)
  implicit INTEGER*4 (i-n)
!!$include 'mpif.h' !MPI
  TYPE(mesh) :: mesh0
  TYPE(vdw_data):: vdw1
  REAL(8) :: x0(3*(mesh0%numnods+mesh0%nedge))
  DIMENSION :: xneigh(12,3), neigh(12), neigh1(12)

  ! Compute the position of all the gauss points
  vdw1%x=0.d0

  !*********************************************************************************************
  !*********************************************************************************************
  vdw1%x1=0.d0
  vdw1%x2=0.d0
  !*********************************************************************************************
  !*********************************************************************************************

  do ielem=1,mesh0%numele
     neigh=mesh0%connect(ielem)%neigh_vert
     xneigh(:,1)=x0(3*neigh(:)-2)
     xneigh(:,2)=x0(3*neigh(:)-1)
     xneigh(:,3)=x0(3*neigh(:))
     do ig=1,vdw1%ngauss_vdw
        igg=(ielem-1)*vdw1%ngauss_vdw+ig
        do inod=1,12
           do idof=1,3
              vdw1%x(igg,idof)= vdw1%x(igg,idof) + vdw1%shapef(ig,inod)*xneigh(inod,idof)
           enddo
        enddo
     enddo

          !        write(*,*)"11111111"
     !*******************************************************************************************************
     !*******************************************************************************************************
!!$     if(mod(ielem,2).eq.1)then
!!$        neigh=mesh0%connect(ielem)%neigh_vert
!!$     else
!!$        neigh1=mesh0%connect(ielem)%neigh_vert
!!$        neigh(1)=  neigh1(6) 
!!$        neigh(2)=  neigh1(10)
!!$        neigh(3)=  neigh1(3)
!!$        neigh(4)=  neigh1(7)
!!$        neigh(5)=  neigh1(11)
!!$        neigh(6)=  neigh1(1)
!!$        neigh(7)=  neigh1(4)
!!$        neigh(8)=  neigh1(8)
!!$        neigh(9)=  neigh1(12)
!!$        neigh(10)= neigh1(2)
!!$        neigh(11)= neigh1(5)
!!$        neigh(12)= neigh1(9)
!!$     endif
!!$     xneigh(:,1)=x0(3*neigh(:)-2)
!!$     xneigh(:,2)=x0(3*neigh(:)-1)
!!$     xneigh(:,3)=x0(3*neigh(:))
!!$
!!$
!!$     do ig=1,2
!!$        igg=(ielem-1)*2+ig
!!$        !    write(*,*)"11111111111111111111111111111111"
!!$        if(igg.gt.vdw1%ngaussfind) goto 100
!!$
!!$        do inod=1,12
!!$           do idof=1,3
!!$
!!$              vdw1%x1(igg,idof)= vdw1%x1(igg,idof) + vdw1%shapeff(ig,inod)*xneigh(inod,idof)
!!$           enddo
!!$        enddo
!!$     enddo
!!$ !                          write(*,*) ielem, vdw1%x1(1,1), vdw1%x1(2,1),
!!$     do ig=1,2
!!$        igg=(ielem-1)*2+ig
!!$
!!$        if(igg.gt.vdw1%ngaussfind) goto 100
!!$
!!$        iielem= vdw1%indexelem(igg)
!!$!         write(*,*)igg,iielem,vdw1%ntub_pos(vdw1%ntubes-1),vdw1%ntub_pos(vdw1%ntubes)
!!$        !              write(*,*)iielem
!!$        !        neigh=mesh0%connect(iielem)%neigh_vert
!!$        !          write(*,*)"2222",iielem 
!!$
!!$        if(mod(iielem,2).eq.1)then
!!$           neigh=mesh0%connect(iielem)%neigh_vert
!!$        else
!!$           neigh1=mesh0%connect(iielem)%neigh_vert
!!$           neigh(1)=  neigh1(6) 
!!$           neigh(2)=  neigh1(10)
!!$           neigh(3)=  neigh1(3)
!!$           neigh(4)=  neigh1(7)
!!$           neigh(5)=  neigh1(11)
!!$           neigh(6)=  neigh1(1)
!!$           neigh(7)=  neigh1(4)
!!$           neigh(8)=  neigh1(8)
!!$           neigh(9)=  neigh1(12)
!!$           neigh(10)= neigh1(2)
!!$           neigh(11)= neigh1(5)
!!$           neigh(12)= neigh1(9)
!!$        endif
!!$
!!$        xneigh(:,1)=x0(3*neigh(:)-2)
!!$        xneigh(:,2)=x0(3*neigh(:)-1)
!!$        xneigh(:,3)=x0(3*neigh(:))
!!$
!!$        do inod=1,12
!!$           do idof=1,3
!!$              vdw1%x2(igg,idof)= vdw1%x2(igg,idof) + vdw1%shapefff(inod,igg)*xneigh(inod,idof)
!!$           enddo
!!$        enddo
        !          write(*,*)  vdw1%x2(igg,3)
!      write(*,*) ielem, vdw1%x1(1,1), vdw1%x1(2,1)

!!$     enddo
     !                     write(*,*)  vdw1%x2(1,1), vdw1%x2(2,1)
     !*******************************************************************************************************
     !*******************************************************************************************************
     !          write(*,*)ielem
100 continue
  enddo
!!$  write(*,*)  vdw1%x2(1,1), vdw1%x2(1,2), vdw1%x2(1,3)
!!$  write(*,*)  vdw1%x2(2,1), vdw1%x2(2,2), vdw1%x2(2,3)
!!$  write(*,*)  vdw1%x2(3,1), vdw1%x2(3,2), vdw1%x2(3,3)
!!$  write(*,*)  vdw1%x2(4,1), vdw1%x2(4,2), vdw1%x2(4,3)
!!$    stop   



END SUBROUTINE xlocate_gauss
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE find_in_range(vdw1,mesh0,x0,J0_val)
USE data_mesh
USE data_vdw
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
include 'mpif.h' !MPI
TYPE(mesh) :: mesh0
TYPE(vdw_data):: vdw1
REAL(8) :: x0(3*(mesh0%numnods+mesh0%nedge))
REAL(8) :: vmax(3), vmix(3), dx(3), vec(3), J0_val
INTEGER(4) ::  nloc(3)
LOGICAL :: ask1, ask2, ask3


call MPI_COMM_RANK( MPI_COMM_WORLD, id, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

! This is previous to binning
vmax(1)=maxval(x0(1:3*mesh0%numnods:3))
vmix(1)=minval(x0(1:3*mesh0%numnods:3))
vmax(2)=maxval(x0(2:3*mesh0%numnods:3))
vmix(2)=minval(x0(2:3*mesh0%numnods:3))
vmax(3)=maxval(x0(3:3*mesh0%numnods:3))
vmix(3)=minval(x0(3:3*mesh0%numnods:3))
rc=vdw1%r_cut
fact=1.05d0

vdw1%nx(:)=max(int((vmax(:)-vmix(:))/(fact*rc)),1)
dx(:)=(vmax(:)-vmix(:))/vdw1%nx(:)

! Check whether we need to enlarge binning variables
if ((vdw1%nx(1).gt.vdw1%nxdim(1)).or. &
    (vdw1%nx(2).gt.vdw1%nxdim(2)).or.(vdw1%nx(3).gt.vdw1%nxdim(3))) then
   vdw1%nxdim(:)=(vdw1%nx(:)+1)
   write(*,*) ' Need to enlarge the bins '
   DEALLOCATE(vdw1%binx)
   ALLOCATE(vdw1%binx(3,0:maxval(vdw1%nxdim(:))),STAT=istat)
   if (istat/=0) STOP '**** Not enough memory **** vdw1'
   DEALLOCATE(vdw1%bin)
   ng_per_bin=int(sqrt(2.)*rc*rc*vdw1%ngauss_vdw/(J0_val/2.)*1.3 *3)
   ALLOCATE(vdw1%bin(vdw1%nxdim(1),vdw1%nxdim(2),vdw1%nxdim(3), &
            0:ng_per_bin),STAT=istat)
   if (istat/=0) STOP '**** Not enough memory **** vdw2'
endif
vdw1%binx=0.
vdw1%bin=0

do idim=1,3
  vdw1%binx(idim,0:vdw1%nx(idim))=vmix(idim)+[(ijk,ijk=0,vdw1%nx(idim))]*dx(idim)
  vdw1%binx(idim,0)=vdw1%binx(idim,0)-1000.*rc
  vdw1%binx(idim,vdw1%nx(idim))=vdw1%binx(idim,vdw1%nx(idim))+1000.*rc
enddo

! This is previous to binning

! bin all the gauss points
do ielem=1,mesh0%numele
  do ig=1,vdw1%ngauss_vdw
    igg=(ielem-1)*vdw1%ngauss_vdw+ig
!   this is binning
    do idof=1,3
      do inx=1,vdw1%nx(idof)
        if ((vdw1%x(igg,idof).gt.vdw1%binx(idof,inx-1)).and. &
            (vdw1%x(igg,idof).le.vdw1%binx(idof,inx))) then
          nloc(idof)=inx
          exit
        endif
      enddo
    enddo
    vdw1%bin(nloc(1),nloc(2),nloc(3),0)=vdw1%bin(nloc(1),nloc(2),nloc(3),0)+1 
!   the 0th component is the number of gp's in this bin
    vdw1%bin(nloc(1),nloc(2),nloc(3),vdw1%bin(nloc(1),nloc(2),nloc(3),0))=igg
    vdw1%gbin(igg,:)=nloc(:)
!   this is binning
  enddo
enddo


! This loop is parallelized_MPI
! find list of neighbors
vdw1%neigh =0 ! it is set to zero
do ielem=1+id,mesh0%numele,nprocs
do ig=1,vdw1%ngauss_vdw
  i=(ielem-1)*vdw1%ngauss_vdw+ig
  nloc=vdw1%gbin(i,:)
  !##This is special case##
  if (vdw1%nneigh.lt.0) then
    ! find to which NT i belongs
    my_nt=0
    do i_nt=1,vdw1%ntubes
      if ((i.gt.vdw1%ntub_pos(i_nt-1)).and.(i.le.vdw1%ntub_pos(i_nt))) then
         my_nt=i_nt
        exit
      endif
    enddo
    if (my_nt.eq.0) STOP ' Guat is my nanotube?'
  endif
  !##This is special case##
  do inx=max(1,nloc(1)-1), min(vdw1%nx(1),nloc(1)+1)
    do iny=max(1,nloc(2)-1), min(vdw1%nx(2),nloc(2)+1)
      do inz=max(1,nloc(3)-1), min(vdw1%nx(3),nloc(3)+1)
         do inbin=1,vdw1%bin(inx,iny,inz,0)
           j=vdw1%bin(inx,iny,inz,inbin)
           if (vdw1%nneigh.gt.0) then
             ask1=(.not.(any( vdw1%near(i,1:vdw1%near(i,0)).eq.j )))
           else
           !##This is special case##
             ! find to which NT j belongs
             myour_nt=0
             do i_nt=1,vdw1%ntubes
               if ((j.gt.vdw1%ntub_pos(i_nt-1)).and.(j.le.vdw1%ntub_pos(i_nt))) then
                 myour_nt=i_nt
                 exit
               endif
             enddo
             if (myour_nt.eq.0) STOP ' Guat is my nanotube?'
!             if (myour_nt.eq.vdw1%ntubes) then ! last shell, look for self-contact
!               jelem=int(j/vdw1%ngauss_vdw)+1
!               ask1=(.not.(any(mesh0%connect(ielem)%neigh_elem(:).eq.jelem))).and.(ielem.ne.jelem)
!             else
!               ask1=(myour_nt.gt.my_nt) ! only look in immediate outer shell
 !              ask1=(myour_nt.eq.my_nt+1) ! only look in immediate outer shell

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
           ask1=(myour_nt.eq.my_nt+1)  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            if((myour_nt.eq.1).and.(my_nt.eq.1))then
               if(i.gt.82*2*vdw1%ngauss_vdw*3)then
                  if((j-i).ge.(82*2*vdw1%ngauss_vdw*3))then
                     !  write(*,*)j-i,202*2*vdw1%ngauss_vdw*3
                     ask1=.true.
                     !        write(*,*)j,i
                  endif
               else
                  if((j-i).ge.(82*2*vdw1%ngauss_vdw*3).and.(j.lt.(82*2*vdw1%ngauss_vdw*7)))then
                     !  write(*,*)j-i,202*2*vdw1%ngauss_vdw*3
                     ask1=.true.
                     !        write(*,*)j,i
                  endif
               endif
            endif
           
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 
!             endif
           !##This is special case##
           endif
           ask2=j.gt.i
           if (ask1.and.ask2) then ! make sure i and j are not bonded
             vec(1)=vdw1%x(i,1)-vdw1%x(j,1)
             vec(2)=vdw1%x(i,2)-vdw1%x(j,2)
             vec(3)=vdw1%x(i,3)-vdw1%x(j,3)
             dist=sqrt( (vec(1))**2 + (vec(2))**2 + (vec(3))**2)
             ask3=(dist.lt.vdw1%r_cut)
             if (ask3) then ! if within range...
                if (vdw1%neigh(i,0)+1.ge.vdw1%ninrange-3) STOP ' Insufficient numax of inrange'
                vdw1%neigh(i,0)=vdw1%neigh(i,0)+1
                vdw1%neigh(i,vdw1%neigh(i,0))=j
             endif
           endif
         enddo
      enddo
    enddo
  enddo
enddo
enddo

END SUBROUTINE find_in_range
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SUBROUTINE vdw_elem(s0,mesh0,ielem,neigh_i,vdw1,J0,W_vdw,f_vdw)
  USE data_vdw
  USE data_mesh
  implicit REAL(8) (a-h,o-z)
  implicit INTEGER*4 (i-n)
  TYPE(vdw_data):: vdw1
  TYPE(mesh) :: mesh0
  REAL(8) :: Vs(2), vec(3),J0(mesh0%numele), f_vdw(3*(mesh0%numnods+mesh0%nedge))
  REAL(8) :: f_i(12,3), f_ii(12,3), f_j(12,3)
  REAL(8) :: fact(3),W_vdwbond
  INTEGER(4) :: neigh_i(12), neigh_j(12), nloc(3),neigh1(12), neigh_ii(12)
  !*************************************************
   REAL(8) :: ss0,density
  density=0.000
  ss0=density*2.0/s0*1.143! *1.143 due to the out most layer
  !*************************************************

  W_vdw=0.d0
  do ig=1,vdw1%ngauss_vdw
     i=(ielem-1)*vdw1%ngauss_vdw+ig
     f_ii=0.d0
     ! this is neighbor inplementation
     do jk=1,vdw1%neigh(i,0)
        j=vdw1%neigh(i,jk)
        vec(1)=vdw1%x(i,1)-vdw1%x(j,1)
        vec(2)=vdw1%x(i,2)-vdw1%x(j,2)
        vec(3)=vdw1%x(i,3)-vdw1%x(j,3)
        dist=sqrt( (vec(1))**2 + (vec(2))**2 + (vec(3))**2)
        jelem=int((j-1)/vdw1%ngauss_vdw)+1
        jg=j-(jelem-1)*vdw1%ngauss_vdw
        neigh_j=mesh0%connect(jelem)%neigh_vert
        ! Here is the energy
        wei=vdw1%weight(ig)*vdw1%weight(jg)*J0(ielem)*J0(jelem) ! THIS WEIGHT TAKES CARE
        ! OF THE FACTOR 4 IN (57) JMPS
        call Vvdw_cut(dist,vdw1%a,vdw1%sig,vdw1%y0,vdw1%r_cut,vdw1%Vcut,Vs)
        W_vdw=W_vdw + Vs(1)*wei
        ! Here are the forces
        fact=Vs(2)/dist*vec*wei/s0/s0
        do inod=1,12
           do jdof=1,3
              f_i(inod,jdof)=  fact(jdof)*vdw1%shapef(ig,inod)
              f_j(inod,jdof)= -fact(jdof)*vdw1%shapef(jg,inod)
           enddo
        enddo
        ! Assemble these forces
        f_ii=f_ii+f_i
        do iq=1,12
           f_vdw(3*neigh_j(iq)-2:3*neigh_j(iq))=f_vdw(3*neigh_j(iq)-2:3*neigh_j(iq))+f_j(iq,:)
        enddo
     enddo
     do iq=1,12
        f_vdw(3*neigh_i(iq)-2:3*neigh_i(iq))=f_vdw(3*neigh_i(iq)-2:3*neigh_i(iq))+f_ii(iq,:)
     enddo
  enddo
  W_vdw=W_vdw/s0/s0


!!$  !****************************************************************************************************************
!!$  !****************************************************************************************************************
!!$
!!$  !write(*,*)'hereeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee'
!!$  !do iii=1,2
!!$  W_vdwbond=0.d0
!!$
!!$  do ig=1,2
!!$     i=(ielem-1)*2+ig
!!$
!!$     if(i.gt.vdw1%ngaussfind) exit
!!$
!!$     ! this is neighbor inplementation
!!$     !      if(ig.eq.1.and.mod(iii,2).eq.1)then
!!$     !         vec(1)=vdw1%x1(i,1)-vdw1%x2(i+1,1)
!!$
!!$     !      elseif(ig.eq.2.and.mod(iii,2).eq.1)then
!!$     !         vec(1)=vdw1%x1(i,1)-vdw1%x2(i-1,1)
!!$     !      else
!!$     vec(1)=vdw1%x1(i,1)-vdw1%x2(i,1)
!!$     !      endif
!!$     vec(2)=vdw1%x1(i,2)-vdw1%x2(i,2)
!!$     !        if(ig.eq.1.and.mod(iii,2).eq.1)then
!!$     !         vec(3)=vdw1%x1(i,3)-vdw1%x2(i+1,3)
!!$
!!$     !      elseif(ig.eq.2.and.mod(iii,2).eq.1)then
!!$     !         vec(3)=vdw1%x1(i,3)-vdw1%x2(i-1,3)
!!$     !      else
!!$     vec(3)=vdw1%x1(i,3)-vdw1%x2(i,3)
!!$     !      endif
!!$
!!$
!!$
!!$
!!$
!!$     !    vec(3)=vdw1%x1(i,3)-vdw1%x2(i,3)
!!$     dist=sqrt( (vec(1))**2 + (vec(2))**2 + (vec(3))**2)
!!$
!!$     !         dista=sqrt(vdw1%x1(i,2)**2+vdw1%x1(i,3)**2)
!!$     !         distb=sqrt(vdw1%x2(i,2)**2+vdw1%x2(i,3)**2)
!!$
!!$
!!$     !        write(*,*)"distanceeeeeeeeeeeeeeeeeeee",dista,distb
!!$     !         write(*,*)vdw1%x1(i,1),vdw1%x1(i,2),vdw1%x1(i,3)
!!$     !         write(*,*)vdw1%x2(i,1),vdw1%x2(i,2),vdw1%x2(i,3)
!!$     !   enddo
!!$     !enddo
!!$     !stop
!!$     !        write(*,*)vdw1%x1(i,1),vdw1%x2(i,1)
!!$     jelem=vdw1%indexelem(i)
!!$     !    jg=j-(jelem-1)*vdw1%ngauss_vdw
!!$
!!$     if(mod(jelem,2).eq.1)then
!!$        neigh_j=mesh0%connect(jelem)%neigh_vert
!!$     else
!!$        neigh1=mesh0%connect(jelem)%neigh_vert
!!$        neigh_j(1)=  neigh1(6) 
!!$        neigh_j(2)=  neigh1(10)
!!$        neigh_j(3)=  neigh1(3)
!!$        neigh_j(4)=  neigh1(7)
!!$        neigh_j(5)=  neigh1(11)
!!$        neigh_j(6)=  neigh1(1)
!!$        neigh_j(7)=  neigh1(4)
!!$        neigh_j(8)=  neigh1(8)
!!$        neigh_j(9)=  neigh1(12)
!!$        neigh_j(10)= neigh1(2)
!!$        neigh_j(11)= neigh1(5)
!!$        neigh_j(12)= neigh1(9)
!!$     endif
!!$
!!$
!!$     !   neigh_j=mesh0%connect(jelem)%neigh_vert
!!$     ! Here is the energy
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$     wei=vdw1%weightt(ig)*J0(ielem)/2 !divided by 2 due to the taking care of 4 of weight funcition
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$     call Vvdw_cutbond(dist,Vs)
!!$     W_vdwbond=W_vdwbond + Vs(1)*wei
!!$     ! Here are the forces
!!$     fact=Vs(2)/dist*vec*wei*ss0
!!$     do inod=1,12
!!$        do jdof=1,3
!!$           f_i(inod,jdof)=  fact(jdof)*vdw1%shapeff(ig,inod)
!!$           f_j(inod,jdof)= -fact(jdof)*vdw1%shapefff(inod,i)
!!$        enddo
!!$     enddo
!!$     ! Assemble these forces
!!$
!!$     do iq=1,12
!!$        f_vdw(3*neigh_j(iq)-2:3*neigh_j(iq))=f_vdw(3*neigh_j(iq)-2:3*neigh_j(iq))+f_j(iq,:)
!!$     enddo
!!$
!!$
!!$     if(mod(ielem,2).eq.1)then
!!$        neigh_ii=mesh0%connect(ielem)%neigh_vert
!!$     else
!!$        neigh1=mesh0%connect(ielem)%neigh_vert
!!$        neigh_ii(1)=  neigh1(6) 
!!$        neigh_ii(2)=  neigh1(10)
!!$        neigh_ii(3)=  neigh1(3)
!!$        neigh_ii(4)=  neigh1(7)
!!$        neigh_ii(5)=  neigh1(11)
!!$        neigh_ii(6)=  neigh1(1)
!!$        neigh_ii(7)=  neigh1(4)
!!$        neigh_ii(8)=  neigh1(8)
!!$        neigh_ii(9)=  neigh1(12)
!!$        neigh_ii(10)= neigh1(2)
!!$        neigh_ii(11)= neigh1(5)
!!$        neigh_ii(12)= neigh1(9)
!!$     endif
!!$
!!$
!!$     do iq=1,12
!!$        f_vdw(3*neigh_ii(iq)-2:3*neigh_ii(iq))=f_vdw(3*neigh_ii(iq)-2:3*neigh_ii(iq))+f_i(iq,:)
!!$     enddo
!!$  enddo
!!$
!!$  W_vdwbond=W_vdwbond*ss0            
!!$!     write(*,*)"wwwwwwwwwwwwwwwwwwwwwww",W_vdwbond
!!$  W_vdw=W_vdw+W_vdwbond
  !stop
  !****************************************************************************************************************
  !****************************************************************************************************************
  !if(ielem.gt.10) stop
END SUBROUTINE vdw_elem


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

SUBROUTINE get_reac(mesh0,forces,nCodeLoad,ylength, &
	        mdofBC,ndofBC,reaction1,reaction2,x0,BCs,iload)
USE data_mesh
USE data_BC
implicit REAL(8) (a-h,o-z)
implicit INTEGER*4 (i-n)
TYPE(mesh) :: mesh0
TYPE(BC_data):: BCs
DIMENSION forces(3*(mesh0%numnods)),mdofBC(ndofBC)
DIMENSION x0(3*(mesh0%numnods))
integer :: iload
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
DIMENSION nntube(BCs%nnodBC)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
reaction1=0.d0
reaction2=0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!
nntube=0
!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=40,file='forcefx.dat',status='unknown',position='append')! out put force informations of all tubes
open(unit=41,file='forcefy.dat',status='unknown',position='append')! out put force informations of all tubes
open(unit=42,file='forcefz.dat',status='unknown',position='append')! out put force informations of all tubes
open(unit=43,file='forcetwist.dat',status='unknown',position='append')! out put force informations of all tubes
open(unit=44,file='forceMz.dat',status='unknown',position='append')! out put force informations of all tubes
open(unit=45,file='forceMy.dat',status='unknown',position='append')! out put force informations of all tubes
!!!!!!!!!!!!!!!!!!!!!!


if (nCodeLoad.eq.3) then
!!$   do i=1,ndofBC/6
!!$     reaction1=reaction1 + forces(mdofBC(3*i))*(x0(mdofBC(3*i-1))-BCs%xc(2)) &
!!$  -forces(mdofBC(3*i-1))*(x0(mdofBC(3*i))-BCs%xc(3))
!!$   !  reaction2=reaction2 + forces(mdofBC(ndofBC/2+3*i-2))
!!$     reaction2=reaction2 +  forces(mdofBC(ndofBC/2+3*i))*(x0(mdofBC(ndofBC/2+3*i-1)) &
!!$-BCs%xc(2))-forces(mdofBC(ndofBC/2+3*i-1))*(x0(mdofBC(ndofBC/2+3*i))-BCs%xc(3))
!!$   end do


   do i=1,BCs%nnodBC
      if (BCs%mnodBC(i,2).eq.1) then
         reaction1=reaction1 + forces(mdofBC(3*i-2))
      else
         reaction2=reaction2 +  forces(mdofBC(3*i-2))
      end if
   enddo


   do i=1,BCs%nnodBC
      dis=sqrt((x0(mdofBC(3*i))**2+x0(mdofBC(3*i-1))**2))
      do ii=1,50  ! assume the largest number of tube is 50
         if((dis.gt.((ii-1)*0.34+0.17)).and.(dis.lt.((ii)*0.34+0.17)))then
            nntube(i)=ii
            goto 1003
         endif
      enddo
1003  continue
      if(nntube(i).eq.0)then
         write(*,*)"finding BC point - tube number wrong"
         stop
      endif
   enddo

   nmax=maxval(nntube(:))
      reactiontot1=0  !@@@@@@@@@@@@@@@  total force @@@@@@@@@@@@@@@@@@@@
      reactiontot2=0 !fy
      reactiontot3=0!fz
      reactiontot4=0!torque
      reactiontot5=0!Mz
      reactiontot6=0!My
      reactiontot11=0  !   fx  !!!!!!other edges
      reactiontot22=0 !fy
      reactiontot33=0!fz
      reactiontot44=0!torque
      reactiontot55=0!Mz
      reactiontot66=0!My

   do iii=1,nmax

      reaction1=0  !   fx
      reaction2=0 !fy
      reaction3=0!fz
      reaction4=0!torque
      reaction5=0!Mz
      reaction6=0!My
      reaction11=0  !   fx  !!!!!!other edges
      reaction22=0 !fy
      reaction33=0!fz
      reaction44=0!torque
      reaction55=0!Mz
      reaction66=0!My

      do iiii=1,BCs%nnodBC
         if(nntube(iiii).eq.iii)then
            if (BCs%mnodBC(iiii,2).eq.1) then
               reaction1=reaction1 + forces(mdofBC(3*iiii-2))
               reaction2=reaction2 + forces(mdofBC(3*iiii-1))
               reaction3=reaction3 + forces(mdofBC(3*iiii))
               reaction4=reaction4 + forces(mdofBC(3*iiii))*(x0(mdofBC(3*iiii-1))-BCs%xc(2)) &
                    -forces(mdofBC(3*iiii-1))*(x0(mdofBC(3*iiii))-BCs%xc(3))
               reaction5=reaction5+forces(mdofBC(3*iiii-2))*(x0(mdofBC(3*iiii-1))-BCs%xc(2))
               reaction6=reaction6+forces(mdofBC(3*iiii-2))*(x0(mdofBC(3*iiii))-BCs%xc(3))

               reactiontot1=reactiontot1 + forces(mdofBC(3*iiii-2))
               reactiontot2=reactiontot2 + forces(mdofBC(3*iiii-1))
               reactiontot3=reactiontot3 + forces(mdofBC(3*iiii))
               reactiontot4=reactiontot4 + forces(mdofBC(3*iiii))*(x0(mdofBC(3*iiii-1))-BCs%xc(2)) &
                    -forces(mdofBC(3*iiii-1))*(x0(mdofBC(3*iiii))-BCs%xc(3))
               reactiontot5=reactiontot5+forces(mdofBC(3*iiii-2))*(x0(mdofBC(3*iiii-1))-BCs%xc(2))
               reactiontot6=reactiontot6+forces(mdofBC(3*iiii-2))*(x0(mdofBC(3*iiii))-BCs%xc(3))

            endif
            if (BCs%mnodBC(iiii,2).eq.2) then
               reaction11=reaction11 + forces(mdofBC(3*iiii-2))
               reaction22=reaction22 + forces(mdofBC(3*iiii-1))
               reaction33=reaction33 + forces(mdofBC(3*iiii))
               reaction44=reaction44 + forces(mdofBC(3*iiii))*(x0(mdofBC(3*iiii-1))-BCs%xc(2)) &
                    -forces(mdofBC(3*iiii-1))*(x0(mdofBC(3*iiii))-BCs%xc(3))
               reaction55=reaction55+forces(mdofBC(3*iiii-2))*(x0(mdofBC(3*iiii-1))-BCs%xc(2))
               reaction66=reaction66+forces(mdofBC(3*iiii-2))*(x0(mdofBC(3*iiii))-BCs%xc(3))

               reactiontot11=reactiontot11 + forces(mdofBC(3*iiii-2))
               reactiontot22=reactiontot22 + forces(mdofBC(3*iiii-1))
               reactiontot33=reactiontot33 + forces(mdofBC(3*iiii))
               reactiontot44=reactiontot44 + forces(mdofBC(3*iiii))*(x0(mdofBC(3*iiii-1))-BCs%xc(2)) &
                    -forces(mdofBC(3*iiii-1))*(x0(mdofBC(3*iiii))-BCs%xc(3))
               reactiontot55=reactiontot55+forces(mdofBC(3*iiii-2))*(x0(mdofBC(3*iiii-1))-BCs%xc(2))
               reactiontot66=reactiontot66+forces(mdofBC(3*iiii-2))*(x0(mdofBC(3*iiii))-BCs%xc(3))
            endif

         endif
      enddo
      
!      write(40,'(I3,3f12.6)')iii,BCs%value*iload/BCs%nloadstep,reaction1,reaction11
!      write(41,'(I3,3f12.6)')iii,BCs%value*iload/BCs%nloadstep,reaction2,reaction22
!      write(42,'(I3,3f12.6)')iii,BCs%value*iload/BCs%nloadstep,reaction3,reaction33
!      write(43,'(I3,3f12.6)')iii,BCs%value*iload/BCs%nloadstep,reaction4,reaction44
!      write(44,'(I3,3f12.6)')iii,BCs%value*iload/BCs%nloadstep,reaction5,reaction55
!      write(45,'(I3,3f12.6)')iii,BCs%value*iload/BCs%nloadstep,reaction6,reaction66

   enddo
      write(40,'(a5,3f12.6)')"total",BCs%value*iload/BCs%nloadstep,reactiontot1,reactiontot11
      write(41,'(a5,3f12.6)')"total",BCs%value*iload/BCs%nloadstep,reactiontot2,reactiontot22
      write(42,'(a5,3f12.6)')"total",BCs%value*iload/BCs%nloadstep,reactiontot3,reactiontot33
      write(43,'(a5,3f12.6)')"total",BCs%value*iload/BCs%nloadstep,reactiontot4,reactiontot44
      write(44,'(a5,3f12.6)')"total",BCs%value*iload/BCs%nloadstep,reactiontot5,reactiontot55
      write(45,'(a5,3f12.6)')"total",BCs%value*iload/BCs%nloadstep,reactiontot6,reactiontot66

 close(40)
 close(41)
 close(42)
 close(43)
 close(44)
 close(45)

else if (nCodeLoad.eq.4) then
   do i=1,(ndofBC-3)/2
      reaction1=reaction1 + forces(mdofBC(i))
      reaction2=reaction2 + forces(mdofBC((ndofBC-3)/2+i))
   end do

else if (nCodeLoad.eq.20) then
   write(*,*)'we are writing the reaction forces'

else
   STOP 'Code of loading not implemented' 
end if

!reaction1=reaction1/ylength
!reaction2=reaction2/ylength

END SUBROUTINE get_reac

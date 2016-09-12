Module  param
implicit none
save
integer, parameter :: imax=150 !matrix dimension
integer :: im
real(8) :: c  !CFL number
end module param
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program wave

!variables

use param
implicit none

integer, parameter :: nmax=150
integer :: method,nm,i1,i2,i3,count,i, order, n, icase
real(8), dimension(imax) :: x,u,u1,u2,u3,ui
real(8), dimension(imax,nmax) :: ut
real(8) :: pi,int,dlength,delx,delt,a,tottime,w1,w2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!input data
tottime=0.150
pi=dacos(-1.0)
dlength=70.0
delx=1.0
delt=0.001250
a=200.0
int=0.0250

	print*, 'please inter a number' 
	print*, 'First upwind differencing: 1'
	print*, 'Lax - Wendroff	: 2'
	print*, 'Eulers BTCS : 3'
	read*, method
	
if ( method<1 .or. method>3) then 
		print*, "wrong number entered!"
		print*, " re-run the program."
		stop
end if


im=idint(dlength/delx)+1
nm=idint(tottime/delt)+1
c= A*delt/delx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!initial and boundary condition
w1= 40.0
w2=120.0
i1=5
i2=15
i3=25

do i=1,im

x(i)=dble(i-1)*delx

	if (i>i1 .and. i<i2+1) then

		u(i)=2*(x(i)-5)
	
	else if (i>i2 .and. i<i3+1) then
		
		u(i)= -2*(x(i)-15) +20
		
	else
		u(i)=0.0
		
	end if
	
ui(i)=u(i)

end do

count=0

do n=2,nm

	if(dabs (delt*dble(n-1)-int) <=1.d-7) then 
		order=1
	elseif(dabs (delt*dble(n-1)-int*2.0) <=1.d-7) then
		order=1
	elseif(dabs (delt*dble(n-1)-int*3.0) <=1.d-7) then
		order=1
	elseif(dabs (delt*dble(n-1)-int*4.0) <=1.d-7) then
		order=1
	elseif(dabs (delt*dble(n-1)-int*5.0) <=1.d-7) then
		order=1
	elseif(dabs (delt*dble(n-1)-int*6.0) <=1.d-7) then
		order=1
	else
		order=0
	end if
	
	
	select case (method)
		case(1)
	call fud(u)
		case(2)
	call lwm(u)
		case(3)
	call eubtcs(u)
	end select 
	
	
	if (order==1) then
		count=count+1
		
		do i=1,im
		
		if (count==1) then
		u1(i)=u(i)
		elseif (count==2) then
		u2(i)=u(i)
		elseif (count==3) then
		u3(i)=u(i)
		end if
		
	
		end do
		
	end if
	
	do i=1,im
	ut(i,n)= u(i)
	end do
	
end do

select case (method)
case(1)
open( unit=9, file='fupwind.dat',status='replace', action='write')
open( unit=19, file='tpfupwind.dat',status='replace', action='write')
case(2)
open( unit=9, file='l_wm.dat',status='replace', action='write')
open( unit=19, file='tpl_wm.dat',status='replace', action='write')
case(3)
open( unit=9, file='eu_btcs.dat',status='replace', action='write')
open( unit=19, file='tpeu_btcs.dat',status='replace', action='write')

end select
!printing results
write(9,*)
write(9,10)
write(9,*)
write(19,*) 'title="time marching" '
write(19,*) 'variable= "x", "t", "u" '
write(19,*) 'zone f=point, I=', im, 'j=', nm

do i=1,im
	write(9,20) x(i), ui(i), u1(i), u2(i), u3(i)
	
end do

do n=1, nm
		do i=1,im
		
			write(19,30) x(i), delt*dble(n-1), ut(i,n)
		end do
end do

10 format (1x,'x', 3x, 't=0.0', 3x, 't=0.025', 3x, 't=0.05', 3x, &
		't=0.075', 3x, 't=0.1', 3x, 't=0.125', 3x, 't=0.15')
			
20 format (1x, f7.3, 3x, f9.4, 6(3x, f9.4))

30 format (1x, d15.8, 3x, d15.8, 3x, d15.8)

close(9)
close(19)

end program wave
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!First upwinding method
subroutine fud(u)
use param
implicit none
integer :: i
real(8), dimension(imax) :: u,uold

do i=1,im
	uold(i)=u(i)
end do

do i=2,im-1
	u(i)= uold(i) - c*(uold(i)-uold(i-1))
end do

end subroutine fud
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!LAX-WENDROFF method
subroutine lwm(u)

use param
implicit none
integer :: i
real(8), dimension(imax) :: u, uold

do  i=1,im
	uold(i)=u(i)
end do

do i=2,im-1
	u(i)=uold(i)- c/2.0*(uold(i+1)- uold(i-1)) &
		+c*c/2.0*(uold(i+1)-2.0*uold(i)+ uold(i-1))
end do
return
end subroutine lwm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Euler's BTCS method
subroutine eubtcs(u)

use param
implicit none
integer :: i
real(8), dimension(imax) :: u, aa, bb, cc, dd

do i=2,im-1

	aa(i)= 0.50*c
	bb(i)= -1.0
	cc(i)= -0.50*c
	dd(i)= -u(i)
end do

call trid(aa,bb,cc,dd,u)
return
end subroutine eubtcs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Tridiogonal matrix inverse for implicit method
subroutine trid(aa,bb,cc,dd,u)
use param
implicit none
integer :: i
real(8), dimension(imax) :: h,g,u,aa,bb,cc,dd

h(1)=0.0
g(1)=u(1)

do i=2,im-1
	h(i)= cc(i)/ (bb(i)-aa(i)*h(i-1))
	g(i)= (dd(i)-aa(i)*g(i-1))/ (bb(i)-aa(i)*h(i-1))
end do

do i=im-1,2,-1
	u(i)=-h(i)*u(i+1)+g(i)
end do
return
end subroutine trid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		
			
	
		
		
	


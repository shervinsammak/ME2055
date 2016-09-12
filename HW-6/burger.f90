module prop2
implicit none
save 

integer, parameter :: imax=150
integer :: im
real(8) :: c

end module prop2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program burger

use prop2
implicit none

integer, parameter :: nmax=100
integer :: opt, nm, i1, i2, count, i, order, n,incase
integer :: dopt, modi, lim, ot

real(8), dimension(imax) :: x, u, u1, u2, u3, u4,u5,u6,u7,ui,e, a
real(8), dimension(imax,imax) :: ut
real(8) :: int,dlength, delx, delt,tottime, w1, w2, eps
integer(double precision) :: t1, t2, cr            ! timing variables

! use only two threads
  !$ call omp_set_num_threads(4)

  ! time the entire program
  call system_clock( t1, cr )
  
  
tottime=2.40
dlength=40.
delx=1.0
delt=0.2
int=0.4


print*, ' The LAX method: 1'
print*, ' The LAX_Wendroff method: 2'
print*, ' The MacCormack method: 3'
print*, ' The Beam and Warming method: 4'
read*, opt


dopt=0
eps=0.
modi=0
lim=0
ot=0
if (opt ==2) then
print*, '2nd order damping term: 1, y'
print*, '					   : 2, N'
read*, dopt

end if
if (dopt ==1) then 
print*, 'Epsilon :?'
read*, eps
end if

im=idint(dlength/delx) +1
nm=idint(tottime/delt) +1
c= delt/delx
w1=20.
i1= idint(w1/delx) +1


!$omp parallel default(none) &
  !$omp shared(dx,eps,L,Tn) private(i,x)
  !$omp do
  

do i=1,im
x(i)=dble(i-1)*delx

if(i> i1) then
u(i) =0.
else
u(i)=5.
end if
ui(i)=u(i)
end do

count=0

do n=2,nm
if (dabs(delt*dble(n-1)-int) <= 1.d-7) then

order=1
elseif (dabs(delt*dble(n-1)-int*2.) <= 1.d-7) then
order=1
elseif (dabs(delt*dble(n-1)-int*3.) <= 1.d-7) then
order=1
elseif (dabs(delt*dble(n-1)-int*4.) <= 1.d-7) then
order=1
elseif (dabs(delt*dble(n-1)-int*5.) <= 1.d-7) then
order=1
elseif (dabs(delt*dble(n-1)-int*6.) <= 1.d-7) then
order=1
elseif (dabs(delt*dble(n-1)-int*7.) <= 1.d-7) then
order=1
else
	order=0
end if

do i=1,im
e(i)=u(i)*u(i)/2.	
a(i)=u(i)
enddo

select case(opt)
case(1)
call lm(e,u)
case(2)
call lwm(eps,e,u)
case(3)
call mmm(e,u)
case(4)
call bwim(a,e,u)

end select

if(order ==1) then
count=count+1

do i=1,im
if (count ==1) then
u1(i)=u(i)
elseif(count==2) then
u2(i)= u(i)
elseif(count==3) then
u3(i) = u(i)
elseif(count==4) then
u4(i) = u(i)
elseif(count==5) then
u5(i) = u(i)
elseif(count==6) then
u6(i) = u(i)
elseif(count==7) then
u7(i) = u(i)
endif

end do
endif

do i=1,im
ut(i,n)=u(i)
end do


end do

call system_clock( t2 )
print *, "wall time in ms: ", ( t2 - t1 )*1000._PREC / cr
select case(opt)

case(1)
open(unit=9, file='lax.dat', status='replace', &
	action='write')
open(unit=19, file='tplax.dat', status='replace', &
	action='write')
case(2)
open(unit=9, file='lw.dat', status='replace', &
	action='write')
open(unit=19, file='tplw.dat', status='replace', &
	action='write')
case(3)
open(unit=9, file='mac.dat', status='replace', &
	action='write')
open(unit=19, file='mactp.dat', status='replace', &
	action='write')
case(4)
open(unit=9, file='beam.dat', status='replace', &
	action='write')
open(unit=19, file='tpbeam.dat', status='replace', &
	action='write')

end select

write(9,*)
write(9,10)
write(9,*) 
write(19,*) 'title= "nonlinear solution"'
write(19,*) 'variables= "x", "t", "u" '
write(19,*) 'zone F=point, I=', im, 'J=', nm

do i=1,im
write(9,20) x(i), ui(i), u1(i), u2(i), u3(i), u4(i), u5(i), u6(i)

enddo

do n=1,nm
do i=1,im
write(19,30) x(i), delt*dble(n-1), ut(i,n)

enddo
enddo

10 format (1x,'x', 3x, 't=0.0',3x,'t=0.4', 3x, 't=0.8', 3x, 't=1.2', 3x, 't=1.6', 3x, 't=2.0', 3x, 't=2.4')
20 format (1x, d15.8,3x, d15.8,6(3x, d15.8))
30 format (1x, d15.8, 3x, d15.8, 3x, d15.8)

close(9)
close(19)
end program burger
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine errh(sta, end)
implicit none
integer :: sta, end

print*, "wrong key entered!"
print*, "re-enter key", sta, "-", end
print*, " "
end subroutine errh


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lax method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lm(e,u)
use prop2
implicit none
integer :: i
real(8), dimension(imax) :: u,uold,e

do i=1,im
uold(i)=u(i)
enddo

do i=2,im-1
u(i)= 0.5 *(uold(i+1) +uold(i-1)) - c/2*(e(i+1)-e(i-1))

enddo

return
end subroutine lm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  LAX-WENDROFF method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lwm(eps,e,u)

use prop2
implicit none

integer :: i
real(8), dimension (imax) :: u, uold, e
real(8) 				  :: eps

do i=1,im
uold(i)=u(i)

enddo

do i=2,im-1

u(i)= uold(i)-(c/2.)*(e(i+1)-e(i-1))  &
	+ c* c/4 *((uold(i+1)+uold(i))*(e(i+1)-e(i)) - &
			(uold(i)+uold(i-1))*(e(i)-e(i-1)) ) + &
	 eps * (uold(i+1)-2*uold(i)+uold(i-1))
enddo

return
end subroutine lwm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MacCormack method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mmm(e,u)
use prop2
implicit none
integer :: i
real(8), dimension (imax) :: ustar, u, e
real(8) :: estari, estarim


do i=1,im-1

ustar(i) = u(i) - c* (e(i+1)-e(i))

enddo

do i=2,im-1
estari= ustar(i)*ustar(i)/2.
estarim= ustar (i-1) *ustar(i-1)/2.
u(i) =0.5 * (( u(i) + ustar(i)) -c * (estari-estarim))

enddo

return
end subroutine mmm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Beam and Warming method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bwim(a,e,u)

use prop2
implicit none

integer :: i
real(8), dimension (imax) :: aa, bb,cc, dd
real(8), dimension(imax)  :: u, a, e

do i=2,im-1
aa(i)= -0.25*c*a(i-1)
bb(i)= 1.
cc(i)= 0.25*c*a(i+1)
dd(i)= u(i) - 0.5*c*( e(i+1) - e(i-1)) + c/4*a(i+1)*u(i+1) - &
		c/4 * a(i-1)* u(i-1)
enddo

call trid(aa,bb,cc,dd,u)
return
end subroutine bwim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tridiagonal sover
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine trid(aa,bb,cc,dd,u)
use prop2
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

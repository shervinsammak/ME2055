! main program
program incompns
implicit none
integer,parameter :: imax=100
integer,parameter :: jmax=300

integer im,jm,k,i,j
real*8 psi(imax,jmax)
real*8 w(imax,jmax)
real*8 u(imax,jmax),v(imax,jmax)
real*8 check(imax,jmax)
real*8 dx,dy,duv
real*8 dt
real*8 nu
real*8 l,h
real*8 congs,conss
integer kmax
integer :: t1, t2,cr           ! timing variables

im=61
jm=61

l=12.0
h=3.0
dx=l/dfloat(im-1)
dy=h/dfloat(jm-1)
dt=0.001
nu=0.0025
congs=0.001
conss=0.002

! use only two threads
  !$ call omp_set_num_threads(4)
  
call system_clock( t1,cr )


do i=1,im
	do j=1,jm
		check(i,j)=0.0
	enddo
enddo

do i=1,10
	do j=1,10
		check(i,j)=1.0
	enddo
enddo


kmax=100
k=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variable intialization
call initial(im,jm,psi,w,u,v)

5	continue

!	do i=1,im
!		do j=1,jm
!			if( check(i,j).eq.0.0) then 
k=k+1

call boundary_condition(im,jm,dx,dy,psi,w,u,v)

call update_omega(im,jm, dx, dy, w, u, v)


call ftcs(im,jm, dx, dy,dt, w, u, v, nu)

call pgs(im,jm, dx, dy, psi, w, congs)

call boundary_condition(im,jm,dx,dy,psi,w,u,v)

call velocity(im,jm, dx, dy, psi, u, v, duv)

!			endif
!		enddo
!	enddo
	
	
call system_clock( t2,cr )
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!display variation of velocity 

write(*,15) k, duv
15	format(2x, 'at iteration', i6, ', duv= ', f10.5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! checking convergence

if(k.gt.kmax) then
	write(*,10) kmax
10	format (2x, 'no convegence reached in', i6, 'iteration')
else if (duv.gt.conss) then
	goto 5
	
endif

write(*,*) 'cpu_time in ms=', (t2-t1)*1000/cr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write reasult in files

call write_in_file(im,jm, dx, dy, psi, w, u, v)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end program incompns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initial subroutine

subroutine initial(im,jm, psi, w, u, v)
implicit none
integer,parameter :: imax=100
integer,parameter :: jmax=300

integer im,jm,k, i, j
real*8 psi(imax,jmax)
real*8 w(imax,jmax)
real*8 u(imax,jmax),v(imax,jmax)
real*8 dx,dy
real*8 dt
real*8 nu
real*8 l,h
real*8 congs,conss
integer kmax

	!$omp parallel default(none) &
    !$omp private(i,j)
     
    !$omp do

do i=1,im
	do j=1, jm
	
		psi(i,j)= 0.0
		w(i,j)= 0.0
		u(i,j)= 0.0
		v(i,j)= 0.0
	enddo
enddo

	!$omp end do
    !$omp end parallel
	
return
end subroutine initial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary condition subroutine

subroutine boundary_condition(im, jm, dx, dy, psi, w, u, v)
implicit none
integer,parameter :: imax=100
integer,parameter :: jmax=300

integer im,jm,k, i, j
real*8 psi(imax,jmax)
real*8 w(imax,jmax)
real*8 u(imax,jmax),v(imax,jmax)
real*8 dx,dy
real*8 dt
real*8 nu
real*8 l,h
real*8 congs,conss
integer kmax
real*8 :: psi_0, psi_20, psi_25, psi_30
integer i1,i2, i3, i4, j1, j2, j3
real*8 dx2, dy2

dx2=dx*dx
dy2=dy*dy

i1=10
i2=11
i3=16
i4=20

j1=20
j2=6
j3=21

psi_0=0.0
psi_20=20.0
psi_25=25.0
psi_30=10.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary condition for psi

do i=1, i1
	psi(i,j1)= psi_0
enddo

do j=1, j1
	psi(i1,j)= psi_0
enddo

do i=i1+1, im
	psi(i,1)= psi_0
enddo

do i=1,im
	psi(i, jm)=psi_30
enddo

do j=j1,jm
	psi(1,j)= psi(2,j)
enddo


do j=1,jm
	psi(im, j)= psi(im-1,j)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary condition for vorticity w
do i=1,i1
	w(i,j1)= 2.0*(psi(i,j1)-psi(i,j1+1))/dy2
enddo

do j=1,j1
	w(1,j)= 2.0*(psi(i1,j)-psi(i1+1,j))/dx2
enddo

do i=i1+1,im
	w(i,j1)= 2.0*(psi(i,1)-psi(i,2))/dy2
enddo



do i=1,im
	w(i, jm)=2.0*(psi(i,jm)-psi(i,jm-1))/dy2
enddo


do j=j1,jm-1
	w(1,j)= 2.0*(psi(1,j)-psi(2,j))/dx2 &
		-(psi(1,j+1)-2.0*psi(1,j)+psi(1,j-1))*dy2
enddo


do j=j3,jm
	w(im, j)= 2.0*(psi(im,j)-psi(im-1,j))/dx2 &
		-(psi(im,j+1)-2.0*psi(im,j)+psi(im,j-1))*dy2
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary condition for velocity 

do i=1,im-1
		u(i,1)=0.0
		v(i,1)=0.0
enddo

do i=1,i1
		u(i,j1)=0.0
		v(i,j1)=0.0
enddo

do j=1,j1
		u(i1,j)=0.0
		v(i1,j)=0.0
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,im
		u(i,jm)=0.0
		v(i,jm)=0.0
enddo


do j=j1,jm-1
		u(i,1)=-(psi(1,j1+1)-psi(1,j1))/dx
		v(i,1)=0.0
enddo



do j=1,jm-1
		u(i,1)=-(psi(im,j+1)-psi(im,j))/dx
		v(i,1)=0.0
enddo

u(im,jm)=0.0
v(im,jm)=0.0

return

end subroutine boundary_condition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine update_omega

subroutine update_omega(im,jm,dx,dy,w,u,v)
implicit none
integer,parameter :: imax=100
integer,parameter :: jmax=300

integer im,jm,k, i, j
real*8 psi(imax,jmax)
real*8 w(imax,jmax)
real*8 u(imax,jmax),v(imax,jmax)
real*8 check(imax,jmax)
real*8 dx,dy
real*8 dt
real*8 nu
real*8 l,h
real*8 congs,conss
integer kmax
real*8 dvdx, dudy

	
	
	!$omp parallel default(none) &
    !$omp private(i,j)
     
    !$omp do
	
do i=2,im-1
	do j=2,jm-1
		if(check(i,j) .eq. 0.0) then
		dvdx=(v(i+1,j)-v(i-1,j))/2.0/dx
		dudy=(u(i,j+1)-u(i,j-1))/2.0/dy
		w(i,j)= dvdx-dudy
		
		endif
	enddo
enddo

	!$omp end do
	!$omp end parallel
	
return
end subroutine update_omega
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine ftcs for vorticity unsteady equation

subroutine ftcs(im,jm,dy,dt,w,u,v,nu)
implicit none
integer,parameter :: imax=100
integer,parameter :: jmax=300

integer im,jm,k, i, j
real*8 psi(imax,jmax)
real*8 w(imax,jmax)
real*8 u(imax,jmax),v(imax,jmax)
real*8 check(imax,jmax)
real*8 dx,dy
real*8 dt
real*8 nu
real*8 l,h
real*8 congs,conss
integer kmax
real*8 dw(imax,jmax)
real*8 dx2,dy2,dx22,dy22
real*8 udwdx, vdwdy, d2wdx2, d2wdy2

dx2=dx*dx
dy2=dy*dy
dx22=2.0*dx
dy22=2.0*dy

	!$omp parallel default(none) &
    !$omp private(i,j)
     
    !$omp do
do i=2,im-1
	do j=2,jm-1
		if(check(i,j).eq. 0.0) then
	
		udwdx= u(i,j)*(w(i+1,j)-w(i-1,j))/dx22
		vdwdy= v(i,j)*(w(i,j+1)-w(i,j-1))/dy22
		
		d2wdx2= (w(i+1,j)-2.0*w(i,j)+w(i-1,j))/dx2
		d2wdy2= (w(i,j+1)-2.0*w(i,j)+w(i,j-1))/dy2
		
		dw(i,j)=dt*(-(udwdx+vdwdy)+nu*(d2wdx2+d2wdy2))
		
		endif
	enddo
enddo

	!$omp end do
	!$omp end parallel
	
	
	!$omp parallel default(none) &
    !$omp private(i,j)
     
    !$omp do
do i=2,im-1
	do j=2,jm-1
		if(check(i,j).eq. 0.0) then
		w(i,j)= w(i,j)+dw(i,j)
		endif
	enddo
enddo

	!$omp end do
	!$omp end parallel
return
end subroutine ftcs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine point gauss-seidel

subroutine pgs(im,jm,dx,dy,psi,w,congs)
implicit none
integer,parameter :: imax=100
integer,parameter :: jmax=300

integer im,jm,k, i, j
real*8 psi(imax,jmax)
real*8 w(imax,jmax)
real*8 u(imax,jmax),v(imax,jmax)
real*8 check(imax,jmax)
real*8 dx,dy
real*8 dt
real*8 nu
real*8 l,h
real*8 congs,conss
integer kmax
real*8 psi_old,error,errinit,err

real*8 :: dx2,dy2,beta2,cxy

dx2=dx*dx
dy2=dy*dy
beta2=(dx2/dy2)
cxy=0.5/(1.0+beta2)
k=0

5 continue 

k=k+1
error=0.0
	
	!$omp parallel default(none) &
    !$omp private(i,j)
     
    !$omp do

do i=2,im-1
	do j=2,jm-1
		if(check(i,j).eq. 0.0) then
		psi_old=psi(i,j)
		psi(i,j)=cxy*(dx2*w(i,j)+psi(i+1,j)+psi(i-1,j) &
			+beta2*(psi(i,j+1)+psi(i,j-1)))
		error=error+abs(psi(i,j)-psi_old)
		
		if(k.eq.1) then
		errinit=error
		endif
		
		endif
	enddo
enddo
	!$omp end do
	!$omp end parallel
	
err=error/errinit

if(k.gt.10000) then 
	write(*,100) 
else
	if(err.gt.congs) then
	goto 5
	endif
endif

100 format (2x, 'pgs did not converge in 10000 iterations')

return

end subroutine pgs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine velocity

subroutine velocity(im,jm,dx,dy,psi,u,v,duv)
implicit none
integer,parameter :: imax=100
integer,parameter :: jmax=300

integer :: im,jm,k, i, j
real*8 psi(imax,jmax)
real*8 w(imax,jmax)
real*8 u(imax,jmax),v(imax,jmax)
real*8 check(imax,jmax)
real*8 dx,dy,duv
real*8 dt
real*8 nu
real*8 l,h
real*8 congs,conss
integer kmax
real*8 	 unew,vnew
real*8 	 num,denum

num=0.0
denum=0.0

	!$omp parallel default(none) &
    !$omp private(i,j)
     
    !$omp do
do i=2,im-1
	do j=2,jm-1
		if(check(i,j).eq. 0.0) then
		
		unew= (psi(i,j+1)-psi(i,j-1))/2.0/dy
		vnew= (psi(i+1,j)-psi(i-1,j))/2.0/dx
		
		num=num+ sqrt((unew-u(1,j))**2.0+ (vnew-v(i,j))**2.0)
		
		u(i,j)=unew
		v(i,j)=vnew
		
		endif
	enddo
enddo
	!$omp end do
	!$omp end parallel
	
	!$omp parallel default(none) &
    !$omp private(i,j)
     
    !$omp do
do i=2,im-1
	do j=2,jm-1
		if(check(i,j).eq. 0.0) then
		denum=denum+ sqrt(u(i,j)**2.0+v(i,j)**2.0)
		
		endif
	enddo
enddo
	!$omp end do
	!$omp end parallel
duv=num/denum
return

end subroutine velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine write_in_file
subroutine write_in_file(im,jm,dx,dy,psi,w,u,v)
implicit none
integer,parameter :: imax=100
integer,parameter :: jmax=300

integer im,jm,k, i, j
real*8 psi(imax,jmax)
real*8 w(imax,jmax)
real*8 u(imax,jmax),v(imax,jmax)
real*8 check(imax,jmax)
real*8 dx,dy
real*8 dt
real*8 nu
real*8 l,h
real*8 congs,conss
integer kmax
open(1,file='pbs8.dat')
write(1,111)im,jm,dx,dy
111 format ('im=', i4, ' ; jm=',  i4, '; dx=', f5.3, ';dy=', f5.3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(1,*) 'stream function distribuation'
write(1,*) 
write(1,100) ((i-1)*dx, i=1,im,5)

do j=jm,1,-1
	write(1,200) (j-1)*dy, (psi(i,j), i=1,im,5)
enddo

write(1,*)
write(1,*)
write(1,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(1,*) 'vorticity distribution'
write(1,*)
write(1,100) ((i-1)*dx, i=1,im,5)

do j=jm,1,-1
	write(1,200) (j-1)*dy, (w(i,j), i=1,im,5)
enddo

write(1,*)
write(1,*)
write(1,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(1,*) 'u component of velocity'
write(1,*)
write(1,100) ((i-1)*dx, i=1,im,5)

do j=jm,1,-1
	write(1,200) (j-1)*dy, (u(i,j), i=1,im,5)
enddo

write(1,*)
write(1,*)
write(1,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(1,*) 'v component of velocity'
write(1,*)
write(1,100) ((i-1)*dx, i=1,im,5)

do j=jm,1,-1
	write(1,200) (j-1)*dy, (v(i,j), i=1,im,5)
enddo

write(1,*)
write(1,*)
write(1,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

100 format(14x, 'y', 11('x=', f3.1, 6x))
200 format(f17.3, 2x, 11('x=', f3.1, 6x))
close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(2, file='pbs8.dat')
write(2,221)
write(2,222)
write(2,223) im, jm
221 format(1x, 'title="psi,w,u,v,velocity"')
222 format(1x, 'variables= "x", "y", "psi", "w", "u", "v","velocity" ')
223 format(1x, 'zone t= "l", i=', i3, 2x, 'j=', i3,2x, 'f=point')

do  j=1,jm
	do i=1,im
		if(check(i,j) .eq. 0.0) then
		write(2,224) (i-1)*dx, (j-1)*dy, psi(i,j), w(i,j), u(i,j), v(i,j), sqrt(u(i,j)**2+v(i,j)**2)
		
		endif
	enddo
enddo

224 format (7(2x,d14.8))
close(2)

return

end subroutine write_in_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
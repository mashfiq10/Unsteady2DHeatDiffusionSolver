!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
! MAE 6263 Computational Fluid Dynamics - Project 1 !
! Sk. Mashfiqur Rahman !
! CWID: A20102717 !
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

program heatdiff
implicit none
integer::i,j,k,nx,ny,nt,nitr,N,ischeme,c,pr
real*8 ::dx,dt,dy,xi,xf,yi,t,yf,alpha,Tmax,gamma
real*8 ::tol,omega,w,t1,t2,probe
real*8,allocatable ::Te(:,:),Ti(:,:),x(:),y(:),b(:,:)

!Domain
xi = 0.0d0 !left
xf = 1.0d0 !right

yi = 0.0d0 !bottom
yf = 1.0d0 !up

open(9,file='project1input.txt')
read(9,*)nx		    !resolution in x direction
read(9,*)ny	    	!resolution in y direction
read(9,*)gamma		!value of gamma to measure time-step size
read(9,*)Tmax		!final time
read(9,*)alpha      !diffusion coefficient
read(9,*)tol        !tolerance for iterative solvers
read(9,*)omega      !relaxation parameter for SOR iterative solver (optimized value = 1.08-1.28)
read(9,*)probe      !value to calculate data on a point in space
read(9,*)ischeme	![0]Explicit, [1]Implicit-Gauss-Seidel(GS) solver, [2]Implicit-SOR, [3]Both explicit and implicit (GS) (for analysis), [4]Both explicit and implicit (SOR)(for analysis), [5]Approximate factorization
close(9)


!Calculating grid spacing (spatial)
dx = (xf-xi)/dfloat(nx)
dy = (yf-yi)/dfloat(ny)
pr = dint(((probe - yi)*dfloat(ny))/(yf-yi))

!Time step
dt = ((dx*dx)*gamma) /alpha

!calculating number of snapshots (ns)
nt = nint(Tmax/dt)

!spatial coordinate points 
allocate(x(0:nx))
do i=0,nx
x(i) = xi + dfloat(i)*dx
end do

allocate(y(0:ny))
do j=0,ny
y(j) = yi + dfloat(j)*dy
end do

!T: temperature variable 
!b: solution vector for Ax=b
allocate(Te(0:nx,0:ny))
allocate(Ti(0:nx,0:ny))
allocate(b(0:nx,0:ny))

!initial and boundary condition
!N: iteration counter
!c: time step counter
t = 0.0d0
c = 0
N = 0

do j=0,ny
do i=0,nx
Te(i,j) = 0.0d0
end do
end do

do j=0,ny
Te(0,j)  = 0.0d0 +y(j)
Te(nx,j) = 1.0d0 +y(j)
end do

do i=0,nx
Te(i,0)  = x(i)+ 0.0d0
Te(i,ny) = x(i)+ 1.0d0
end do

do j=0,ny
do i=0,nx
    Ti(i,j) = Te(i,j)
end do
end do


!Plot initial condition
open(18,file='Tin.plt')
write(18,*) 'variables ="x","y","T"'
write(18,100)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
do j=0,ny
do i=0,nx
write(18,*) x(i),y(j),Te(i,j)
end do
end do
close(18)

!$$$$$$ open(20,file='nitr.plt')
!$$$$$$ write(20,*) 'variables ="w","N-ITR"'

if(ischeme.eq.1 .or. ischeme.eq.3) then
  w = 1.0d0
else
  w = omega
end if
  
call cpu_time(t1)

!-----------------------------------------------------------------------------!
!Time integration (for both explicit and implicit solvers)
!-----------------------------------------------------------------------------!
if (ischeme.eq.0) then !explicit solver
  
do k=1,nt   
  
    c = c+1 !time step counter
    call FTCS(nx,ny,Te,gamma)
    t = t + dt !time update
end do

else if(ischeme.eq.3 .or. ischeme.eq.4) then

do k=1,nt   
  
    c = c+1 !time step counter
    call FTCS(nx,ny,Te,gamma)
    t = t + dt !time update
end do

do k=1,nt

    c = c+1 !time step counter
    
    !solution vector update
    do j=0,ny
    do i=0,nx
       b(i,j) = Ti(i,j)
    end do
    end do

    call BTCS(nx,ny,Ti,b,w,gamma,tol,nitr)
    t = t + dt !time update for plotting
    N = N + nitr  
    
end do

else if(ischeme.eq.5) then !Approximate factorization
  
do k=1,nt

    c = c+1 !time step counter
    !solution vector update
    do j=0,ny
    do i=0,nx
       b(i,j) = Ti(i,j)
    end do
    end do

    call BTCSAF(nx,ny,gamma,Ti,b)
    t = t + dt !time update
    
end do

else  !implicit solver
  
do k=1,nt

    c = c+1 !time step counter
    !solution vector update
    do j=0,ny
    do i=0,nx
       b(i,j) = Ti(i,j)
    end do
    end do

    call BTCS(nx,ny,Ti,b,w,gamma,tol,nitr)
    t = t + dt !time update
    N = N + nitr  
    
end do
end if

!call OptOm(nx,ny,Th,b,gamma,tol,omega)  


!-----------------------------------------------------------------------------!
!output to .plt file for Tecplot
!-----------------------------------------------------------------------------!

call cpu_time(t2)

open(4,file='cpu.txt')
write(4,*)"cpu time (sec)=",(t2-t1)
close(4)

print*,"------------------------"
print*,"Total CPU time (seconds) = ", (t2-t1)
print*,"------------------------"

if (ischeme.eq.0) then    !FTCS
  
open(12, file="explicit.plt")
write(12,*)'variables ="x","y","T"'
write(12,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
do j=0,ny
do i=0,nx
	write(12,*) x(i),y(j),Te(i,j)
end do
end do
close(12) 

open(22, file="explicit.txt")
do j=0,ny
do i=0,nx
	write(22,*) x(i),y(j),Te(i,j)
end do
end do
close(22) 

print*,"-------------------------------"
print*,"Total number of time iterations = ", c
print*,"-------------------------------"

else if(ischeme.eq.3 .or. ischeme.eq.4) then   !GS + SOR

open(19, file="explicit.plt")
write(19,*)'variables ="x","y","T"'
write(19,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
do j=0,ny
do i=0,nx
	write(19,*) x(i),y(j),Te(i,j)
end do
end do
close(19)

open(29, file="explicit.txt")
do j=0,ny
do i=0,nx
	write(29,*) x(i),y(j),Te(i,j)
end do
end do
close(29)

open(14, file="implicit.plt")
write(14,*)'variables ="x","y","T"'
write(14,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
do j=0,ny
do i=0,nx
	write(14,*) x(i),y(j),Ti(i,j)
end do
end do
close(14)

open(24, file="implicit.txt")
do j=0,ny
do i=0,nx
	write(24,*) x(i),y(j),Ti(i,j)
end do
end do
close(24)

open(12, file="error_iso.plt")
write(12,*)'variables ="x","y","error"'
write(12,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
do j=0,ny
do i=0,nx
	write(12,*) x(i),y(j),dabs(Te(i,j)-Ti(i,j))
end do
end do
close(12)

open(22, file="error_iso.txt")
do j=0,ny
do i=0,nx
	write(22,*) x(i),y(j),dabs(Te(i,j)-Ti(i,j))
end do
end do
close(22)

open(3, file="impact_analysis_ex.plt")
write(3,*)'variables ="x","T"'

do i=0,nx
	write(3,*) x(i),Te(i,pr)
end do
close(3)

open(23, file="impact_analysis_ex.txt")
write(23,*)'variables ="x","T"'

do i=0,nx
	write(23,*) x(i),Te(i,pr)
end do
close(23)

open(4, file="impact_analysis_im.plt")
write(4,*)'variables ="x","T"'

do i=0,nx
	write(4,*) x(i),Ti(i,pr)
end do
close(4)

open(34, file="impact_analysis_im.txt")
write(34,*)'variables ="x","T"'

do i=0,nx
	write(34,*) x(i),Ti(i,pr)
end do
close(34)

else if(ischeme.eq.5) then    !Approximate Factorization

!$$$$$$ print*,"------------------------"
!$$$$$$ print*,"Total CPU time (seconds) = ", (t2-t1)
!$$$$$$ print*,"------------------------"

open(12, file="implicit.plt")
write(12,*)'variables ="x","y","T"'
write(12,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
do j=0,ny
do i=0,nx
	write(12,*) x(i),y(j),Ti(i,j)
end do
end do
close(12)

open(22, file="implicit.txt")
do j=0,ny
do i=0,nx
	write(22,*) x(i),y(j),Ti(i,j)
end do
end do
close(22)

open(4, file="impact_analysis_im.plt")
write(4,*)'variables ="x","T"'

do i=0,nx
	write(4,*) x(i),Ti(i,10)
end do
close(4)


open(34, file="impact_analysis_im.txt")
write(34,*)'variables ="x","T"'

do i=0,nx
	write(34,*) x(i),Ti(i,pr)
end do
close(34)

print*,"-------------------------------"
print*,"Total number of time iterations = ", c
print*,"-------------------------------"

else           !implicits

open(12, file="implicit.plt")
write(12,*)'variables ="x","y","T"'
write(12,*)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
do j=0,ny
do i=0,nx
	write(12,*) x(i),y(j),Ti(i,j)
end do
end do
close(12)

open(22, file="implicit.txt")
do j=0,ny
do i=0,nx
	write(22,*) x(i),y(j),Ti(i,j)
end do
end do
close(22)

open(4, file="impact_analysis_im.plt")
write(4,*)'variables ="x","T"'

do i=0,nx
	write(4,*) x(i),Ti(i,pr)
end do
close(4)


open(34, file="impact_analysis_im.txt")
write(34,*)'variables ="x","T"'

do i=0,nx
	write(34,*) x(i),Ti(i,pr)
end do
close(34)

print*,"----------------------------------------------"
print*,"Total number of iterations in iterative solver = ", N
print*,"----------------------------------------------"


print*,"-------------------------------"
print*,"Total number of time iterations = ", c
print*,"-------------------------------"
!$$$$$$ 
!$$$$$$ print*,"------------------------"
!$$$$$$ print*,"Total CPU time (seconds) = ", (t2-t1)
!$$$$$$ print*,"------------------------"

end if

100 format(a16,i8,a4,i8,a10,f10.4,a3)
end

!-----------------------------------------------------------------------------!
!Explicit Solver (spatial treatment)
!-----------------------------------------------------------------------------!
subroutine FTCS(nx,ny,Te,gamma)
implicit none
integer::nx,ny,i,j
real*8 ::gamma
real*8 ::u(0:nx,0:ny),Te(0:nx,0:ny)

!previous step (t=n)
do j=0,ny
do i=0,nx
u(i,j) = Te(i,j)
end do
end do

!update (t=n+1)
do j=1,ny-1
do i=1,nx-1
Te(i,j) = u(i,j) + gamma*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4.0d0*u(i,j))
end do
end do

end 

!-----------------------------------------------------------------------------!
!Iterative solver (for implicit scheme) (spatial treatment)
!-----------------------------------------------------------------------------!
subroutine BTCS(nx,ny,Ti,b,w,gamma,tol,nitr)  
implicit none
integer::nx,ny,i,j,nitr
real*8 ::tol,err,w,gamma
real*8 ::Ti(0:nx,0:ny),e(0:nx,0:ny),v(0:nx,0:ny),b(0:nx,0:ny)

err=1.0d0
nitr = 0
   
do while(err.gt.tol)

  nitr = nitr + 1
  err=0.0d0
    do j=0,ny
	do i=0,nx
	v(i,j) = Ti(i,j)
	end do
	end do

	!update
    do j=1,ny-1
	do i=1,nx-1
	Ti(i,j) = v(i,j) + ((w*gamma)/(1 + (4.0d0*gamma)))*(v(i+1,j)+Ti(i-1,j)+v(i,j+1)+Ti(i,j-1)+ &
              (b(i,j)/gamma)-((1 + (4.0d0*gamma))*v(i,j)/gamma))
	end do
	end do

    !compute L1 norm
    do j=0,ny
    do i=0,nx
    e(i,j) = dabs(Ti(i,j)-v(i,j))
    err =  e(i,j) + err
    end do
    end do 

end do
end 

!-----------------------------------------------------------------------------!
!Optimized Omega calculation
!-----------------------------------------------------------------------------!
subroutine OptOm(nx,ny,Ti,b,gamma,tol,omega)  
implicit none
integer::i,j,nx,ny,nitr,p
real*8::dom,om1,nom,gamma,tol,omega
real*8 ::b(0:nx,0:ny),Ti(0:nx,0:ny)


om1=0.5d0
dom = 0.05d0
nom = 20

do p=0,nom
    do j=0,ny
    do i=0,nx
       b(i,j) = Ti(i,j)
    end do
    end do

  	omega = om1 + dfloat(p)*dom
    
    call SOR(nx,ny,Ti,b,omega,gamma,tol,nitr)

       
	!write error
    write(*,*)'omega and total interval:'
    write(20,*) omega, nitr
    write(*,*) omega, nitr    
end do
close(20)
end


!-----------------------------------------------------------------------------!
!Approximate Factorization
!-----------------------------------------------------------------------------!
subroutine BTCSAF(nx,ny,gamma,Ti,b)
implicit none
integer::nx,ny,i,j
real*8 ::gamma,bx,by,Tia,Tib
real*8 ::Ti(0:nx,0:ny)
real*8,allocatable ::t(:),d(:),m(:),r(:),q(:)
real*8 ::b(0:nx,0:ny),z(0:nx,0:ny)

do j=0,ny
do i=0,nx
z(i,j) = 0.0d0
end do
end do

bx = gamma
by = gamma

!x-sweep to compute intermediate values:
do j=1,ny-1
 
	!Build coefficient matrix:
	allocate(t(1:nx-1),d(1:nx-1),m(1:nx-1),r(1:nx-1),q(1:nx-1))

	do i=1,nx-1
	t(i) = -bx
	d(i) = (1.0d0+2.0d0*bx)
	m(i) = -bx
	r(i) = b(i,j) 
	end do
    
	!apply boundary conditions
    Tia = Ti(0,j)  - by*(Ti(0,j+1)-2.0d0*Ti(0,j)+Ti(0,j-1))
    Tib = Ti(nx,j) - by*(Ti(nx,j+1)-2.0d0*Ti(nx,j)+Ti(nx,j-1))
    
	r(1)   = r(1) - t(1)*Tia        !b.c.
	r(nx-1) = r(nx-1) - m(nx-1)*Tib !b.c.
    
	call tdma(t,d,m,r,q,1,nx-1)

	!assign solutions for as z
	do i=1,nx-1
	z(i,j)=q(i)
	end do
    z(0,j) =Tia
    z(nx,j)=Tib

    deallocate(t,d,m,r,q)

end do


!y-sweep to compute final solution:
do i=1,nx-1

	!Build coefficient matrix:
	allocate(t(1:ny-1),d(1:ny-1),m(1:ny-1),r(1:ny-1),q(1:ny-1))

	do j=1,ny-1
	t(j) = -by
	d(j) = (1.0d0+2.0d0*by)
	m(j) = -by
	r(j) = z(i,j) 
	end do
    
	!apply boundary conditions
    Tia = Ti(i,0)  
    Tib = Ti(i,ny) 
    
	r(1)   = r(1) - t(1)*Tia        !b.c.
	r(ny-1) = r(ny-1) - m(ny-1)*Tib !b.c.
    
	call tdma(t,d,m,r,q,1,ny-1)

	!assign solutions for as z
	do j=1,ny-1
	Ti(i,j)=q(j)
	end do

    deallocate(t,d,m,r,q)
end do

end 

!------------------------------------------------------------------!
!TDMA
!------------------------------------------------------------------!
subroutine tdma(t,d,m,r,x,k,l)
implicit none
integer k,l,i
real*8, dimension(k:l) ::t,d,m,r,x    

! forward elimination phase
do i=k+1,l
d(i) = d(i) -t(i)/d(i-1)*m(i-1)
r(i) = r(i) - t(i)/d(i-1)*r(i-1)
end do
! backward substitution phase 
x(l) = r(l)/d(l)
do i=l-1,k,-1
x(i) = (r(i)-m(i)*x(i+1))/d(i)
end do

return
end

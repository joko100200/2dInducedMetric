!Whenever you edit this code and update it please put your name infront of the code and what you did and why you did it :) pretty
!please. :)))
!
!People:
!Jonathan Sar-Shalom
!
!




!This is the main code
program main
implicit none

!Variable declerations
!mesh is a 2d array of size t,12
real*8, dimension(:,:), allocatable :: mesh,uv,xyz1
integer :: error, i, ut, vt, t, a, b,switch
real*8 :: u,v,hu,hv,umax,umin,vmax,vmin
real*8 :: xpu,ypu,zpu,xpv,ypv,zpv,invdetg
namelist /inputs/ hu,hv,umax,umin,vmax,vmin,switch

!------Variable Explinations------!
!error is just to check if a file breaks
!i is the index for the mesh
!ut is U if you look at the word doc
!vt is V if you look at the word doc
!t is T if you look at the word doc
!u,v are the paramiterization of our manifold
!umax,umin,vmax,vmin look at word doc
!hu,hv look at word doc
!switch determines what set of functions we will be using


!Reads input data from file
open(10,file = 'parameters.inp',status = 'old',&
        action = 'read', position = 'rewind', iostat = error)
if(error .ne. 0) stop
read(10,inputs)
close(10)

if(umax < umin) then
        write(*,*) "umax smaller than umin please fix"
        stop
end if
if(vmax < vmin) then
        write(*,*) "vmax smaller than vmin please fix"
        stop
end if


!Calculations for our integer values
ut = int((umax-umin)/hu)
vt = int((vmax-vmin)/hv)
t = ut*vt

!Allocating memory for the size of mesh
allocate(mesh(t,12))
allocate(uv(t,2))
allocate(xyz1(t,3))

write(*,*) "Total columbs are: ",t

do i = 1, t , 1
call IndexToUV(i,ut,hu,hv,umin,vmin,u,v)
uv(i,1) = u
uv(i,2) = v
end do

open(11,file = 'uv.out',status = 'replace', &
        action = 'write', position = 'append', iostat = error)
if(error .ne. 0) stop


!writing the uv array to file uv
do i = 1, t, 1
write(11,*)uv(i,1),uv(i,2)
end do
close(11)

!calling xyz to use 2d array (uv) such to implant it into 3d
!size (t,3)
call xyz(switch,t,"uv.out",6,"xyzuv.out",9)


!Here we are writing out the 3d coords for partial derivative
!calculations
open(26,file='xyzuv.out',status='old',&
        action='read',position='rewind',iostat=error)
if(error .ne. 0) stop

do i = 1,t,1
read(26,*) xyz1(i,1),xyz1(i,2),xyz1(i,3)
enddo

close(26)



!first loop to calculate metric and inverse metric
!for all points i=[1,t]
do i = 1,t,1

!here the partial derivative will be defined all 6 and then they
!will be used to calculate the 3 metric components below :)

!!!!!!!!!!!!!!PROBLEM AREA!!!!!!!!!!!!!!!!

!if this condition is false we have made it to the end of th mesh
!in the U direction we will then skip the calculation and keep the
!last value of the derivatives
if(MOD(i,ut) .ne. 0) then
xpu = (xyz1(i+1,1)-xyz1(i,1))/hu   !dx/du
ypu = (xyz1(i+1,2)-xyz1(i,2))/hu   !dy/du
zpu = (xyz1(i+1,3)-xyz1(i,3))/hu   !dz/du
end if

!if this condition is true we are on the upper part of the mesh
!all derivatives are equal to the last ones since there are no
!more points "above"
if (i>(vt-1)*ut+1) then
xpv= (xyz1(i,1)-xyz1(i-ut,1))/hv   !dx/dv
ypv= (xyz1(i,2)-xyz1(i-ut,2))/hv   !dy/dv
zpv= (xyz1(i,3)-xyz1(i-ut,3))/hv   !dz/dv
else
xpv= (xyz1(i+ut,1)-xyz1(i,1))/hv   !dx/dv
ypv= (xyz1(i+ut,2)-xyz1(i,2))/hv   !dy/dv
zpv= (xyz1(i+ut,3)-xyz1(i,3))/hv   !dz/dv
end if
!!!!!!!!!!!!!!End PROBLEM AREA!!!!!!!!!!!!
mesh(i,1) = xpu*xpu+ypu*ypu+zpu*zpu !g11
mesh(i,2) = xpu*xpv+ypu*ypv+zpu*zpv !g21 g12
mesh(i,3) = xpv*xpv+ypv*ypv+zpv*zpv !g22
!now to calculate the inverse metric components
invdetg = 1/(mesh(i,1)*mesh(i,3)-mesh(i,2)*mesh(i,2)) 
mesh(i,4) = invdetg*mesh(i,3)  !g^11
mesh(i,5) = -invdetg*mesh(i,2) !g^21 g^12
mesh(i,6) = invdetg*mesh(i,1)  !g^22
enddo
deallocate(xyz1)

!second loop to calculate all connections for i=[1,t]

do i = 1, t, 1

!here we need to calculate the partial derivatives of the metric
!For simplisity I will reuse the previous variables for the metric
!derivatives

!!!!!!!!!!!!!!Problem Area!!!!!!!!!!!!!
if(MOD(i,ut) .ne. 0) then
xpu = (mesh(i+1,1)-mesh(i,1))/hu   !g11;1
ypu = (mesh(i+1,2)-mesh(i,2))/hu   !g12;1
zpu = (mesh(i+1,3)-mesh(i,3))/hu   !g22;1
end if

!if this condition is true we are on the upper part of the mesh
!all derivatives are equal to the last ones since there are no
!more points "above"
if (i>(vt-1)*ut+1) then
xpv= (mesh(i,1)-mesh(i-ut,1))/hv   !g11;2
ypv= (mesh(i,2)-mesh(i-ut,2))/hv   !g21;2
zpv= (mesh(i,3)-mesh(i-ut,3))/hv   !g22;2
else
xpv= (mesh(i+ut,1)-mesh(i,1))/hv   !g11;2
ypv= (mesh(i+ut,2)-mesh(i,2))/hv   !g12;2
zpv= (mesh(i+ut,3)-mesh(i,3))/hv   !g22;2
end if

!!!!!!!!!!!!!End Problem Area!!!!!!!!!!

mesh(i,7) = 0.5*(mesh(i,4)*xpu+mesh(i,5)*(2*ypu-xpv)) !111
mesh(i,8) = 0.5*(mesh(i,4)*xpv+mesh(i,5)*zpv) !121 & 112
mesh(i,9) = 0.5*(mesh(i,4)*(2*ypv-zpu)+mesh(i,5)*zpv) !122
mesh(i,10) = 0.5*(mesh(i,5)*xpu+mesh(i,6)*(2*ypu-xpv)) !211
mesh(i,11) = 0.5*(mesh(i,5)*xpv+mesh(i,6)*zpv) !212 & 221
mesh(i,12) = 0.5*(mesh(i,5)*(2*ypv-zpu)+mesh(i,6)*zpv) !222

enddo

open(12,file='connections.out',status='replace',&
        action='write',position='append',iostat=error)
if(error .ne. 0) stop


!third loop to write mesh(i,7-12) for i=[1,t]
do i = 1, t, 1
write(12,*) mesh(i,7),mesh(i,8),mesh(i,9),mesh(i,10),&
        mesh(i,11),mesh(i,12)
enddo
close(12)

open(37,file='metric.out',status='replace',&
        action='write',position='append',iostat=error)
if(error .ne. 0) stop

do i =1,t,1
write(37,*)mesh(i,1),mesh(i,2),mesh(i,3),mesh(i,4),&
        mesh(i,5),mesh(i,6)
enddo
close(37)

deallocate(uv)
deallocate(mesh)
end program main


!this subrutine takes in integer index i and returns
!the 2 real valued u,v paramiterization
subroutine IndexToUV(i,ut,hu,hv,umin,vmin,u,v)
        real*8, intent(out) :: u,v
        integer :: a,b
        integer, intent(in) :: ut
        real*8, intent(in) :: hu,hv,umin,vmin
        
        a = MOD(i-1,ut)
        b = int((i-1)/ut)
        u = a*hu+umin
        v = b*hv+vmin

        return
end subroutine IndexToUV

!This subrutine takes in the 2 u,v coords and outputs
!the index i
subroutine UVToIndex (i,u,v,ut,umin,vmin,hu,hv)
        integer, intent(out) :: i
        integer :: a,b
        integer, intent(in) :: ut
        real*8, intent(in) :: u,v,umin,vmin,hu,hv
        
        a = int((u-umin)/hu)
        b = int((v-vmin)/hv)

        i = b*ut+a+1

        return
end subroutine UVToIndex










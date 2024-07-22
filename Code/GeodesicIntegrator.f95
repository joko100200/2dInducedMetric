!This is the code where we are going to integrate

program integrator
implicit none
integer :: switch,time,i,error,T,l,j
integer :: ut,vt
real*8 :: u,v,hu,hv,umax,umin,vmax,vmin,ui,vi,uii,vii,utemp,vtemp
real*8 :: Au,Av,dt,f,maxvar,f1,f2,f3
real*8, dimension(2) :: Vu,Vv
real*8, dimension(:,:), allocatable :: uvpath, connection
namelist /inputs/ hu,hv,umax,umin,vmax,vmin,switch
namelist /geodesics/ ui,vi,uii,vii,time,dt,maxvar

!!!!!!!!!!!!Variable explination!!!!!!!!
!switch shouldn't be here lel
!time is the dummy variable for the integration loop
!u,v are dummy uv coordinates for the integration loop
!hu,hv encode how far away nodes are in the u and v respectivly
!umax,umin,vmax,vmin: are the total size of the grid in the u and v
!ui,vi are the initial coordinates for the integration
!uii,vii are the initial velocites for the integration
!Vu,Vv: are the velocity dummy variables
!Au,Av: are the acceleration dummy variables
!j,f are dummy variables


open(10,file = 'parameters.inp',status='old',&
action='read',position = 'rewind',iostat = error)
if(error .ne. 0) stop
read(10,inputs)
read(10,geodesics)
close(10)

ut = int((umax-umin)/hu)
vt = int((vmax-vmin)/hv)
T = ut*vt

allocate(connection(T,6))
allocate(uvpath(time,2))

open(1,file='connections.out',status='old',&
      action='read',position = 'rewind',iostat=error)
if(error .ne. 0) then
write(*,*) "No file of that name"
stop
end if

do i = 1, T, 1
read(1,*) connection(i,1:6)
enddo

close(1)


!Now there are many diffrent ways to integrate the one way that we know how is through
!the euler method. but is that really worth while to do here it seems like it is highly
!prone to error
!
!the whole grid is already calculated so we will be taking
!steps and reading what is at that location by finding it's
!clostest position and shoving it in there

!the way this is be done is by solving the geodesic
!for uv and then shoving those values into uvtoindex
!then taking that index and finding the metric and connection
!associated with that point and making the new uv point
!whatever indextouv says it is



! this is the main integration loop
!These are defining our initial conditions
Vu(1) = uii
Vv(1) = vii
!This procedure will correct most error build up
!hopefully
call UVtoIndex(i,ui,vi,ut,umin,vmin,hu,hv) 
call IndextoUV(i,ut,hu,hv,umin,vmin,ui,vi)
uvpath(1,1) = ui
uvpath(1,2) = vi




!!!!!!!!!!!!!!!!Problem we don't have a condition that stops the code from running off mesh
!!!!!!!!!!!!!!!!Please implement that 

do l = 2, time, 1
    !We don't need to have this line here since the i carries over this is a redudant calculations
    !call UVtoIndex(i,uvpath(l-1,1),uvpath(l-1,2),ut,umin,vmin,hu,hv)

    !l here is our time step variable where time is our max time which we input into parameters
    Au = -(connection(i,1)*Vu(1)*Vu(1)+2*connection(i,2)*Vu(1)*Vv(1)+connection(i,3)*Vv(1)*Vv(1))    
    Av = -(connection(i,4)*Vu(1)*Vu(1)+2*connection(i,5)*Vu(1)*Vv(1)+connection(i,6)*Vv(1)*Vv(1))
    !The second index is for the n+1 value while the first one is for the old velocity value
    Vu(2) = Au*dt+Vu(1) 
    Vv(2) = Av*dt+Vv(1)
    !we are using these dummy variables becasue we need to put them through the two call procedures
    u = Vu(2)*dt+uvpath(l-1,1)
    v = Vv(2)*dt+uvpath(l-1,2)
    call UVtoIndex(i,u,v,ut,umin,vmin,hu,hv)
    call IndextoUV(i,ut,hu,hv,umin,vmin,u,v)
    if (i .GT. T .or. i .lt. 1 ) exit !This is the boundry condition. if this is true it has left the surface.
    uvpath(l,1) = u
    uvpath(l,2) = v

    Vu(1) = Vu(2)
    Vv(1) = Vv(2)
    
enddo


!now we must write all of uvpath into a file or process it here idk your choice
open(15,file='uvpath.out',status='unknown',action='write')

!this loop iterates for all written uvpath 
do j = 1,l-1,1
write(15,*) uvpath(j,1),uvpath(j,2)
enddo
close(15)
deallocate(uvpath)
deallocate(connection)
call XYZ(switch,l-1,"uvpath.out",10,"xyzpath.out",11)  !This actually outputs the file






!please input a function here that adds a small variation to the uv path and inputs into xyz.f95
!This is so that we are able to test weather or not the distance we get really is a local minimum
!and if that is the case we know that the line we are getting is a geodesic

!Add a new section to the input parameters file where there is a new switch to 1 or 0 to calculate
!The distnace of the uvpath.out on uv using the metric.out but it also checks a small variation in
!the path 



!here is where we will be adding that function -maxvar*(4/time**2)*(-l**2+time*l)
open(55, file='uv+path.out',status='unknown',action='write')
open(56,file='uv-path.out',status='unknown',action='write')
open(15,file='uvpath.out',status='old',action='read')

do j = 0, l-2,1
read(15,*) u,v
f = maxvar*(4.0/(l-1)**2)*(-(j)**2+(l-2)*(j)) 
utemp = u
vtemp = v
u = u + f
v = v + f
call UVtoIndex(i,u,v,ut,umin,vmin,hu,hv)
call IndextoUV(i,ut,hu,hv,umin,vmin,u,v)
write(55,*) u,v
u = utemp - f
v = vtemp - f
call UVtoIndex(i,u,v,ut,umin,vmin,hu,hv)
call IndextoUV(i,ut,hu,hv,umin,vmin,u,v)
write(56,*) u,v
enddo

close(15)
close(55)
close(56)

call distance("uvpath.out",10,l-1,T,ut,umin,vmin,hu,hv,f2)
call distance("uv+path.out",11,l-1,T,ut,umin,vmin,hu,hv,f1)
call distance("uv-path.out",12,l-1,T,ut,umin,vmin,hu,hv,f3)

write(*,*) f1,f2,f3



end program integrator

subroutine IndexToUV(i,ut,hu,hv,umin,vmin,u,v)
        implicit none
        real*8, intent(out) :: u,v
        integer :: a,b,i
        integer, intent(in) :: ut
        real*8, intent(in) :: hu,hv,umin,vmin
        
        a = MOD(i-1,ut)
        b = int((i-1)/ut)
        u = a*hu+umin
        v = b*hv+vmin

        return
end subroutine IndexToUV

subroutine UVToIndex (i,u,v,ut,umin,vmin,hu,hv)
        implicit none
        integer, intent(out) :: i
        integer :: a,b
        integer, intent(in) :: ut
        real*8, intent(in) :: u,v,umin,vmin,hu,hv
        
        a = int((u-umin)/hu)
        b = int((v-vmin)/hv)

        i = b*ut+a+1

        return
end subroutine UVToIndex


subroutine distance(filename,il,l,t,ut,umin,vmin,hu,hv,S)
character (len = 30), intent(in) :: filename
character (len = :), allocatable :: infile
real*8, dimension(:,:), allocatable :: metric
real*8, intent(in) :: umin,vmin,hu,hv
integer, intent(in) :: il,l,t,ut
integer :: i,j
real*8 :: u,v,u1,v1,u2,v2
real*8, intent(out) :: S
!Variable defs
!filename is obviouse
!il is the length of the file name
!l is the number of rows in the file it should be l-1
!t is the total number of nodes

allocate(character (len = il) :: infile)
infile = filename(1:il)
allocate(metric(t,6))

open(98,file=infile,status='old',action='read')
open(100,file='metric.out',status='old',action='read')
do i = 1, t, 1
read(100,*) metric(i,1:6) 
enddo
close(100)

S = 0

do j = 1,l-1,1

if(j .eq. 1) then
    read(98,*) u,v
    read(98,*) u1,v1
else
    u = u1
    v = v1
    read(98,*) u1,v1
endif


u2 = u1-u
v2 = v1-v

!calling for for i for the integration
call UVtoIndex(i,u,v,ut,umin,vmin,hu,hv)

!This has been tested with straight line lengths and it works very well :)))))) This is finished lets goo
S = S + sqrt(metric(i,1)*u2*u2+2*metric(i,2)*u2*v2+metric(i,3)*v2*v2)

enddo

!So the way this function will work is that we will input a file and the length of the file name
!Then it use the inputted dt value to 


close(98)
deallocate(metric)
deallocate(infile)
end subroutine distance




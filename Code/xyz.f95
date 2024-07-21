!This is where the code for the
!x=f(u,v)
!y=g(u,v)
!z=k(u,v)
!This code in theory will be ran once at the begining giving
!us a 2d array of size (t,3).
!this gives us the xyz coords of each point i on our u,v mesh
!We will use this array in conbination with the back and forth
!Subroutines defined below main to calculate all needed partial
!derivatves for metric
!
!
!we will then store it in a file (then deallocate the array in memeory)
!to be graphed as a cool 3d shape :)
subroutine XYZ(switch,t,inputfile,il,outputfile,ol)
!these are buffer strings they aren't real
character (len = 30),intent(in) :: inputfile,outputfile
character (len = :),allocatable :: infile,outfile
integer,intent(in) :: switch,t,il,ol
integer :: error,i
real*8 :: u,v,x,y,z

allocate(character (len = il) :: infile)
allocate(character (len = ol) :: outfile)
infile = inputfile(1:il)
outfile = outputfile(1:ol)

open(12,file=outfile,status='replace', action='write',position='append',iostat=error)
open(11,file=infile,status='old', action='read',position='rewind',iostat=error)
if(error .ne. 0) stop
!make this code read the parameters for the constant from the parameters inp file
!if you can make it such that switch is also read from the file and not needed to
!copy it from the main mesh

select case(switch)

case(1)
        !Sphere
        do i = 1,t,1
        read(11,*)u,v
        x=cos(u)*sin(v/2)
        y=sin(u)*sin(v/2)
        z=cos(v/2)
        write(12,*)x,y,z
        end do
case(2)
        !Dounut
        do i = 1,t,1
        read(11,*)u,v
	x= (3.0+1.0*cos(u))*cos(v) 
        y= (3.0+1.0*cos(u))*sin(v)
        z= 1.0*sin(u)
        write(12,*)x,y,z
        end do
case(3)
        !gaberiles horn
        do i = 1,t,1
        read(11,*)u,v
	x = 2.0*sin(v)/(u+0.3)
        y = 2.0*cos(v)/(u+0.3)
        z = -u
        write(12,*)x,y,z
        enddo
case(4)
        do i = 1,t,1
        read(11,*)u,v
	x = u
        y = v
        z = exp(-((u-2)**2+(v-2)**2))
        write(12,*) x,y,z
        enddo
case(5)
        do i = 1,t,1
	read(11,*)u,v
	x = u
	y = v
	z = sin(u)*cos(v)
	write(12,*)x,y,z
	enddo
case default
        write(*,*)"bruh how do you not have a proper switch condition"
        stop
end select

close(11)
close(12)
deallocate(infile)
deallocate(outfile)
return
end subroutine XYZ

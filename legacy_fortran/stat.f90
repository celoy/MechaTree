program main

implicit none
integer 		:: i, j
integer			:: data(1000)=0
character(80)	:: cfile
character(80)	:: lfile

call system('ls ls -1 -p T* > list')
open(15,file="list")
i = 0
do 
	read(15,'(A)',END=10) cfile
	if ( cfile(1:5)=='Trial' ) then
		i = i +1
		j = 0
		lfile = cfile
	end if
	if ( cfile(1:5)=='ZZZ00' ) then
		j = j + 1
		data(i) = j
		if ( j==101 ) print *, lfile
	end if  
! 	print *, cfile, 'Trial', i , 'File', j 
end do 
10 continue
print *, data(1:136)
end program main
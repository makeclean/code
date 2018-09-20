module mero_fortran
  abstract interface 
     ! send array
     subroutine send_array_int(array, array_size, block_size, & 
          block_count, idhi, idlo) bind(C)
       use iso_c_binding
       implicit none
       integer(kind=C_INT) :: array(*)
       integer(kind=C_INT),value :: array_size
       integer(kind=C_INT),value :: block_size
       integer(kind=C_INT)  :: block_count
       integer(kind=C_INT64_T) :: idhi
       integer(kind=C_INT64_T) :: idlo
     end subroutine send_array_int

     ! send array
     subroutine send_array_real(array, array_size, block_size, & 
          block_count, idhi, idlo) bind(C)
       use iso_c_binding
       implicit none
       real(kind=C_FLOAT) :: array(*)
       integer(kind=C_INT),value :: array_size
       integer(kind=C_INT),value :: block_size
       integer(kind=C_INT)  :: block_count
       integer(kind=C_INT64_T) :: idhi
       integer(kind=C_INT64_T) :: idlo
     end subroutine send_array_real

     ! send array
     subroutine send_array_long(array, array_size, block_size, & 
          block_count, idhi, idlo) bind(C)
       use iso_c_binding
       implicit none
       integer(kind=C_LONG) :: array(*)
       integer(kind=C_INT),value :: array_size
       integer(kind=C_INT),value :: block_size
       integer(kind=C_INT)  :: block_count
       integer(kind=C_INT64_T) :: idhi
       integer(kind=C_INT64_T) :: idlo
     end subroutine send_array_long

     ! send array
     subroutine send_array_double(array, array_size, block_size, & 
          block_count, idhi, idlo) bind(C)
       use iso_c_binding
       implicit none
       real(kind=C_DOUBLE) :: array(*)
       integer(kind=C_INT),value :: array_size
       integer(kind=C_INT),value :: block_size
       integer(kind=C_INT)  :: block_count
       integer(kind=C_INT64_T) :: idhi
       integer(kind=C_INT64_T) :: idlo
     end subroutine send_array_double

     subroutine recieve_array(array_recieved, array_length, block_size, &
          block_count, idhi, idlo ) bind(C)
       use iso_c_binding
       implicit none
       type(C_PTR) :: array_recieved 
       integer(kind=C_INT),value :: array_length
       integer(kind=C_INT),value :: block_size
       integer(kind=C_INT),value :: block_count
       integer(kind=C_INT64_T),value :: idhi
       integer(kind=C_INT64_T),value :: idlo
     end subroutine recieve_array
     
     subroutine start_clovis() bind(C)
     end subroutine start_clovis 

     subroutine finish_clovis() bind(C)
     end subroutine finish_clovis
  end interface
end module

program p121_get_ensi
! program to retrieve the arrays from mero and to write
! an ensight gold file
use iso_c_binding
use mero_fortran

implicit none
real(kind=C_DOUBLE), allocatable :: array(:) ! data to fill
type(C_PTR) :: p
real(C_DOUBLE), pointer :: array_recieved(:)
integer(C_INT) :: i ! loop variable
integer(C_INT) :: array_size
integer(C_INT) :: block_size = 4096
integer(C_INT) :: block_count 
integer(C_INT64_T) :: idhi
integer(C_INT64_T) :: idlo

! mero procedures
procedure(start_clovis)  :: mero_start  ! start mero
procedure(finish_clovis) :: mero_finish ! end mero
procedure(recieve_array) :: mero_recieve_array_double ! send array

write(*,*) 'enter array size, blockcount, idhi, idlo'
! read the data
do i = 1,3
  read(*,*) array_size,block_count,idhi,idlo
enddo

stop 

call mero_start()

allocate(array_recieved(array_size))

call mero_recieve_array_double(p,array_size,block_size,block_count,idhi,idlo)
call C_F_POINTER(p,array_recieved,[array_size])

! write the ensight gold mesh
open(12,file="p121_medium.ensi.DISPL-000001",status='replace',       &
     action='write')
write(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
write(12,'(A/A/A)') "part", "     1","coordinates"
do i  = 1, array_size
  write(12,*) array_recieved(i)
enddo
close(12) 

call mero_finish()

deallocate(array_recieved)

end program p121_get_ensi


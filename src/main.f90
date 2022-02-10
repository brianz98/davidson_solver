program main

   use const
   use error_handling, only: check_allocate
   use davidson, only: davidson_solver
   use linalg, only: eigs_wrapper

   implicit none

   real(p), allocatable :: A(:,:), eigval(:), eigval_exact(:)
   integer :: dimA, nEig, nTrial, maxSize, ierr, i, maxiter
   integer, parameter :: iunit = 6
   integer(kind=8) :: t0, t1, count_rate, count_max
   real(p) :: tol

   dimA = 2000
   nEig = 4
   nTrial = nEig*2
   maxSize = 100
   tol = 1e-3
   maxiter = 100

   allocate(A(dimA,dimA),source=0.0_p,stat=ierr)
   call check_allocate('A',dimA**2,ierr)

   allocate(eigval(nEig))
   allocate(eigval_exact(dimA))

   call random_number(A)

   ! Distributed -0.05 - 0.05
   A = (A*2 - 1.0_p)/20

   do i = 1, dimA
      A(i,i) = A(i,i) + i
   end do

   ! Real symmetric (hence also Hermitian) matrix
   A = (transpose(A) + A)/2

   call system_clock(t0, count_rate, count_max)
   call davidson_solver(A, eigval, nTrial, maxSize, tol, maxiter)
   call system_clock(t1)
   write(iunit, '(1X,A,I0,A,I0,A,F14.6)') 'Time taken for Davidson for a (',dimA,',',dimA,') matrix: ',&
      real((t1-t0),kind=p)/real(count_rate, kind=p)
   write(iunit, '(1X,A,I0,A)') 'The lowest ',nEig,' eigenvalues:'
   do i = 1, nEig
      write(iunit, '(1X, F14.6)') eigval(i)
   end do

   call system_clock(t0)
   call eigs_wrapper(A, dimA, eigval_exact, ierr)
   call system_clock(t1)

   write(iunit, '(1X,A,I0,A,I0,A,F14.6)') 'Time taken for dsyev for a (',dimA,',',dimA,') matrix: ',&
      real((t1-t0),kind=p)/real(count_rate, kind=p)
   write(iunit, '(1X,A,I0,A)') 'The lowest ',nEig,' eigenvalues:'
   do i = 1, nEig
      write(iunit, '(1X, F14.6)') eigval_exact(i)
   end do

end program main
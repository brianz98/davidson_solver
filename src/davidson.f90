module davidson
   use const

   implicit none

   contains

   subroutine davidson_solver(A, eigval, nTrial, maxSize, tol)
      ! The Davidson's method
      use linalg, only: qr_wrapper, dgemm_wrapper, eigs_wrapper
      use error_handling, only: check_allocate

      real(p), intent(inout) :: A(:,:)
      real(p), intent(out) :: eigval(:)
      real(p), intent(in) :: tol
      integer, intent(in) :: nTrial, maxsize

      integer :: nEig, dimA, maxEig, ierr
      integer :: i, j
      real(p), allocatable :: V(:,:), theta(:), theta_old(:), tmp(:,:), tmpV(:), eye(:,:), T(:,:), w(:)

      nEig = size(eigval); dimA = size(A, dim=1)
      allocate(theta_old(nEig), source=0.0_p)

      ! Integer division always rounds towards zero, i.e 8*50/8 = 48
      maxEig = nTrial*(maxSize/nTrial)
      allocate(theta(maxEig))
      allocate(V(dimA,maxEig+nTrial),source=0.0_p,stat=ierr)
      call check_allocate('V',dimA*(maxEig+nTrial),ierr)

      allocate(T(maxEig,maxEig),source=0.0_p,stat=ierr)
      call check_allocate('T',maxEig**2,ierr)

      allocate(tmp,source=V,stat=ierr)
      call check_allocate('tmp',dimA*maxEig,ierr)

      allocate(tmpV(dimA),source=0.0_p,stat=ierr)
      call check_allocate('tmpV',dimA,ierr)

      allocate(w(dimA),source=0.0_p,stat=ierr)
      call check_allocate('w',dimA,ierr)

      allocate(eye(dimA,dimA),source=0.0_p,stat=ierr)
      call check_allocate('eye',dimA**2,ierr)
      do i = 1, dimA
         eye(i,i) = 1.0_p
      end do

      ! Initial guesses are nTrial lowest unit vectors
      do i = 1, nTrial
         V(i,i) = 1.0_p
      end do

      do i = nTrial*2, maxSize, nTrial
         ! V is only passed in by reference, and (dimA, i) specifies the slice that should be worked on
         call qr_wrapper(V, dimA, i)

         ! Likewise, tmp is designed to hold the largest tmp matrix but only the relevant slice will be written to
         call dgemm_wrapper('N', 'N', dimA, i, dimA, A, V, tmp)
         call dgemm_wrapper('T', 'N', i, i, dimA, V, tmp, T)

         call eigs_wrapper(T, i, theta, ierr)

         do j = 1, nTrial
            call dgemm_wrapper('N','N',dimA,1,i,V,T,tmpV)
            call dgemm_wrapper('N','N',dimA,1,dimA,(A-theta(j)*eye),tmpV,w)
            w = w/(theta(j)-A(j,j))
            V(:,(i+j)) = w
         end do

         if (sqrt(sum((theta(1:nEig)-theta_old)**2)) < tol) exit
         theta_old = theta(1:nEig)
      end do

      eigval = theta(1:nEig)

   end subroutine

end module davidson
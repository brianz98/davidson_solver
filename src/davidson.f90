module davidson
   use const

   implicit none

   contains

   subroutine davidson_solver(A, eigval, ntrial, maxSize, tol, maxiter)
      ! The Davidson's method
      
#ifdef GPU
      use cublas_v2
#endif

      use linalg, only: qr_wrapper, dgemm_wrapper, eigs_wrapper, dgemv_wrapper
      use error_handling, only: check_allocate

      real(p), intent(inout) :: A(:,:)
      real(p), intent(out) :: eigval(:)
      real(p), intent(in) :: tol
      integer, intent(in) :: ntrial, maxsize, maxiter

      integer :: nEig, dimA, maxguess, ierr, nactive
      integer :: i, j
      integer, parameter :: iunit = 6
      real(p), allocatable :: V(:,:), theta(:), theta_old(:), tmp(:,:), tmpV(:), T(:,:), w(:), V_coll(:,:)
      logical, allocatable :: normconv(:)
      logical :: conv
      real(p) :: norm

#ifdef GPU
      real(p), device, allocatable :: A_d(:,:), V_d(:,:), tmp_d(:,:), tmpV_d(:), w_d(:)
      type(cublasHandle) :: h
      integer(4) :: istat
#endif

      nEig = size(eigval); dimA = size(A, dim=1)
      allocate(theta_old(nEig), source=0.0_p)
      allocate(normconv(ntrial))

      ! Integer division always rounds towards zero, i.e 8*50/8 = 48
      maxguess = ntrial*(maxSize/ntrial)
      allocate(theta(maxguess))
      allocate(V(dimA,maxguess+ntrial),source=0.0_p,stat=ierr)
      call check_allocate('V',dimA*(maxguess+ntrial),ierr)
      allocate(V_coll(dimA,ntrial),source=0.0_p,stat=ierr)

      allocate(T(maxguess,maxguess),source=0.0_p,stat=ierr)
      call check_allocate('T',maxguess**2,ierr)

      allocate(tmp,source=V,stat=ierr)
      call check_allocate('tmp',dimA*maxguess,ierr)

      allocate(tmpV(dimA),source=0.0_p,stat=ierr)
      call check_allocate('tmpV',dimA,ierr)

      allocate(w(dimA),source=0.0_p,stat=ierr)
      call check_allocate('w',dimA,ierr)

#ifdef GPU
      allocate(A_d(dimA, dimA))
      allocate(V_d(dimA,maxguess+ntrial))
      allocate(tmp_d(dimA,maxguess+ntrial))
      allocate(tmpV_d(dimA))
      allocate(w_d(dimA))

      A_d = A
#endif

      ! Initial guesses are ntrial lowest unit vectors
      do i = 1, ntrial
         V(i,i) = 1.0_p
      end do

      nactive = ntrial
      conv = .false.

      do i = 1, maxiter
         ! Orthonormalise the current set of guess vectors
         call qr_wrapper(V, dimA, nactive)

         ! Form the subspace Hamiltonian / Rayleigh matrix
         ! T = V^T A V
#ifdef GPU
         V_d = V
         istat = cublasDgemm(h, 'N', 'N', dimA, nactive, dimA, 1.0_p, A_d, dimA, V_d, dimA, 0.0_p, tmp_d, dimA)
         if (istat .ne. 0) write(iunit,*) 'cublasDgemm error! Returned ',istat
         tmp = tmp_d
#else
         call dgemm_wrapper('N', 'N', dimA, nactive, dimA, A, V, tmp)
#endif
         call dgemm_wrapper('T', 'N', nactive, nactive, dimA, V, tmp, T)

         ! Diagonalise the subspace Hamiltonian
         ! T C = C t (T gets overwritten by C in syev, theta is the diagonal of t)
         call eigs_wrapper(T, nactive, theta, ierr)

         ! We don't use this norm for convergence checking as after each subspace collapse the change in norm is 
         ! essentially zero, but we report it nonetheless as during each restart-block it is still 
         ! a useful measure of convergence.
         norm = sqrt(sum((theta(1:nEig)-theta_old)**2))
         write(iunit,'(1X, A, I3, A, I3, A, ES15.6)') 'Iteration ', i, ', basis size ', nactive, ', rmsE ', norm
         theta_old = theta(1:nEig)

         if (all(normconv)) then
            write(iunit,'(1X, A, ES10.4, A)') 'Residue tolerance of ', tol,' reached, printing results...'
            conv = .true.
            exit
         end if

         if (nactive <= (maxguess-ntrial)) then
            ! If the number of guess vectors can be grown by at least another lot of ntrial
            do j = 1, ntrial
               ! Residue vector w = (A-theta(j)*I) V T(:,j)
               ! Storing a diagonal matrix as large as A is obviously a bad idea, so we use a tmp vector
               ! Technically speaking, tmpV is the 'Ritz vector' and w is the residue vector
               call dgemv_wrapper(dimA, nactive, V, T(:,j), tmpV)
#ifdef GPU
                  istat = cublasDgemv(h, 'N', dimA, dimA, 1.0_p, A_d, dimA, tmpV, 1, 0.0_p, w_d, 1)
                  if (istat .ne. 0) write(iunit,*) 'cublasDgemv error! Returned ',istat
                  w = w_d
#else
               call dgemv_wrapper(dimA, dimA, A, tmpV, w)
#endif
               w = w - theta(j)*tmpV
               if (sqrt(sum(w**2)) < tol) normconv(j) = .true.
               ! Precondition the residue vector to form the correction vector,
               ! if preconditioner = 1, we recover the Lanczos algorithm.
               if (abs((theta(j)-A(j,j))) < depsilon) then
                   w = w/(theta(j) - A(j,j) + 0.01)
               else
                   w = w/(theta(j)-A(j,j))
               end if
               V(:,(nactive+j)) = w
            end do
            nactive = nactive + ntrial
         else
            ! We need to collapse the subspace into ntrial best guesses and restart the iterations
            ! V holds the approximate eigenvectors and T holds the CI coefficients,
            ! so one call to gemm gives us the actual guess vectors.
            write(iunit, '(1X, A)') 'Collapsing subspace...'
            call dgemm_wrapper('N', 'N', dimA, ntrial, maxguess, V, T, V_coll)
            V(:,:ntrial) = V_coll
            nactive = ntrial
         end if
      end do

      if (conv .eqv. .false.) write(iunit, *) 'Davidson did not converge'
      eigval = theta(1:nEig)

#ifdef GPU
      istat = cublasDestroy(h)
      if (istat .ne. 0) write(iunit, *) 'cublasDestroy error! Returned ',istat
#endif

   end subroutine

end module davidson
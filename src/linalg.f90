module linalg
   use const

   implicit none

   contains
      subroutine eigs_wrapper(A, lda, eigval, ierr)
         ! Wrapper around LAPACK dsyev, so named to be reminiscent of the NumPy package 
         ! Finds the eigenvalues and eigenvectors of a double-precision symmetric matrix
         ! This further automates finding the optimal workspace by first passing lwork = -1
         ! to LAPACK dsyev, where work(1) contains the optimised workspace dimension,
         ! then we pass that to dsyev again for actual computation of eigen quantities

         real(p), intent(inout) :: A(*)
         integer, intent(in) :: lda
         real(p), intent(out) :: eigval(:)
         integer, intent(out) :: ierr

         real(p), allocatable :: work(:)
         integer :: lwork, i

         lwork = -1
         do i = 1, 2
            allocate(work(abs(lwork)))
            call dsyev('V', 'U', lda, A, lda, eigval, work, lwork, ierr)
            lwork = nint(work(1))
            deallocate(work)
         end do
      end subroutine eigs_wrapper

      subroutine dgemm_wrapper(transA, transB, outer_row, outer_col, inner_dim, A, B, C)
         ! Wraps around dgemm with fewer arguments

         character(1), intent(in) :: transA, transB
         integer, intent(in) :: outer_row, outer_col, inner_dim
         real(p), intent(in) :: A(*), B(*)
         real(p), intent(inout) :: C(*)

         integer :: LDA, LDB

         if (transA == 'T') then
            LDA = inner_dim
         else
            LDA = outer_row
         end if

         if (transB == 'T') then
            LDB = outer_col
         else
            LDB = inner_dim
         end if
         
         call dgemm(transA, transB, outer_row, outer_col, inner_dim, 1.0_p, &
                    A, LDA, B, LDB, 0.0_p, C, outer_row)
      end subroutine dgemm_wrapper

      subroutine qr_wrapper(A, rowA, colA)
         ! Wrapper around dgeqrf and dorgqr with only the information we need
         ! On exit the matrix A is the orthonormalised matrix Q
         ! Not specifying A as A(:,:) (assumed shape) makes it possible to pass in 
         ! pointers to the start of array slices without having to create array temporaries
         ! And specifying rowA and colA makes it possible to pass in a large matrix A 
         ! but only letting dgeqrf do work on the specified slice.
         real(p), intent(inout) :: A(*)
         integer, intent(in) :: rowA, colA
         real(p), allocatable :: tau(:), work(:)
         integer :: lwork, i, ierr

         allocate(tau(min(rowA,colA)))

         lwork = -1
         do i = 1, 2
            allocate(work(abs(lwork)))
            call dgeqrf(rowA,colA,A,rowA,tau,work,lwork,ierr)
            lwork = nint(work(1))
            deallocate(work)
         end do

         lwork = -1
         do i = 1, 2
            allocate(work(abs(lwork)))
            call dorgqr(rowA,colA,min(rowA,colA),A,rowA,tau,work,lwork,ierr)
            lwork = nint(work(1))
            deallocate(work)
         end do
      end subroutine
      
end module linalg

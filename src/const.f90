module const
   implicit none

   ! Integer types
   integer, parameter :: int_32 = selected_int_kind(6)
   integer, parameter :: int_64 = selected_int_kind(15)

   ! Real types
   integer, parameter :: sp = selected_real_kind(6,37)
   integer, parameter :: dp = selected_real_kind(15,307)

   integer, parameter :: p = dp

end module const
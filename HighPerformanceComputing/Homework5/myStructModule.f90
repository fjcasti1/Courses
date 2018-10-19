module myStructModule
  use precision
  implicit none
  type paramStruct
    integer  :: k
    real(DP) :: a, b, L
  end type paramStruct
end myStructModule

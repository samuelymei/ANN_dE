module atom_m
  use precision_m
  implicit none
  type atom_t
    character(len=2) :: element
    integer(kind=4) :: atomic_number
    real(kind=fp_kind) :: atomic_mass
    character(len=4) :: name
    real(kind=fp_kind) :: crd(3)
    real(kind=fp_kind) :: eta, xi
    real(kind=fp_kind) :: g1, g2
    real(kind=fp_kind) :: dE
    contains
      procedure :: translate
  end type atom_t

  type, extends(atom_t) :: mm_atom_t
    real(kind=fp_kind) :: charge
    real(kind=fp_kind) :: dipole
  end type mm_atom_t

  contains
    subroutine translate(this,r)
      class( atom_t ) :: this
      real(kind=fp_kind), intent(in) :: r(3)
      this%crd = this%crd + r
    end subroutine translate

end module atom_m

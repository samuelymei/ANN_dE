module atom_m
  use precision_m
  use ann_m
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
    type(NeuralNetwork_t), pointer :: ann
    contains
      procedure :: translate
      procedure :: computedE
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

    subroutine computedE(this)
      use ann_m
      class( atom_t ) :: this
      real(kind=fp_kind) :: ann_input(2)
      real(kind=fp_kind), allocatable :: ann_output(:)
      allocate(ann_output(this%ann%nOutput))
      ann_input(1) = this%g1
      ann_input(2) = this%g2
      ann_output = this%ann%propagate(ann_input)
      this%dE = ann_output(1)
      deallocate(ann_output)
    end subroutine

end module atom_m

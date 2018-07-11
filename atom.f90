module atom_m
  use precision_m
  use ann_m
  implicit none
  type atom_t
    character(len=2) :: element=''
    integer(kind=4) :: atomic_number=0
    real(kind=fp_kind) :: atomic_mass=0.d0
    character(len=4) :: name=''
    real(kind=fp_kind) :: crd(3)=0.d0
    real(kind=fp_kind) :: eta=0.d0, xi=0.d0
    real(kind=fp_kind) :: g1=0.d0, g2=0.d0
    real(kind=fp_kind) :: dE=0.d0
    type(NeuralNetwork_t), pointer :: ann => null()
    contains
      procedure :: resetatom
      procedure :: translate
      procedure :: computedE
      procedure :: copyfromatom
  end type atom_t

  type, extends(atom_t) :: mm_atom_t
    real(kind=fp_kind) :: charge
    real(kind=fp_kind) :: dipole
  end type mm_atom_t

  contains
    subroutine resetatom(this)
      class( atom_t ) :: this
      this%element=''
      this%atomic_number=0
      this%atomic_mass=0.d0
      this%name=''
      this%crd=0.d0
      this%eta=0.d0
      this%xi=0.d0
      this%g1=0.d0
      this%g2=0.d0
      this%dE=0.d0
!      if(associated(this%ann))nullify(this%ann)
    end subroutine resetatom

    subroutine copyfromatom(this,that)
      class( atom_t ) :: this
      type(atom_t) :: that
      this%element = that%element
      this%atomic_number = that%atomic_number
      this%atomic_mass = that%atomic_mass
      this%name = that%name
      this%crd = that%crd
      this%eta = that%eta
      this%xi = that%xi
      this%g1 = that%g1
      this%g2 = that%g2
      this%dE = that%dE
      if(associated(that%ann))this%ann=that%ann
    end subroutine copyfromatom

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

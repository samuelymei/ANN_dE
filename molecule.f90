module molecule_m
  use precision_m
  use atom_m
  implicit none
  type molecule_t
    integer(kind=4) :: num_atoms
    integer(kind=4) :: num_elements
    character(len=2), allocatable :: elements(:)
    real(kind=fp_kind) :: molecular_mass
    type(atom_t), allocatable :: atoms_p(:)
    real(kind=fp_kind) :: charge  ! may allow fractional charge
    real(kind=fp_kind) :: dipole(3) 
    integer(kind=4) :: spin_multiplicity
    contains
      procedure :: whole_translate
      procedure :: read_in_from_xyz
      procedure :: count_elements
      procedure :: gfactor
  end type molecule_t
  contains
    subroutine whole_translate(this,r)
      class(molecule_t) :: this
      real(kind=fp_kind) :: r(3)
      integer(kind=4) :: i, j, k
      do i = 1, this%num_atoms
        call this%atoms_p(i)%translate(r)
      end do
    end subroutine whole_translate

    subroutine read_in_from_xyz(this,f_xyz)
      implicit none
      class(molecule_t) :: this
      integer(kind=4) :: f_xyz
      integer(kind=4) :: idx_atom
      read(f_xyz,*) this%num_atoms
      read(f_xyz,*)
      allocate(this%atoms_p(this%num_atoms))
      do idx_atom = 1, this%num_atoms
        read(f_xyz,*)this%atoms_p(idx_atom)%element, this%atoms_p(idx_atom)%crd(1:3)
      end do
    end subroutine read_in_from_xyz

    subroutine count_elements(this)
      implicit none
      class(molecule_t) :: this
      integer(kind=4) :: nelement
      character(len=2), allocatable :: elements(:)
      logical :: is_new_element
      integer(kind=4) :: idx_atom
      integer(kind=4) :: i, j, k
      if(this%num_atoms>0)then
        allocate(elements(this%num_atoms))
        nelement = 1
        elements(nelement) = this%atoms_p(1)%element
        do idx_atom = 2, this%num_atoms
          is_new_element = .true.
          do i = 1, nelement
            if(elements(i) == this%atoms_p(idx_atom)%element)then
              is_new_element = .false.
              exit
            end if
          end do
          if(is_new_element)then
            nelement = nelement + 1
            elements(nelement) = this%atoms_p(idx_atom)%element
          end if
        end do
        this%num_elements = nelement
        allocate(this%elements(this%num_elements))
        this%elements = elements
!        print*,'Number of elements:', this%num_elements
        deallocate(elements)
      end if

    end subroutine count_elements

    subroutine gfactor(this)
      implicit none
      class(molecule_t) :: this
      real(kind=fp_kind) :: rs = 0.0
      real(kind=fp_kind) :: rc = 6.0
      real(kind=fp_kind), allocatable :: r_matrix(:,:)
      real(kind=fp_kind), allocatable :: fc_matrix(:,:)
      real(kind=fp_kind), allocatable :: cos_angle_matrix(:,:,:)
      integer(kind=4) :: idx_atom, jdx_atom, kdx_atom
      do idx_atom = 1, this%num_atoms
        if(this%atoms_p(idx_atom)%element == 'C ')then
          this%atoms_p(idx_atom)%eta = 0.8d0
          this%atoms_p(idx_atom)%xi = 0.8d0
        else if(this%atoms_p(idx_atom)%element == 'H ')then
          this%atoms_p(idx_atom)%eta = 0.2d0
          this%atoms_p(idx_atom)%xi = 0.6d0
        else if(this%atoms_p(idx_atom)%element == 'O ')then
          this%atoms_p(idx_atom)%eta = 0.2d0
          this%atoms_p(idx_atom)%xi = 0.9d0
        else if(this%atoms_p(idx_atom)%element == 'N ')then
          this%atoms_p(idx_atom)%eta = 0.8d0
          this%atoms_p(idx_atom)%xi = 1.0d0
        else if(this%atoms_p(idx_atom)%element == 'Cl')then
          this%atoms_p(idx_atom)%eta = 0.09d0
          this%atoms_p(idx_atom)%xi = 0.09d0
        else
          print*,'Element not in the library |', this%atoms_p(:)%element,'|'
          stop
        end if
      end do
      allocate(r_matrix(this%num_atoms,this%num_atoms))
      allocate(fc_matrix(this%num_atoms,this%num_atoms))
      allocate(cos_angle_matrix(this%num_atoms,this%num_atoms,this%num_atoms))
      do idx_atom = 1, this%num_atoms
        do jdx_atom = 1, this%num_atoms
          if(jdx_atom == idx_atom) cycle
          r_matrix(idx_atom,jdx_atom) = distance2(this%atoms_p(idx_atom)%crd,this%atoms_p(jdx_atom)%crd)
          fc_matrix(idx_atom,jdx_atom) = fc(r_matrix(idx_atom,jdx_atom), rc)
          do kdx_atom = 1, this%num_atoms
            if(kdx_atom == jdx_atom .or. kdx_atom == idx_atom) cycle
            cos_angle_matrix(idx_atom,jdx_atom,kdx_atom) = cos_angle(this%atoms_p(idx_atom)%crd,this%atoms_p(jdx_atom)%crd,this%atoms_p(kdx_atom)%crd)
          end do
        end do
      end do
      this%atoms_p(:)%g1=0.d0
      this%atoms_p(:)%g2=0.d0
      do idx_atom = 1, this%num_atoms
        do jdx_atom = idx_atom + 1, this%num_atoms
          if(jdx_atom == idx_atom) cycle
          this%atoms_p(idx_atom)%g1 = this%atoms_p(idx_atom)%g1 + exp(-this%atoms_p(idx_atom)%eta*(r_matrix(idx_atom,jdx_atom)-rs)**2)*fc_matrix(idx_atom,jdx_atom)
          this%atoms_p(jdx_atom)%g1 = this%atoms_p(jdx_atom)%g1 + exp(-this%atoms_p(jdx_atom)%eta*(r_matrix(idx_atom,jdx_atom)-rs)**2)*fc_matrix(idx_atom,jdx_atom)
        end do
      end do
      do idx_atom = 1, this%num_atoms
        do jdx_atom = 1, this%num_atoms
          if(jdx_atom == idx_atom) cycle
          do kdx_atom = 1, this%num_atoms
            if(kdx_atom == jdx_atom .or. kdx_atom == idx_atom) cycle
            this%atoms_p(idx_atom)%g2 = this%atoms_p(idx_atom)%g2 + &
                                      & 2.0d0**(1-this%atoms_p(idx_atom)%xi)*(1+cos_angle_matrix(idx_atom,jdx_atom,kdx_atom))**this%atoms_p(idx_atom)%xi &
                                      & *exp(-this%atoms_p(idx_atom)%eta*(r_matrix(idx_atom,jdx_atom)**2+r_matrix(jdx_atom,kdx_atom)**2+r_matrix(idx_atom,kdx_atom)**2)) &
                                      & *fc_matrix(idx_atom,jdx_atom)*fc_matrix(jdx_atom,kdx_atom)*fc_matrix(idx_atom,kdx_atom)
          end do
        end do
      end do
      
      deallocate(r_matrix)
      deallocate(fc_matrix)
      deallocate(cos_angle_matrix)
      contains
        function fc(rij, rc)
          implicit none
          real(kind=fp_kind) :: fc
          real(kind=fp_kind) :: rij, rc
          real(kind=fp_kind) :: pi
          pi = atan(1.0d0)*4.d0
          fc = 0.d0
          if (rij < rc)then
            fc = 0.5*(cos(pi*rij/rc)+1)
          end if
        end function fc

        function distance2(crd1,crd2)
          implicit none
          real(kind=fp_kind) :: distance2
          real(kind=fp_kind), intent(in) :: crd1(3), crd2(3)
          integer(kind=4) :: i
          distance2=sqrt((crd1(1)-crd2(1))**2+(crd1(2)-crd2(2))**2+(crd1(3)-crd2(3))**2)
        end function distance2

        function cos_angle(crd1,crd2,crd3)
          implicit none
          real(kind=fp_kind) :: cos_angle
          real(kind=fp_kind), intent(in) :: crd1(3), crd2(3), crd3(3)
          real(kind=fp_kind) :: r12, r23, r13
          r12 = distance2(crd1,crd2)
          r23 = distance2(crd2,crd3)
          r13 = distance2(crd1,crd3)
          cos_angle = (r12**2+r23**2-r23**2)/(2*r12*r23)
        end function cos_angle
    end subroutine gfactor
end module molecule_m

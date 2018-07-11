program ANN_dE
  use precision_m
  use molecule_m
  use ann_m
  use random_m
  implicit none

  integer(kind=4) :: n_mol_training, n_mol_test
  type(molecule_t), allocatable :: molecules(:)
  type(molecule_t), allocatable :: mol_training(:), mol_test(:)
  integer(kind=4) :: idx_molecule, jdx_molecule, idx_atom
  integer(kind=4), allocatable :: mol_seq(:)
  integer(kind=4) :: nchsn, ichsn
  logical, allocatable :: chosen(:)
  integer(kind=4) :: id_f_xyz
  character(len=20) :: f_xyz

  integer(kind=4) :: nann
  type(NeuralNetwork_t), allocatable, target :: ann(:)
  character(len=2), allocatable :: ann_elements(:)
  integer(kind=4) :: nInput, nHiddenLayer, nOutput ! inputs for the initiation of the network
  integer(kind=4), allocatable :: nNeurons(:) ! inputs for the initiation of the network
  integer(kind=4) :: idx_ann, idx_layer
  integer(kind=4), allocatable :: idx_ann_for_atom(:)
  character(len=20) :: actvfunc

! variables for training
  integer(kind=4) :: irest = 0
  integer(kind=4) :: nvariables, ivariable, nv_tmp
  real(kind=fp_kind), allocatable :: variables(:)
  integer(kind=4) :: iter
  real(kind=fp_kind) :: fret
 
! read in molecules
  n_mol_training = 600
  n_mol_test = 200
  allocate(molecules(n_mol_training+n_mol_test))
  allocate(mol_training(n_mol_training))
  allocate(mol_test(n_mol_test))
  do idx_molecule = 1, n_mol_training + n_mol_test
     write(f_xyz,'(I3)') idx_molecule
     f_xyz='Mgas'//trim(adjustl(f_xyz))//'.xyz'
     id_f_xyz = 11   
     open(id_f_xyz, file = trim(f_xyz), form='formatted')
!     call molecules(idx_molecule)%reset
     call molecules(idx_molecule)%read_in_from_xyz(id_f_xyz)
     close(id_f_xyz)
     call molecules(idx_molecule)%count_elements()
     call molecules(idx_molecule)%gfactor()
  end do

! read in target data
  open(12, file = "../e-diff-shif.dat")
  do idx_molecule = 1, n_mol_training + n_mol_test
    read(12,*) jdx_molecule, molecules(idx_molecule)%target_dE
  end do
  close(12)
  
!shuffle the molecules
  allocate(mol_seq(n_mol_training + n_mol_test))
  allocate(chosen(n_mol_training + n_mol_test))
  chosen=.false.
  nchsn = 0
  do while(nchsn<n_mol_training)
    ichsn = 1 + int(MyUniformRand()*(n_mol_training + n_mol_test))
    if(chosen(ichsn))cycle
    chosen(ichsn) = .true.
    nchsn = nchsn + 1
    mol_seq(nchsn) = ichsn
  end do
  do idx_molecule = 1, n_mol_training + n_mol_test
    if(chosen(idx_molecule))cycle
    nchsn = nchsn + 1
    mol_seq(nchsn) = idx_molecule
  end do

! copy the molecules to training set and test set
  do idx_molecule = 1, n_mol_training
    call mol_training(idx_molecule)%copyfrom(molecules(mol_seq(idx_molecule)))
  end do
  do idx_molecule = 1, n_mol_test
    call mol_test(idx_molecule)%copyfrom(molecules(mol_seq(idx_molecule+n_mol_training)))
  end do

  print*,'Is it a restart job? 0. No. 1. Yes.'
  read*, irest
  if(irest == 1) print*, 'restarting'

! initial ANN
! TODO: build an input section for these control variables
  nInput = 2
  nHiddenLayer = 3
  allocate(nNeurons(nHiddenLayer))
  nNeurons = (/15,30,15/)
  nOutput = 1  ! this program so far can deal with only nOutput = 1
  actvfunc = 'sigmoid'

  nann = mol_training(1)%num_elements ! one neural network for each element
  print*,'Number of artificial Neural Networks', nann
  allocate(ann(nann))
  allocate(ann_elements(nann))
  ann_elements = mol_training(1)%elements

  do idx_ann = 1, nann
    ann(idx_ann) = constructor(nInput,nHiddenLayer,nNeurons,nOutput,actvfunc)
    ann(idx_ann)%ctag = ann_elements(idx_ann)
  end do

! assign each atom in each molecule an ANN
  allocate(idx_ann_for_atom(mol_training(1)%num_atoms))
  do idx_atom = 1, mol_training(1)%num_atoms
    idx_ann_for_atom(idx_atom) = -1
    do idx_ann = 1, nann
      if(mol_training(1)%atoms_p(idx_atom)%element == ann_elements(idx_ann))then
        idx_ann_for_atom(idx_atom) = idx_ann
        exit
      end if
    end do
    print*,'ANN for atom ', idx_atom, ' is ', idx_ann_for_atom(idx_atom), 'with element ', ann(idx_ann_for_atom(idx_atom))%ctag
    if(idx_ann_for_atom(idx_atom) <= 0)then
      print*,'Something weird just happened. Quit!'
      stop
    end if
  end do 

! assign network id to each atom in each molecule
  do idx_molecule = 1, n_mol_training
    do idx_atom = 1, mol_training(idx_molecule)%num_atoms
      mol_training(idx_molecule)%atoms_p(idx_atom)%ann => ann(idx_ann_for_atom(idx_atom))
    end do
  end do
  do idx_molecule = 1, n_mol_test
    do idx_atom = 1, mol_test(idx_molecule)%num_atoms
      mol_test(idx_molecule)%atoms_p(idx_atom)%ann => ann(idx_ann_for_atom(idx_atom))
    end do
  end do

  nvariables = 0
  do idx_ann = 1, nann
    do idx_layer = 0, ann(idx_ann)%nHiddenLayer
      nvariables = nvariables + ann(idx_ann)%nNeurons(idx_layer)*ann(idx_ann)%nNeurons(idx_layer+1)
      nvariables = nvariables + ann(idx_ann)%nNeurons(idx_layer+1)
    end do
  end do
  print*,'There are altogether ',nvariables, ' variables in the networks'

  allocate(variables(nvariables))
  if(irest == 1)then
    print*,'reading restart file'
    read(18)nv_tmp
    if(nv_tmp .ne. nvariables)then
      print*,'this is not a correct restart file. Quiting...'
      stop
    end if
    read(18)variables(1:nvariables)
  else
    do ivariable = 1, nvariables
      variables(ivariable) = (MyUniformRand() - 0.5)*5.0
    end do
  end if

  call dfpmin(variables,nvariables,1.D-6,iter,fret,ann,nann,mol_training,n_mol_training,mol_test,n_mol_test,actvfunc)

!  call gdmin(variables,nvariables,1.D-6,iter,ann,nann,mol_training,n_mol_training,mol_test,n_mol_test)

  deallocate(molecules)
  deallocate(mol_training)
  deallocate(mol_test)
  deallocate(mol_seq)
  deallocate(chosen)
  deallocate(nNeurons)
  deallocate(ann)
  deallocate(ann_elements)
  deallocate(idx_ann_for_atom)
  deallocate(variables)
end program ANN_dE


subroutine loss_func(variables,nvariables,ann,nann,mol_training,n_mol_training,loss_training,mol_test,n_mol_test,loss_test)
  use precision_m
  use molecule_m
  use ann_m
  implicit none
  real(kind=fp_kind), intent(in) :: variables(nvariables)
  integer(kind=4), intent(in) :: nvariables

  type(NeuralNetwork_t), intent(in out) :: ann(nann)
  integer(kind=4), intent(in) :: nann

  type(molecule_t), intent(in out) :: mol_training(n_mol_training)
  integer(kind=4), intent(in) :: n_mol_training
  real(kind=fp_kind), intent(out) :: loss_training

  type(molecule_t), intent(in out) :: mol_test(n_mol_test)
  integer(kind=4), intent(in) :: n_mol_test
  real(kind=fp_kind), intent(out) :: loss_test

  integer(kind=4) :: idx_molecule, idx_atom, idx_ann, idx_layer
  integer(kind=4) :: jdx_neuron, kdx_neuron
  integer(kind=4) :: ivariable

! setup the networks 
  ivariable = 0
  do idx_ann = 1, nann
    do idx_layer = 0, ann(idx_ann)%nHiddenLayer
      do jdx_neuron = 1, ann(idx_ann)%nNeurons(idx_layer)
        do kdx_neuron = 1, ann(idx_ann)%nNeurons(idx_layer+1)
          ivariable = ivariable + 1
          ann(idx_ann)%w(jdx_neuron,kdx_neuron,idx_layer) = variables(ivariable)
        end do
      end do
    end do
    do idx_layer = 1, ann(idx_ann)%nHiddenLayer+1
      do jdx_neuron = 1, ann(idx_ann)%nNeurons(idx_layer)
        ivariable = ivariable + 1
        ann(idx_ann)%bias(jdx_neuron,idx_layer) = variables(ivariable)
      end do
    end do
  end do

  loss_training = 0.d0
  do idx_molecule = 1, n_mol_training
    do idx_atom = 1, mol_training(idx_molecule)%num_atoms
      call mol_training(idx_molecule)%atoms_p(idx_atom)%computedE
    end do
    call mol_training(idx_molecule)%computeMoledE
    loss_training = loss_training + (mol_training(idx_molecule)%dE - mol_training(idx_molecule)%target_dE)**2
  end do
  loss_training = loss_training / n_mol_training

  loss_test = 0.d0
  do idx_molecule = 1, n_mol_test
    do idx_atom = 1, mol_test(idx_molecule)%num_atoms
      call mol_test(idx_molecule)%atoms_p(idx_atom)%computedE
    end do
    call mol_test(idx_molecule)%computeMoledE
    loss_test = loss_test + (mol_test(idx_molecule)%dE - mol_test(idx_molecule)%target_dE)**2
  end do
  loss_test = loss_test / n_mol_test

end subroutine loss_func

subroutine dloss_func(variables,nvariables,ann,nann,mol_training,n_mol_training,loss_training,mol_test,n_mol_test,loss_test,actvfunc,dfdv)
  use precision_m
  use molecule_m
  use ann_m
  implicit none
  real(kind=fp_kind), intent(in) :: variables(nvariables)
  integer(kind=4), intent(in) :: nvariables

  type(NeuralNetwork_t), intent(in out) :: ann(nann)
  integer(kind=4), intent(in) :: nann

  type(molecule_t), intent(in out) :: mol_training(n_mol_training)
  integer(kind=4), intent(in) :: n_mol_training
  real(kind=fp_kind), intent(out) :: loss_training

  type(molecule_t), intent(in out) :: mol_test(n_mol_test)
  integer(kind=4), intent(in) :: n_mol_test
  real(kind=fp_kind), intent(out) :: loss_test

  character(len=*), intent(in) :: actvfunc

  real(kind=fp_kind), intent(out) :: dfdv(nvariables)

  real(kind=fp_kind), allocatable :: dLoss_dOutput(:,:)
  real(kind=fp_kind) :: ann_input(2), ann_output(ann(1)%nOutput)
  integer(kind=4) :: idx_ann, idx_atom, idx_molecule
  integer(kind=4) :: ivariable
  
  type(NeuralNetwork_t) :: working_ann

  integer(kind=4) :: idx_layer, jdx_neuron, kdx_neuron

  call loss_func(variables,nvariables,ann,nann,mol_training,n_mol_training,loss_training,mol_test,n_mol_test,loss_test)

!! clear the gradient matrix
  do idx_ann = 1, nann
    ann(idx_ann)%wg=0.d0
    ann(idx_ann)%zg=0.d0
    ann(idx_ann)%ag=0.d0
    ann(idx_ann)%biasg=0.d0
  end do

  allocate(dLoss_dOutput(mol_training(1)%num_atoms,n_mol_training))
  do idx_molecule = 1, n_mol_training
    dLoss_dOutput(:,idx_molecule) = &
            & 2.d0*(mol_training(idx_molecule)%dE-mol_training(idx_molecule)%target_dE)/n_mol_training
  end do

! it seems that I have to repropagate the networks for all the elements in all
! the snapshots, in order to save memory occupation
  
  working_ann = constructor(ann(1)%nInput,ann(1)%nHiddenLayer,ann(1)%nNeurons(1:),ann(1)%nOutput,actvfunc)
  do idx_molecule = 1, n_mol_training
    do idx_atom = 1, mol_training(idx_molecule)%num_atoms
      working_ann%wg = 0.d0
      working_ann%biasg = 0.d0
      working_ann%ag = 0.d0
      working_ann%zg = 0.d0
      working_ann%z = 0.d0
      working_ann%a = 0.d0
      working_ann%w = mol_training(idx_molecule)%atoms_p(idx_atom)%ann%w
      working_ann%bias = mol_training(idx_molecule)%atoms_p(idx_atom)%ann%bias

      ann_input(1) = mol_training(idx_molecule)%atoms_p(idx_atom)%g1
      ann_input(2) = mol_training(idx_molecule)%atoms_p(idx_atom)%g2
      ann_output = working_ann%propagate(ann_input)
      call working_ann%backpropagate(dLoss_dOutput(idx_atom,idx_molecule))
      mol_training(idx_molecule)%atoms_p(idx_atom)%ann%wg = &
         mol_training(idx_molecule)%atoms_p(idx_atom)%ann%wg + working_ann%wg
      mol_training(idx_molecule)%atoms_p(idx_atom)%ann%biasg = &
         mol_training(idx_molecule)%atoms_p(idx_atom)%ann%biasg + working_ann%biasg
    end do
  end do
  call working_ann%destroy()

  ivariable = 0
  do idx_ann = 1, nann
    do idx_layer = 0, ann(idx_ann)%nHiddenLayer
      do jdx_neuron = 1, ann(idx_ann)%nNeurons(idx_layer)
        do kdx_neuron = 1, ann(idx_ann)%nNeurons(idx_layer+1)
          ivariable = ivariable + 1
          dfdv(ivariable) = ann(idx_ann)%wg(jdx_neuron,kdx_neuron,idx_layer)
        end do
      end do
    end do
    do idx_layer = 1, ann(idx_ann)%nHiddenLayer+1
      do jdx_neuron = 1, ann(idx_ann)%nNeurons(idx_layer)
        ivariable = ivariable + 1
        dfdv(ivariable) = ann(idx_ann)%biasg(jdx_neuron,idx_layer)
      end do
    end do
  end do
  deallocate(dLoss_dOutput)
end subroutine dloss_func

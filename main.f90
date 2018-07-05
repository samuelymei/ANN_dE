program ANN_dE
  use precision_m
  use molecule_m
  use ann_m
  use random_m
  implicit none

  integer(kind=4) :: n_mol_training, n_mol_test
  type(molecule_t), allocatable :: mol_training(:), mol_test(:)
  real(kind=fp_kind), allocatable :: target_training(:), target_test(:)
  real(kind=fp_kind), allocatable :: predicted_training(:), predicted_test(:)
  real(kind=fp_kind) :: training_min, training_max
  real(kind=fp_kind) :: loss
  integer(kind=4) :: idx_molecule, jdx_molecule, idx_atom
  integer(kind=4) :: id_f_xyz
  character(len=20) :: f_xyz

  integer(kind=4) :: nann
  type(NeuralNetwork_t), allocatable :: ann(:)
  character(len=2), allocatable :: ann_elements(:)
  integer(kind=4) :: nInput, nHiddenLayer, nOutput ! inputs for the initiation of the network
  integer(kind=4), allocatable :: nNeurons(:) ! inputs for the initiation of the network
  integer(kind=4) :: idx_ann, idx_layer
  integer(kind=4), allocatable :: idx_ann_for_atom(:)

! variables for training
  integer(kind=4) :: irest = 0
  integer(kind=4) :: nvariables, ivariable, nv_tmp
  integer(kind=4) :: jdx_neuron, kdx_neuron
  real(kind=fp_kind), allocatable :: variables(:)
  integer(kind=4) :: iter
  real(kind=fp_kind) :: fret
 
! read in molecules
  n_mol_training = 600
  n_mol_test = 200
  allocate(mol_training(n_mol_training))
  allocate(mol_test(n_mol_test))
  do idx_molecule = 1, n_mol_training
     write(f_xyz,'(I3)') idx_molecule
     f_xyz='Mgas'//trim(adjustl(f_xyz))//'.xyz'
     id_f_xyz = 11   
     open(id_f_xyz, file = trim(f_xyz), form='formatted')
     call mol_training(idx_molecule)%read_in_from_xyz(id_f_xyz)
     close(id_f_xyz)
     call mol_training(idx_molecule)%count_elements()
     call mol_training(idx_molecule)%gfactor()
  end do

  do idx_molecule = n_mol_training + 1, n_mol_training + n_mol_test
     write(f_xyz,'(I3)') idx_molecule
     f_xyz='Mgas'//trim(adjustl(f_xyz))//'.xyz'
     id_f_xyz = 11
     open(id_f_xyz, file = trim(f_xyz), form='formatted')
     call mol_test(idx_molecule)%read_in_from_xyz(id_f_xyz)
     close(id_f_xyz)
     call mol_test(idx_molecule)%count_elements()
     call mol_test(idx_molecule)%gfactor()
  end do


! read in target data
  allocate(target_training(n_mol_training))
  allocate(target_test(n_mol_test))
  allocate(predicted_training(n_mol_training))
  allocate(predicted_test(n_mol_test))
  open(12, file = "../e-diff-shif.dat")
  do idx_molecule = 1, n_mol_training
    read(12,*) jdx_molecule, target_training(idx_molecule)
  end do
  do idx_molecule = 1, n_mol_test
    read(12,*) jdx_molecule, target_test(idx_molecule)
  end do
  close(12)

! initial ANN
  nInput = 2
  nHiddenLayer = 2
  allocate(nNeurons(nHiddenLayer))
  nNeurons = (/60,30/)
  nOutput = 1  ! this program so far can deal with only nOutput = 1

  nann = mol_training(1)%num_elements ! one neural network for each element
  print*,'Number of artificial Neural Networks', nann
  allocate(ann(nann))
  allocate(ann_elements(nann))
  ann_elements = mol_training(1)%elements

  do idx_ann = 1, nann
    ann(idx_ann) = constructor(nInput,nHiddenLayer,nNeurons,nOutput,"tanh")
    ann(idx_ann)%ctag = ann_elements(idx_ann)
  end do

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
    read(18)nv_tmp
    if(nv_tmp .ne. nvariables)then
      print*,'this is not a correct restart file. Quiting...'
      stop
    end if
    read(18)variables
  else
    do ivariable = 1, nvariables
      variables(ivariable) = (MyUniformRand() - 0.5)*1.0
    end do
  end if

  call dfpmin(variables,nvariables,1.D-6,iter,fret,mol_training,n_mol_training,ann,nann,idx_ann_for_atom,target_training,predicted_training,mol_test,n_mol_test,target_test,predicted_test)
  write(18)nvariables
  write(18,*)variables
  write(50,'(8F10.6)')predicted_training

  deallocate(mol_training)
  deallocate(mol_test)
  deallocate(target_training)
  deallocate(target_test)
  deallocate(predicted_training)
  deallocate(predicted_test)
  deallocate(nNeurons)
  deallocate(ann)
  deallocate(ann_elements)
  deallocate(idx_ann_for_atom)
  deallocate(variables)
end program ANN_dE


subroutine loss_func(variables,nvariables,mol_training,n_mol_training,target_training,ann,nann,idx_ann_for_atom,predicted_training,loss)
  use precision_m
  use molecule_m
  use ann_m
  implicit none
  real(kind=fp_kind), intent(in) :: variables(nvariables)
  integer(kind=4), intent(in) :: nvariables
  type(molecule_t), intent(in) :: mol_training(n_mol_training)
  integer(kind=4), intent(in) :: n_mol_training
  real(kind=fp_kind), intent(in) :: target_training(n_mol_training)
  type(NeuralNetwork_t), intent(inout) :: ann(nann)
  integer(kind=4), intent(in) :: nann
  integer(kind=4), intent(in) :: idx_ann_for_atom(mol_training(1)%num_atoms)
  real(kind=fp_kind), intent(out) :: predicted_training(n_mol_training)
  real(kind=fp_kind), intent(out) :: loss
  real(kind=fp_kind) :: mse
  real(kind=fp_kind) :: ann_input(2), ann_output(ann(1)%nOutput)
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
!          write(22,*)variables(ivariable)
        end do
      end do
    end do
    do idx_layer = 1, ann(idx_ann)%nHiddenLayer+1
      do jdx_neuron = 1, ann(idx_ann)%nNeurons(idx_layer)
        ivariable = ivariable + 1
        ann(idx_ann)%bias(jdx_neuron,idx_layer) = variables(ivariable)
!        write(22,*)variables(ivariable)
      end do
    end do
  end do
!  write(22,*)'-----------------------------------'
!  print*,ivariable, ' variables have been assigned to the networks'
  loss = 0.d0
  mse = 0.d0
  predicted_training = 0.d0
  do idx_molecule = 1, n_mol_training
    do idx_atom = 1, mol_training(idx_molecule)%num_atoms
      idx_ann = idx_ann_for_atom(idx_atom)
      ann_input(1) = mol_training(idx_molecule)%atoms_p(idx_atom)%g1
      ann_input(2) = mol_training(idx_molecule)%atoms_p(idx_atom)%g2
      ann_output = ann(idx_ann)%propagate(ann_input)
      predicted_training(idx_molecule) = predicted_training(idx_molecule) + ann_output(1)    !pay attention here
    end do
!    write(23,*) predicted_training(idx_molecule)
    loss = loss + (predicted_training(idx_molecule) - target_training(idx_molecule))**2
    mse = mse + predicted_training(idx_molecule) - target_training(idx_molecule)
  end do
!  write(23,*)'-------------'
  loss = loss / n_mol_training
  mse = mse / n_mol_training
  write(*,'(A,G12.5,1X,G12.5)') 'Loss function is', loss, mse
end subroutine loss_func

subroutine dloss_func(variables,nvariables,mol_training,n_mol_training,target_training,ann,nann,idx_ann_for_atom,predicted_training,dfdv)
  use precision_m
  use molecule_m
  use ann_m
  implicit none
  real(kind=fp_kind), intent(in) :: variables(nvariables)
  integer(kind=4), intent(in) :: nvariables
  type(molecule_t) :: mol_training(n_mol_training)
  integer(kind=4), intent(in) :: n_mol_training
  real(kind=fp_kind), intent(in) :: target_training(n_mol_training)
  type(NeuralNetwork_t), intent(inout) :: ann(nann)
  integer(kind=4), intent(in) :: nann
  integer(kind=4), intent(in) :: idx_ann_for_atom(mol_training(1)%num_atoms)
  real(kind=fp_kind), intent(in) :: predicted_training(n_mol_training)
  real(kind=fp_kind), intent(out) :: dfdv(nvariables)
  real(kind=fp_kind), allocatable :: dLoss_dOutput(:,:)
  real(kind=fp_kind) :: ann_input(2), ann_output(ann(1)%nOutput)
  integer(kind=4) :: idx_ann, idx_atom, idx_molecule
  integer(kind=4) :: ivariable
  real(kind=fp_kind) :: dfdv_norm
  real(kind=fp_kind) :: loss
  
  type(NeuralNetwork_t) :: working_ann

  integer(kind=4) :: idx_layer, jdx_neuron, kdx_neuron

  call loss_func(variables,nvariables,mol_training,n_mol_training,target_training,ann,nann,idx_ann_for_atom,predicted_training,loss)

!! clear the gradient matrix
  print*,'clear the gradient matrix'
  do idx_ann = 1, nann
    ann(idx_ann)%wg=0.d0
    ann(idx_ann)%zg=0.d0
    ann(idx_ann)%ag=0.d0
    ann(idx_ann)%biasg=0.d0
  end do

  allocate(dLoss_dOutput(mol_training(1)%num_atoms,n_mol_training))
  do idx_molecule = 1, n_mol_training
    dLoss_dOutput(:,idx_molecule) = &
            & 2.d0*(predicted_training(idx_molecule)-target_training(idx_molecule))/n_mol_training
  end do

! it seems that I have to repropagate the networks for all the elements in all
! the snapshots, in order to save memory occupation
  
  working_ann = constructor(ann(1)%nInput,ann(1)%nHiddenLayer,ann(1)%nNeurons(1:),ann(1)%nOutput,"tanh")
  do idx_molecule = 1, n_mol_training
    do idx_atom = 1, mol_training(idx_molecule)%num_atoms
      idx_ann = idx_ann_for_atom(idx_atom)
      working_ann%wg = 0.d0
      working_ann%biasg = 0.d0
      working_ann%ag = 0.d0
      working_ann%zg = 0.d0
      working_ann%z = 0.d0
      working_ann%a = 0.d0
      working_ann%w = ann(idx_ann)%w
      working_ann%bias = ann(idx_ann)%bias
      ann_input(1) = mol_training(idx_molecule)%atoms_p(idx_atom)%g1
      ann_input(2) = mol_training(idx_molecule)%atoms_p(idx_atom)%g2
      ann_output = working_ann%propagate(ann_input)
      call working_ann%backpropagate(dLoss_dOutput(idx_atom,idx_molecule))
      ann(idx_ann)%wg = ann(idx_ann)%wg + working_ann%wg
      ann(idx_ann)%biasg = ann(idx_ann)%biasg + working_ann%biasg
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
!          write(24,'(A,I3,A,I3,A,I4,I4,f12.6)')'wg ann ', idx_ann, ' layer ',idx_layer, ' Neuron', jdx_neuron, kdx_neuron, dfdv(ivariable)
        end do
      end do
    end do
    do idx_layer = 1, ann(idx_ann)%nHiddenLayer+1
      do jdx_neuron = 1, ann(idx_ann)%nNeurons(idx_layer)
        ivariable = ivariable + 1
        dfdv(ivariable) = ann(idx_ann)%biasg(jdx_neuron,idx_layer)
!        write(25,'(A,I3,A,I4,A,I4,f12.6)')'bias ann ', idx_ann, ' layer ',idx_layer, ' Neuron', jdx_neuron, dfdv(ivariable)
      end do
    end do
  end do
  dfdv_norm = dot_product(dfdv,dfdv)
  print*,'Norm of dfdv:', dfdv_norm
!  dfdv=dfdv/sqrt(dfdv_norm)
  deallocate(dLoss_dOutput)
end subroutine dloss_func

module ann_m
  use precision_m
  implicit none
  private
  public :: NeuralNetwork_t, constructor, randomize
  type NeuralNetwork_t
    integer(kind=4) :: nInput
    integer(kind=4) :: nHiddenLayer
    integer(kind=4) :: nOutput
    character(len=2) :: ctag
    integer(kind=4), allocatable :: nNeurons(:)
    real(kind=fp_kind), allocatable :: output(:)
    real(kind=fp_kind), allocatable :: w(:,:,:)
    real(kind=fp_kind), allocatable :: a(:,:)
    real(kind=fp_kind), allocatable :: z(:,:)
    real(kind=fp_kind), allocatable :: bias(:,:)
    real(kind=fp_kind), allocatable :: wg(:,:,:)
    real(kind=fp_kind), allocatable :: ag(:,:)
    real(kind=fp_kind), allocatable :: zg(:,:)
    real(kind=fp_kind), allocatable :: biasg(:,:)
    procedure(tanh), nopass, pointer:: ActivateFunc => NULL()
    procedure(tanhgradient), nopass, pointer:: DeactivateFunc => NULL()
    contains
      procedure :: initann
      procedure :: destroy
      procedure :: randomize
      procedure :: readinput
      procedure :: propagate
      procedure :: backpropagate
  end type NeuralNetwork_t
  interface NeuralNetwork_t
      procedure :: constructor

  end interface
  interface
      function sigmoid(x)
        use precision_m
        implicit none
        real(kind=fp_kind), intent(in) :: x
        real(kind=fp_kind) :: sigmoid
      end function sigmoid

      function sigmoidgradient(x)
        use precision_m
        implicit none
        real(kind=fp_kind), intent(in) :: x
        real(kind=fp_kind) :: sigmoidgradient
        real(kind=fp_kind) :: sigmoid
      end function sigmoidgradient
      
      function tanh(x)
        use precision_m
        implicit none
        real(kind=fp_kind), intent(in) :: x
        real(kind=fp_kind) :: tanh
      end function tanh
      
      function tanhgradient(x)
        use precision_m
        implicit none
        real(kind=fp_kind), intent(in) :: x
        real(kind=fp_kind) :: tanhgradient
      end function tanhgradient

  end interface

  contains
    function constructor(nInput,nHiddenLayer,nNeurons,nOutput,actvfunc)
      implicit none
      type ( NeuralNetwork_t ) :: constructor
      integer(kind=4), intent(in) :: nInput, nHiddenLayer
      integer(kind=4), intent(in) :: nNeurons(nHiddenLayer)
      integer(kind=4), intent(in) :: nOutput
      character(len=*), intent(in) :: actvfunc
      call constructor%initann(nInput,nHiddenLayer,nNeurons,nOutput,actvfunc)
    end function constructor

    subroutine initann(this,nInput,nHiddenLayer,nNeurons,nOutput,actvfunc)
      implicit none
      class(NeuralNetwork_t) :: this
      integer(kind=4), intent(in) :: nInput, nHiddenLayer
      integer(kind=4), intent(in) :: nNeurons(nHiddenLayer)
      integer(kind=4), intent(in) :: nOutput
      character(len=*), intent(in) :: actvfunc
      integer(kind=4) :: maxNeuronsInEachLayer
      this%nInput=nInput
      this%nOutput=nOutput
      this%nHiddenLayer=nHiddenLayer
      allocate(this%nNeurons(0:nHiddenLayer+1))
      this%nNeurons(0)=this%nInput
      this%nNeurons(1:nHiddenLayer)=nNeurons
      this%nNeurons(nHiddenLayer+1)=this%nOutput
      maxNeuronsInEachLayer=maxval(nNeurons)
      if(nInput>maxNeuronsInEachLayer)maxNeuronsInEachLayer=nInput
      allocate(this%a(maxNeuronsInEachLayer,nHiddenLayer+1))
      allocate(this%z(maxNeuronsInEachLayer,0:nHiddenLayer))
      allocate(this%bias(maxNeuronsInEachLayer,nHiddenLayer+1))
      allocate(this%w(maxNeuronsInEachLayer,maxNeuronsInEachLayer,0:nHiddenLayer))

      allocate(this%ag(maxNeuronsInEachLayer,nHiddenLayer+1))
      allocate(this%zg(maxNeuronsInEachLayer,0:nHiddenLayer))
      allocate(this%biasg(maxNeuronsInEachLayer,nHiddenLayer+1))
      allocate(this%wg(maxNeuronsInEachLayer,maxNeuronsInEachLayer,0:nHiddenLayer))
      if(actvfunc == 'tanh') then
        this%ActivateFunc => tanh
        this%DeactivateFunc => tanhgradient
      else if(actvfunc == 'sigmoid') then
        this%ActivateFunc => sigmoid
        this%DeactivateFunc => sigmoidgradient
      else
        print*, 'Wrong input for activation function'
        stop
      end if
    end subroutine initann

    subroutine destroy(this)
      implicit none
      class(NeuralNetwork_t) :: this
      deallocate(this%wg,this%biasg,this%ag,this%zg)
      deallocate(this%w,this%bias,this%a,this%z)
      deallocate(this%nNeurons)
    end subroutine destroy

    subroutine randomize(this)
      use random_m
      implicit none
      class(NeuralNetwork_t) :: this
      integer(kind=4) :: i, j, k
      do i = 0, this%nHiddenLayer
        do j = 1, this%nNeurons(i)
          do k = 1, this%nNeurons(i+1)
            this%w(j,k,i) = MyUniformRand()
          end do
        end do
      end do 
      do i = 1, this%nHiddenLayer+1
        do j = 1, this%nNeurons(i)
          this%bias(j,i) = MyUniformRand()
        end do
      end do
    end subroutine randomize

    subroutine readinput(this,input)
      implicit none
      class(NeuralNetwork_t) :: this
      real(kind=fp_kind), intent(in) :: input(this%nInput)
      this%z(1:this%nInput,0)=input(1:this%nInput)
    end subroutine readinput

    function propagate(this,input) result(output)
      implicit none
      class(NeuralNetwork_t) :: this
      real(kind=fp_kind), intent(in) :: input(this%nInput)
      real(kind=fp_kind) :: output(this%nOutput)
      integer(kind=4) :: iLayer, iNeuron
      integer(kind=4) :: i, j, k 
      integer(kind=4) :: idx_layer, jdx_neuron, kdx_neuron
      this%a=0.d0
      this%z=0.d0
      this%z(1:this%nInput,0)=input(1:this%nInput)
      do idx_layer = 1, this%nHiddenLayer
        do kdx_neuron = 1, this%nNeurons(idx_layer)
          do jdx_neuron = 1, this%nNeurons(idx_layer-1)
            this%a(kdx_neuron,idx_layer) = this%a(kdx_neuron,idx_layer) &
                                         + this%z(jdx_neuron,idx_layer-1)*this%w(jdx_neuron,kdx_neuron,idx_layer-1)  
          end do
          this%a(kdx_neuron,idx_layer) = this%a(kdx_neuron,idx_layer) + this%bias(kdx_neuron,idx_layer)
          this%z(kdx_neuron,idx_layer) = this%ActivateFunc(this%a(kdx_neuron,idx_layer))
!          write(30,*)this%a(kdx_neuron,idx_layer), this%bias(kdx_neuron,idx_layer), this%z(kdx_neuron,idx_layer)
        end do
      end do  
      idx_layer = this%nHiddenLayer + 1
      do kdx_neuron = 1, this%nNeurons(idx_layer)
        do jdx_neuron = 1, this%nNeurons(idx_layer-1)
            this%a(kdx_neuron,idx_layer) = this%a(kdx_neuron,idx_layer) &
                                         + this%z(jdx_neuron,idx_layer-1)*this%w(jdx_neuron,kdx_neuron,idx_layer-1)
        end do
        this%a(kdx_neuron,idx_layer) = this%a(kdx_neuron,idx_layer) + this%bias(kdx_neuron,idx_layer)
!        write(30,*)this%a(kdx_neuron,idx_layer), this%bias(kdx_neuron,idx_layer)
      end do
      output = this%a(1:this%nOutput,this%nHiddenLayer + 1)
!      write(30,*)'------'
    end function propagate

    subroutine backpropagate(this,dLoss_dOutput)
      implicit none
      class(NeuralNetwork_t) :: this
      real(kind=fp_kind) :: dLoss_dOutput
      integer(kind=4) :: idx_atom, idx_ss
      integer(kind=4) :: idx_layer
      integer(kind=4) :: jdx_neuron, kdx_neuron

      this%zg=0.d0
      this%ag=0.d0

      this%ag(1,this%nHiddenLayer+1) = dLoss_dOutput
      this%biasg(1,this%nHiddenLayer+1) = this%ag(1,this%nHiddenLayer+1)

!      write(70,*)this%ctag,this%biasg(1,this%nHiddenLayer+1)

      idx_layer = this%nHiddenLayer
      do jdx_neuron = 1, this%nNeurons(idx_layer)
!        Basicly, there should be this cycle, and dLoss_dOutput should be an array
!        dLoss_dOutput(this%nOutput)
        kdx_neuron = 1
!        do kdx_neuron = 1, this%nNeurons(idx_layer+1) 
          this%wg(jdx_neuron,kdx_neuron,idx_layer) = &
                                                    this%ag(kdx_neuron,idx_layer+1)*this%z(jdx_neuron,idx_layer) 
          this%zg(jdx_neuron,idx_layer) = this%zg(jdx_neuron,idx_layer) &
                                         + this%ag(kdx_neuron,idx_layer+1)*this%w(jdx_neuron,kdx_neuron,idx_layer)
!        end do
        this%ag(jdx_neuron,idx_layer) = this%zg(jdx_neuron,idx_layer) * this%DeactivateFunc(this%a(jdx_neuron,idx_layer))

        this%biasg(jdx_neuron,idx_layer) = this%ag(jdx_neuron,idx_layer)

!        write(77,'(3I4,3G12.4)')idx_layer,jdx_neuron,kdx_neuron,this%wg(jdx_neuron,kdx_neuron,idx_layer),this%ag(kdx_neuron,idx_layer+1),this%z(jdx_neuron,idx_layer)
!        write(26,'(2I4,3G12.4)')idx_layer,jdx_neuron,this%ag(jdx_neuron,idx_layer), this%zg(jdx_neuron,idx_layer), this%DeactivateFunc(this%a(jdx_neuron,idx_layer))

      end do

      do idx_layer = this%nHiddenLayer-1, 0, -1
        do jdx_neuron = 1, this%nNeurons(idx_layer)
          do kdx_neuron = 1, this%nNeurons(idx_layer+1)
            this%wg(jdx_neuron,kdx_neuron,idx_layer) = this%wg(jdx_neuron,kdx_neuron,idx_layer) &
                                                     + this%ag(kdx_neuron,idx_layer+1) * this%z(jdx_neuron,idx_layer)

!            write(27,'(3I4,3G12.4)')idx_layer,jdx_neuron,kdx_neuron,this%wg(jdx_neuron,kdx_neuron,idx_layer),this%ag(kdx_neuron,idx_layer+1),this%z(jdx_neuron,idx_layer)

            if(idx_layer>0)then
              this%zg(jdx_neuron,idx_layer) = this%zg(jdx_neuron,idx_layer) &
                                            + this%ag(kdx_neuron,idx_layer+1)*this%w(jdx_neuron,kdx_neuron,idx_layer)
!              write(29,'(3I4,3G12.4)')idx_layer,jdx_neuron,kdx_neuron,this%zg(jdx_neuron,idx_layer),this%ag(kdx_neuron,idx_layer+1), this%w(jdx_neuron,kdx_neuron,idx_layer)
            end if
          end do
          if(idx_layer>0)then
            this%ag(jdx_neuron,idx_layer) = this%zg(jdx_neuron,idx_layer)*this%DeactivateFunc(this%a(jdx_neuron,idx_layer))
            this%biasg(jdx_neuron,idx_layer) = this%ag(jdx_neuron,idx_layer)
!            write(28,'(2I4,3G12.4)')idx_layer,jdx_neuron,this%ag(jdx_neuron,idx_layer), this%zg(jdx_neuron,idx_layer), this%DeactivateFunc(this%a(jdx_neuron,idx_layer))

          end if
        end do
      end do
    end subroutine backpropagate

end module ann_m

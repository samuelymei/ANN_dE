function sigmoid(x)
  use precision_m
  implicit none
  real(kind=fp_kind), intent(in) :: x
  real(kind=fp_kind) :: sigmoid
  sigmoid = 1.0d0/(1+exp(-x))
end function sigmoid

function sigmoidgradient(x)
  use precision_m
  implicit none
  real(kind=fp_kind), intent(in) :: x
  real(kind=fp_kind) :: sigmoidgradient
  real(kind=fp_kind) :: sigmoid
  sigmoidgradient = sigmoid(x)*(1-sigmoid(x))
end function sigmoidgradient

function tanh(x)
  use precision_m
  implicit none
  real(kind=fp_kind), intent(in) :: x
  real(kind=fp_kind) :: tanh
  tanh = (1-exp(-2*x))/(1+exp(-2*x))
end function tanh

function tanhgradient(x)
  use precision_m
  implicit none
  real(kind=fp_kind), intent(in) :: x
  real(kind=fp_kind) :: tanhgradient
  tanhgradient = 4*exp(-2*x)/(1+exp(-2*x))**2
end function tanhgradient



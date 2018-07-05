SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,mol_training,n_mol_training,ann,nann,idx_ann_for_atom,target_training)
  use precision_m
  use molecule_m
  use ann_m
  implicit none
  INTEGER(kind=4) :: n
  type(molecule_t), intent(in) :: mol_training(n_mol_training)
  integer(kind=4), intent(in) :: n_mol_training
  type(NeuralNetwork_t), intent(inout) :: ann(nann)
  integer(kind=4), intent(in) :: nann
  integer(kind=4), intent(in) :: idx_ann_for_atom(mol_training(1)%num_atoms)
  real(kind=fp_kind), intent(in) :: target_training(n_mol_training)
  real(kind=fp_kind) :: predicted_training(n_mol_training)

  LOGICAL :: check
  REAL(kind=fp_kind) :: f,fold,stpmax,g(n),p(n),x(n),xold(n),ALF,TOLX
  PARAMETER (ALF=1.e-4,TOLX=1.e-7)
  INTEGER(kind=4) :: i
  REAL(kind=fp_kind) :: a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam 

  check=.false.
  sum=0.
  do i=1,n
    sum=sum+p(i)*p(i)
  end do
  sum=sqrt(sum)
  if(sum.gt.stpmax)then
    do i=1,n
      p(i)=p(i)*stpmax/sum
    end do
  endif
  slope=0.
  do i=1,n
    slope=slope+g(i)*p(i)
  end do
  test=0.
  do i=1,n
    temp=abs(p(i))/max(abs(xold(i)),1.)
    if(temp.gt.test)test=temp
  end do
  alamin=TOLX/test
  alam=1.
  do while(.true.)
    do i=1,n
      x(i)=xold(i)+alam*p(i)
    end do
    call loss_func(x,n,mol_training,n_mol_training,target_training,ann,nann,idx_ann_for_atom,predicted_training,f)
    if(alam.lt.alamin)then
      do i=1,n
        x(i)=xold(i)
      end do
      check=.true.
      return
    else if(f.le.fold+ALF*alam*slope)then
      return
    else
      if(alam.eq.1.)then
        tmplam=-slope/(2.*(f-fold-slope))
      else
        rhs1=f-fold-alam*slope
        rhs2=f2-fold2-alam2*slope
        a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
        b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
        if(a.eq.0.)then
          tmplam=-slope/(2.*b)
        else
          disc=b*b-3.*a*slope
          if(disc.lt.0.) pause 'roundoff problem in lnsrch'
          tmplam=(-b+sqrt(disc))/(3.*a)
        endif
        if(tmplam.gt..5*alam)tmplam=.5*alam
      endif
    endif
    alam2=alam
    f2=f
    fold2=fold
    alam=max(tmplam,.1*alam)
  end do
END


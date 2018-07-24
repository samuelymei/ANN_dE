SUBROUTINE dfpmin(p,n,gtol,iter,fret,ann,nann,mol_training,n_mol_training,mol_test,n_mol_test,actvfunc)
  use precision_m
  use molecule_m
  use ann_m
  implicit none

  real(kind=fp_kind), intent(in out) :: p(n)
  integer(kind=4), intent(in) :: n
  real(kind=fp_kind), intent(in) :: gtol
  integer(kind=4), intent(out) :: iter
  real(kind=fp_kind), intent(out) :: fret

  type(NeuralNetwork_t), intent(in out) :: ann(nann)
  integer(kind=4), intent(in) :: nann

  type(molecule_t), intent(in out) :: mol_training(n_mol_training)
  integer(kind=4), intent(in) :: n_mol_training

  type(molecule_t), intent(in out) :: mol_test(n_mol_test)
  integer(kind=4), intent(in) :: n_mol_test

  character(len=20), intent(in) :: actvfunc
  
  real(kind=fp_kind) :: loss_training, loss_test
  real(kind=fp_kind) :: pre_loss_test
  real(kind=fp_kind) :: min_loss_test

  integer(kind=4), parameter :: ITMAX = 2000
  real(kind=fp_kind), parameter :: EPS = 3.E-8, STPMX = 10, TOLX = 4.*EPS
  INTEGER(kind=4) ::  i,its,j
  LOGICAL :: check
  REAL(kind=fp_kind) :: den,fac,fad,fae,fp,stpmax,sum,sumdg,sumxi,temp,test,dg(n), &
       &g(n),hdg(n),hessin(n,n),pnew(n),xi(n)
  real(kind=fp_kind) :: fp2

  integer(kind=4), parameter :: MAX_coupon = 10
  integer(kind=4) :: n_coupon = MAX_coupon
  
  call dloss_func(p,n,ann,nann,mol_training,n_mol_training,fp,mol_test,n_mol_test,fp2,actvfunc,g)
  write(*,'(A,I,A,G12.5,1X,G12.5)') 'Iteration ', 0, ' Loss functions ', sqrt(fp), sqrt(fp2)
  pre_loss_test = fp2
  min_loss_test = fp2
  sum=0.
  hessin=0.d0
  do i=1,n
    hessin(i,i)=1.
    xi(i)=-g(i)
    sum=sum+p(i)**2
  end do
  stpmax=STPMX*max(sqrt(sum),float(n))
  do its=1,ITMAX
    iter=its
    call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check,ann,nann,mol_training,n_mol_training,mol_test,n_mol_test)
    fp=fret
    do i=1,n
      xi(i)=pnew(i)-p(i)
      p(i)=pnew(i)
    end do
    test=0.
    do i=1,n
      temp=abs(xi(i))/max(abs(p(i)),1.)
      if(temp.gt.test)test=temp
    end do
    if(test.lt.TOLX)return
    do i=1,n
      dg(i)=g(i)
    end do
    call dloss_func(p,n,ann,nann,mol_training,n_mol_training,loss_training,mol_test,n_mol_test,loss_test,actvfunc,g)
    write(*,'(A,I,A,G12.5,1X,G12.5,A,1X,G12.5,A,I2)') 'Iteration ', iter, ' Loss functions ', sqrt(loss_training), sqrt(loss_test), 'Norm of dfdv: ', dot_product(g,g), ' # coupon ', n_coupon
    if(loss_test < min_loss_test)then
      min_loss_test = loss_test
    end if
    if(loss_test>loss_training)then
!      if(loss_test > min_loss_test*1.05d0**2)return
      if(loss_test <= pre_loss_test)then
        pre_loss_test = loss_test
        if(n_coupon < MAX_coupon) n_coupon = n_coupon + 1
      else
        pre_loss_test = loss_test
        write(*,*)'Loss function for test set increased. Coupon reduced by 1'
        if(n_coupon > 0 )then
          n_coupon = n_coupon - 1
        else
          write(*,*)'Coupon used up. Training stops.'
!          return
        end if
      end if
    end if
    test=0.
    den=max(fret,1.)
    do i=1,n
      temp=abs(g(i))*max(abs(p(i)),1.)/den
      if(temp.gt.test)test=temp
    end do
    if(test.lt.gtol)return
    do i=1,n
      dg(i)=g(i)-dg(i)
    end do
    do i=1,n
      hdg(i)=0.
      do j=1,n
        hdg(i)=hdg(i)+hessin(i,j)*dg(j)
      end do
    end do
    fac=0.
    fae=0.
    sumdg=0.
    sumxi=0.
    do i=1,n
      fac=fac+dg(i)*xi(i)
      fae=fae+dg(i)*hdg(i)
      sumdg=sumdg+dg(i)**2
      sumxi=sumxi+xi(i)**2
    end do
    if(fac**2.gt.EPS*sumdg*sumxi)then
      fac=1./fac
      fad=1./fae
      do i=1,n
        dg(i)=fac*xi(i)-fad*hdg(i)
      end do
      do i=1,n
        do j=1,n
          hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j)
        end do
      end do
    endif
    do i=1,n
      xi(i)=0.
      do j=1,n
        xi(i)=xi(i)-hessin(i,j)*g(j)
      end do
    end do
    rewind(18)
    write(18)n
    write(18)p
  end do
  pause 'too many iterations in dfpmin'
  return
END


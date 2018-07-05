SUBROUTINE dfpmin(p,n,gtol,iter,fret,mol_training,n_mol_training,ann,nann,idx_ann_for_atom,target_training,predicted_training)
  use precision_m
  use molecule_m
  use ann_m
  implicit none
  type(molecule_t), intent(in) :: mol_training(n_mol_training)
  integer(kind=4), intent(in) :: n_mol_training
  type(NeuralNetwork_t), intent(inout) :: ann(nann)
  integer(kind=4), intent(in) :: nann
  integer(kind=4), intent(in) :: idx_ann_for_atom(mol_training(1)%num_atoms)
  real(kind=fp_kind), intent(in) :: target_training(n_mol_training)
  real(kind=fp_kind), intent(out) :: predicted_training(n_mol_training)
  
  INTEGER(kind=4) :: iter,n,ITMAX
  REAL(kind=fp_kind) :: fret,gtol,p(n),EPS,STPMX,TOLX
  PARAMETER (ITMAX=1000000,STPMX=100.,EPS=3.e-8,TOLX=4.*EPS)
  INTEGER(kind=4) ::  i,its,j
  LOGICAL :: check
  REAL(kind=fp_kind) :: den,fac,fad,fae,fp,stpmax,sum,sumdg,sumxi,temp,test,dg(n), &
       &g(n),hdg(n),hessin(n,n),pnew(n),xi(n)
  
  call loss_func(p,n,mol_training,n_mol_training,target_training,ann,nann,idx_ann_for_atom,predicted_training,fp) !function called here
  call dloss_func(p,n,mol_training,n_mol_training,target_training,ann,nann,idx_ann_for_atom,predicted_training,g)
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
    call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check,mol_training,n_mol_training,ann,nann,idx_ann_for_atom,target_training)
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
    call dloss_func(p,n,mol_training,n_mol_training,target_training,ann,nann,idx_ann_for_atom,predicted_training,g)
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
  end do
  pause 'too many iterations in dfpmin'
  return
END


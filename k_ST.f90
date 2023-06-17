program rate
! program to calculate ISC and RISC rates
  implicit none
  integer::i,j,k,n
  double precision::est,lambdam,hoc,kbt,rate_const
  double precision::rho,hbar,pi
  double precision::lambda,lambdai,fcwd,omega,S
  double precision::fact
  character::coupling
  write(*,*)
  write(*,*)'SINGLET-TRIPLET GAP IS NEGETIVE FOR RISC'
  write(*,*)
  write(*,*)'singlet-triplet gap (in eV), lambda_M(in eV) and HOC (in cm^-1)'
  write(*,*)'(e.g., -0.15 0.1 0.09 for RISC of CC2TA with lambda_M=0.1 eV)'
  read(*,*)est,lambdam,hoc
  omega=0.20
  S=0.00d0
  kbt=0.02570d0
  hbar=6.582119*(10**(-16.0))
  pi=4.0d0*atan(1.0d0)
  hoc=hoc*0.00012390d0
!---------------------------------------------
  write(*,*)'electron vibration coupling (yes/no)'
  read(*,'(a1)')coupling
  if(coupling.eq.'y')then
     write(*,*)'Effective S and omega(eV)'
     read(*,*)S,omega
  endif
!---------------------------------------------
  fcwd=0.0d0
  do i=0,10
  rho=exp(-(((-est+i*omega+lambdam)**2)/(4.0d0*lambdam*kbt)))
  rho=rho*(exp(-S))*(S**i)/fact(i)
  fcwd=fcwd+rho
  write(*,*)'i=',i,'rho=',rho
  enddo
  fcwd=fcwd/sqrt(4.0d0*pi*lambdam*kbt)
  rate_const=((2.0d0*pi)/hbar)*fcwd*hoc*hoc
  write(*,"(E10.4,2x,A8)")rate_const,'(in s-1)'
end program rate

function fact(n)
integer ::n
double precision::fact,p
p = 1.0d0
do i = 1, n
p = p * real(i)
end do
fact = p
end

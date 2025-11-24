! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

module expm_cheb_mod

#include "../preproc.inc"

use constants
use expm_utils
use lapack
use profiling

use standalone, only: eyef_r, eyef_c

implicit none

private
public :: expm_cheb, expm_cheb_act, expm_cheb_clen, expm_cheb_clen_act

interface expm_cheb
  module procedure expm_cheb_r, expm_cheb_c, expm_cheb_cr
end interface expm_cheb

interface expm_cheb_act
  module procedure expm_cheb_act_r, expm_cheb_act_c, expm_cheb_act_cr
end interface expm_cheb_act

interface expm_cheb_clen
  module procedure expm_cheb_clen_r, expm_cheb_clen_c, expm_cheb_clen_cr
end interface expm_cheb_clen

interface expm_cheb_clen_act
  module procedure expm_cheb_clen_act_r, expm_cheb_clen_act_c, expm_cheb_clen_act_cr
end interface expm_cheb_clen_act

real(wp), dimension(0:50), parameter :: cheb_coeff =[1.2660658777520082_wp, 1.13031820798497_wp, 0.2714953395340767_wp, &
  0.04433684984866381_wp, 0.005474240442093732_wp, 0.0005429263119139442_wp, 4.497732295429519e-05_wp, &
  3.1984364624019926e-06_wp, 1.9921248066727963e-07_wp, 1.103677172551734e-08_wp, 5.50589607967375e-10_wp, &
  2.4979566169849818e-11_wp, 1.039152230678574e-12_wp, 3.991263356414398e-14_wp, 1.4237580108256627e-15_wp, &
  4.7409261025615024e-17_wp, 1.4801800572082964e-18_wp, 4.349919494944147e-20_wp, 1.2074289272797564e-21_wp, &
  3.175356737059443e-23_wp, 7.933671971638028e-25_wp, 1.8879484042289102e-26_wp, 4.288673876592579e-28_wp, &
  9.318985281777475e-30_wp, 1.9406469749017366e-31_wp, 3.879802249226021e-33_wp, 7.458502887391588e-35_wp, &
  1.3807477824110445e-36_wp, 2.4648623717710954e-38_wp, 4.248542192505979e-40_wp, 7.079001176212945e-42_wp, &
  1.141486778254092e-43_wp, 1.7831510375432693e-45_wp, 2.7011422638737947e-47_wp, 3.9714338657414944e-49_wp, &
  5.672351695612143e-51_wp, 7.876788130533234e-53_wp, 1.0642416282087371e-54_wp, 1.4000817885528946e-56_wp, &
  1.794689085390564e-58_wp, 2.2430194826630055e-60_wp, 2.7349926016641725e-62_wp, 3.2554929833033565e-64_wp, &
  3.784956894021595e-66_wp, 4.3005444470977184e-68_wp, 4.777805756025809e-70_wp, 5.192666745513436e-72_wp, &
  5.523501533951769e-74_wp, 5.75303598890998e-76_wp, 5.869845981477718e-78_wp, 5.8692706170236985e-80_wp]

real(wp), dimension(0:120), parameter :: cheb_coeff_f = [1.2660658777520082_wp, 1.13031820798497_wp, 0.24019372387008986_wp, &
  0.1633061176105341_wp, 0.12346931414340684_wp, 0.09917838239971265_wp, 0.0828424078319937_wp, 0.0711122017122309_wp, &
  0.06228433267599542_wp, 0.05540200939493705_wp, 0.04988683481551001_wp, 0.045368757071292155_wp, 0.04160009119505143_wp, &
  0.03840884173253493_wp, 0.03567186335969355_wp, 0.03329867903473396_wp, 0.03122132733536094_wp, 0.02938777261428817_wp, &
  0.027757500539564803_wp, 0.026298498117096426_wp, 0.024985135934632183_wp, 0.02379665318881484_wp, 0.02271605445883036_wp, &
  0.021729293366511606_wp, 0.02082465972659615_wp, 0.019992313385191928_wp, 0.019223925365988638_wp, 0.018512398577268958_wp, &
  0.017851648238514522_wp, 0.01723642764465281_wp, 0.016662188711929527_wp, 0.016124969467299247_wp, 0.015621302598621465_wp, &
  0.015148140606166961_wp, 0.01470279414326713_wp, 0.014282880912466296_wp, 0.013886283067790625_wp, 0.013511111516169362_wp, &
  0.013155675848814727_wp, 0.012818458893358867_wp, 0.01249809507909764_wp, 0.01219335196508002_wp, 0.011903114404486773_wp, &
  0.011626370916582317_wp, 0.01136220191540491_wp, 0.011109769506626484_wp, 0.01086830861418843_wp, 0.01063711923882691_wp, &
  0.01041555968355818_wp, 0.010203040608111804_wp, 0.009999019796335653_wp, 0.009802997538795636_wp, 0.009614512547818168_wp, &
  0.009433138334693494_wp, 0.009258479989174357_wp, 0.009090171310108076_wp, 0.008927872243325821_wp, 0.008771266589089357_wp, &
  0.008620059946575684_wp, 0.008473977867293855_wp, 0.00833276419307036_wp, 0.008196179557440956_wp, 0.008064000031995927_wp, &
  0.007936015901597766_wp, 0.007812030554365281_wp, 0.00769185947408736_wp, 0.0075753293241921146_wp, 0.007462277113721257_wp, &
  0.007352549436865402_wp, 0.007246001778606176_wp, 0.007142497879848231_wp, 0.00704190915618933_wp, 0.006944114165105254_wp, &
  0.006848998116918449_wp, 0.006756452425399983_wp, 0.006666374294319106_wp, 0.006578666336614472_wp, 0.006493236223243149_wp, &
  0.006409996359022867_wp, 0.006328863583091979_wp, 0.006249758891833616_wp, 0.006172607182307156_wp, 0.006097337014457916_wp, &
  0.0060238803905090016_wp, 0.005952172550105096_wp, 0.005882151779915934_wp, 0.005813759236520822_wp, 0.005746938781499922_wp, &
  0.005681636827773729_wp, 0.005617802196294478_wp, 0.005555385982292489_wp, 0.005494341430336254_wp, 0.005434623817542903_wp, &
  0.005376190344319191_wp, 0.005319000032076422_wp, 0.005263013627402863_wp, 0.005208193512227507_wp, 0.005154503619534835_wp, &
  0.0051019093542435156_wp, 0.005050377518875764_wp, 0.004999876243690263_wp, 0.004950374920959055_wp, 0.004901844143114129_wp, &
  0.004854255644498935_wp, 0.0048075822464684476_wp, 0.004761797805641684_wp, 0.0047168771650693345_wp, 0.00467279610814952_wp, &
  0.0046295313150958835_wp, 0.0045870603218011435_wp, 0.004545361480956779_wp, 0.004504413925254491_wp, 0.004464197532575004_wp, &
  0.004424692893019043_wp, 0.004385881277668148_wp, 0.004347744608985669_wp, 0.004310265432740268_wp, 0.004273426891373037_wp, &
  0.004237212698715943_wp, 0.004201607115988343_wp, 0.004166594928987698_wp]

real(wp), dimension(0:120), parameter :: cheb_coeff_g = [1.2660658777520082_wp, 1.13031820798497_wp, 0.2144401364138621_wp, &
  0.03922510451964104_wp, 0.02016329433679518_wp, 0.012245486852745053_wp, 0.008216176002874422_wp, 0.005891106016075633_wp, &
  0.004429176028767081_wp, 0.0034506771840728825_wp, 0.002763830891132558_wp, 0.0022633036898005534_wp, 0.0018873444315718883_wp,&
  0.0015978113187697502_wp, 0.0013701149540870814_wp, 0.0011878259285853223_wp, 0.001039628957978550_wp, 0.000917525268447847_wp, &
  0.0008157311141977117_wp, 0.0007299805756750481_wp, 0.0006570715503323226_wp, 0.000594562614711837_wp, 0.000540566069775017_wp, &
  0.0004936038114655789_wp, 0.0004525051404569871_wp, 0.0004163331233940955_wp, 0.000384330740410385_wp, 0.000355880968594853_wp, &
  0.0003304768274525821_wp, 0.0003076986432009493_wp, 0.0002871966101347241_wp, 0.000268677284238242_wp, 0.000251893027442214_wp, &
  0.0002366336882153993_wp, 0.0002227199929857386_wp, 0.0002099982578287913_wp, 0.000198336127374051_wp, 0.000187619119074014_wp, &
  0.0001777478034639118_wp, 0.0001686354900823856_wp, 0.0001602063180167039_wp, 0.000152393672192472_wp, 0.000145138863414521_wp, &
  0.0001383900231290771_wp, 0.0001321011738975995_wp, 0.0001262314443678985_wp, 0.000120744403630517_wp, 0.000115607494653492_wp, &
  0.0001107915502931266_wp, 0.0001062703784075562_wp, 0.0001020204050233265_wp, 9.802036645384723e-05_wp, &
  9.425104284298124e-05_wp, 9.069502688421516e-05_wp, 8.733652250687324e-05_wp, 8.416116917280247e-05_wp, &
  8.11558881265906e-05_wp, 7.830874751954202e-05_wp, 7.560884380534668e-05_wp, 7.30461972020286e-05_wp, &
  7.061165934545697e-05_wp, 6.829683153621926e-05_wp, 6.609399221344822e-05_wp, 6.399603248440457e-05_wp, &
  6.199639870321048e-05_wp, 6.008904123145452e-05_wp, 5.826836863161891e-05_wp, 5.6529206644820335e-05_wp, &
  5.486676139022481e-05_wp, 5.327658629681655e-05_wp, 5.175455234107112e-05_wp, 5.0296821218166136e-05_wp, &
  4.8899821120878716e-05_wp, 4.7560224840472614e-05_wp, 4.627492993861357e-05_wp, 4.504104076947642e-05_wp, &
  4.385585215730916e-05_wp, 4.2716834557535406e-05_wp, 4.162162054926398e-05_wp, 4.0567992524372e-05_wp, &
  3.955387145363105e-05_wp, 3.857730662342018e-05_wp, 3.76364662483902e-05_wp, 3.672962887571774e-05_wp, &
  3.585517550550404e-05_wp, 3.5011582359967453e-05_wp, 3.4197414241103657e-05_wp, 3.3411318422664886e-05_wp, &
  3.265201902793104e-05_wp, 3.1918311849614845e-05_wp, 3.12090595725863e-05_wp, 3.052318736401889e-05_wp, &
  2.9859678799018144e-05_wp, 2.9217572092881252e-05_wp, 2.8595956613882726e-05_wp, 2.799396965297447e-05_wp, &
  2.7410793429004542e-05_wp, 2.6845652310014528e-05_wp, 2.6297810232986832e-05_wp, 2.576656830601342e-05_wp, &
  2.5251262578294302e-05_wp, 2.475126196466324e-05_wp, 2.426596631252221e-05_wp, 2.3794804600165807e-05_wp, &
  2.333723325631233e-05_wp, 2.289273459167537e-05_wp, 2.2460815334108522e-05_wp, 2.2041005259555325e-05_wp, &
  2.1632855911736376e-05_wp, 2.1235939404012196e-05_wp, 2.0849847297540126e-05_wp, 2.0474189550137092e-05_wp, &
  2.0108593530817584e-05_wp, 1.975270309541777e-05_wp, 1.9406177718923535e-05_wp, 1.906869168063287e-05_wp, &
  1.873993329849378e-05_wp, 1.8419604209227902e-05_wp, 1.810741869116003e-05_wp, 1.7803103026881076e-05_wp, &
  1.7506394903075658e-05_wp]

contains

!******************************************************************************** 
!
! Chebyshev expansion to approximate matrix exponential
!
! INPUT:  
!         f - prefactor in exponential
!         A - matrix to exponentiate
!         kmax - truncation order of the Chebyshev expansion
!         tol - stopping criterion
!
! OUTPUT:  
!         expm_A - exp(f*A)
!
!******************************************************************************** 
subroutine expm_cheb_r(f, A, B, kmax, tol)
  real(wp), intent(in)           :: f
  real(wp), intent(in)           :: A(:,:)
  real(wp), intent(out)          :: B(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local
  integer                        :: i, n, kmax_
  real(wp)                       :: res
  real(wp)                       :: coeff1, coeff2
  real(wp), allocatable          :: X(:,:), T(:,:), TT(:,:), TTT(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(A, 1) 
  allocate(X(n,n), T(n,n), TT(n,n), TTT(n,n))

  X = f * A
  TTT = cheb_coeff(0) * eyef_r(n)
  TT = cheb_coeff(1) * X
  B = TTT + TT
  
  do i = 2, kmax_
    T = TTT
    coeff1 = 2.0_wp * cheb_coeff_f(i)
    coeff2 = -cheb_coeff_g(i)
    call gemm("n", "n", n, n, n, coeff1, X, n, TT, n, coeff2, T, n)
    B = B + T
    if (present(tol)) then
      res = maxval(abs(T))
      if (res < tol) return
    end if
    TTT = TT 
    TT = T
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_cheb_r did not converge: res tol = ", res, tol
end subroutine expm_cheb_r

subroutine expm_cheb_c(f, A, B, kmax, tol)
  complex(wp), intent(in)        :: f
  complex(wp), intent(in)        :: A(:,:)
  complex(wp), intent(out)       :: B(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local
  integer                        :: i, n, kmax_
  real(wp)                       :: res
  complex(wp)                    :: coeff1, coeff2
  complex(wp), allocatable       :: X(:,:), T(:,:), TT(:,:), TTT(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(A, 1) 
  allocate(X(n,n), T(n,n), TT(n,n), TTT(n,n))

  X = f * A
  TTT = cheb_coeff(0) * eyef_c(n)
  TT = cheb_coeff(1) * X
  B = TTT + TT
  
  do i = 2, kmax_
    T = TTT
    coeff1 = (2.0_wp,0.0_wp) * cheb_coeff_f(i)
    coeff2 = -(1.0_wp,0.0_wp) * cheb_coeff_g(i)
    call gemm("n", "n", n, n, n, coeff1, X, n, TT, n, coeff2, T, n)
    B = B + T
    if (present(tol)) then
      res = maxval(abs(T))
      if (res < tol) return
    end if
    TTT = TT 
    TT = T
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_cheb_c did not converge: res tol = ", res, tol
end subroutine expm_cheb_c

subroutine expm_cheb_cr(f, A, B, kmax, tol)
  complex(wp), intent(in)        :: f
  real(wp), intent(in)           :: A(:,:)
  complex(wp), intent(out)       :: B(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local
  integer                        :: i, n, kmax_
  real(wp)                       :: res
  complex(wp)                    :: coeff1, coeff2
  complex(wp), allocatable       :: X(:,:), T(:,:), TT(:,:), TTT(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(A, 1) 
  allocate(X(n,n), T(n,n), TT(n,n), TTT(n,n))

  X = f * A
  TTT = cheb_coeff(0) * eyef_c(n)
  TT = cheb_coeff(1) * X
  B = TTT + TT
  
  do i = 2, kmax_
    T = TTT
    coeff1 = (2.0_wp,0.0_wp) * cheb_coeff_f(i)
    coeff2 = -(1.0_wp,0.0_wp) * cheb_coeff_g(i)
    call gemm("n", "n", n, n, n, coeff1, X, n, TT, n, coeff2, T, n)
    B = B + T
    if (present(tol)) then
      res = maxval(abs(T))
      if (res < tol) return
    end if
    TTT = TT 
    TT = T
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_cheb_cr did not converge: res tol = ", res, tol
end subroutine expm_cheb_cr


!******************************************************************************** 
!
! Chebyshev expansion to approximate action of the matrix exponetial
!
! INPUT:
!         f - prefactor in exponential
!         A - matrix to exponentiate
!         C - matrix to act on
!         kmax - truncation order of the Taylor series
!         tol - stopping criterion
!
! OUTPUT:
!         C - exp(f*A) C 
!
!******************************************************************************** 
subroutine expm_cheb_act_r(f, A, C, kmax, tol)
  real(wp), intent(in)           :: f
  real(wp), intent(in)           :: A(:,:)
  real(wp), intent(inout)        :: C(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local
  integer                        :: i, n, m, kmax_
  real(wp)                       :: res
  real(wp)                       :: coeff1, coeff2
  real(wp), allocatable          :: X(:,:), T(:,:), TT(:,:), TTT(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(C, 1) 
  m = size(C, 2)
  allocate(X(n,n), T(n,m), TT(n,m), TTT(n,m))

  X = f * A
  TTT = cheb_coeff(0) * C
  call gemm("n", "n", n, m, n, cheb_coeff(1), X, n, C, n, zeror, TT, n)
  C = TTT + TT
  
  do i = 2, kmax_
    T = TTT
    coeff1 = 2.0_wp * cheb_coeff_f(i)
    coeff2 = -cheb_coeff_g(i)
    call gemm("n", "n", n, m, n, coeff1, X, n, TT, n, coeff2, T, n)
    C = C + T
    if (present(tol)) then
      res = maxval(abs(T))
      if (res < tol) return
    end if
    TTT = TT 
    TT = T
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_cheb_act_r did not converge: res tol = ", res, tol
end subroutine expm_cheb_act_r

subroutine expm_cheb_act_c(f, A, C, kmax, tol)
  complex(wp), intent(in)        :: f
  complex(wp), intent(in)        :: A(:,:)
  complex(wp), intent(inout)     :: C(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local
  integer                        :: i, n, m, kmax_
  real(wp)                       :: res
  complex(wp)                    :: coeff1, coeff2
  complex(wp), allocatable       :: X(:,:), T(:,:), TT(:,:), TTT(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(C, 1) 
  m = size(C, 2)
  allocate(X(n,n), T(n,m), TT(n,m), TTT(n,m))

  X = f * A
  TTT = cheb_coeff(0) * C
  coeff1 = onec * cheb_coeff(1)
  call gemm("n", "n", n, m, n, coeff1, X, n, C, n, zeroc, TT, n)
  C = TTT + TT
  
  do i = 2, kmax_
    T = TTT
    coeff1 = (2.0_wp,0.0_wp) * cheb_coeff_f(i)
    coeff2 = -onec * cheb_coeff_g(i)
    call gemm("n", "n", n, m, n, coeff1, X, n, TT, n, coeff2, T, n)
    C = C + T
    if (present(tol)) then
      res = maxval(abs(T))
      if (res < tol) return
    end if
    TTT = TT 
    TT = T
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_cheb_act_c did not converge: res tol = ", res, tol
end subroutine expm_cheb_act_c

subroutine expm_cheb_act_cr(f, A, C, kmax, tol)
  complex(wp), intent(in)        :: f
  real(wp), intent(in)           :: A(:,:)
  complex(wp), intent(inout)     :: C(:,:)
  integer, optional, intent(in)  :: kmax
  real(wp), optional, intent(in) :: tol
  !local
  integer                        :: i, n, m, kmax_
  real(wp)                       :: res
  complex(wp)                    :: coeff1, coeff2
  complex(wp), allocatable       :: X(:,:), T(:,:), TT(:,:), TTT(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  if (present(tol)) kmax_ = max_expm_order

  n = size(C, 1) 
  m = size(C, 2)
  allocate(X(n,n), T(n,m), TT(n,m), TTT(n,m))

  X = f * A
  TTT = cheb_coeff(0) * C
  coeff1 = onec * cheb_coeff(1)
  call gemm("n", "n", n, m, n, coeff1, X, n, C, n, zeroc, TT, n)
  C = TTT + TT
  
  do i = 2, kmax_
    T = TTT
    coeff1 = (2.0_wp,0.0_wp) * cheb_coeff_f(i)
    coeff2 = -onec * cheb_coeff_g(i)
    call gemm("n", "n", n, m, n, coeff1, X, n, TT, n, coeff2, T, n)
    C = C + T
    if (present(tol)) then
      res = maxval(abs(T))
      if (res < tol) return
    end if
    TTT = TT 
    TT = T
  end do  

  if (present(tol)) write(*,"(1x,a,2es12.4)") "expm_cheb_act_cr did not converge: res tol = ", res, tol
end subroutine expm_cheb_act_cr


!******************************************************************************** 
!
! Chebyshev expansion to approximate matrix exponential
!
! Clenshaw algorithm used to evaluate Chebyshev polynomials arXiv:1507.03917v4
!
! INPUT:  
!         f - prefactor in exponential
!         A - matrix to exponentiate
!         kmax - truncation order of the Taylor series
!         tol - stopping criterion
!
! OUTPUT:  
!         expm_A - exp(f*A)
!
!******************************************************************************** 
subroutine expm_cheb_clen_r(f, A, B, kmax)
  real(wp), intent(in)          :: f
  real(wp), intent(in)          :: A(:,:)
  real(wp), intent(out)         :: B(:,:)
  integer, optional, intent(in) :: kmax
  !local
  integer                       :: i, n, kmax_
  real(wp), allocatable         :: T(:,:), X(:,:), mat(:,:) , mat2(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  n = size(A, 1) 
  allocate(T(n,n), X(n,n), mat(n,n), mat2(n,n))

  mat = zeror
  mat2 = zeror
  T = f * A
  
  do i = kmax_, 1, -1
    X = mat
    call gemm("n", "n", n, n, n, 2.0_wp, T, n, mat, n, zeror, B, n)
    mat = B - mat2 + cheb_coeff(i)*eyef_r(n)
    mat2 = X
  end do  

  call gemm("n", "n", n, n, n, oner, T, n, mat, n, zeror, B, n)
  B = B - mat2 + cheb_coeff(0)*eyef_r(n)
end subroutine expm_cheb_clen_r

subroutine expm_cheb_clen_c(f, A, B, kmax)
  complex(wp), intent(in)       :: f
  complex(wp), intent(in)       :: A(:,:)
  complex(wp), intent(out)      :: B(:,:)
  integer, optional, intent(in) :: kmax
  !local
  integer                       :: i, n, kmax_
  complex(wp), allocatable      :: T(:,:), X(:,:), mat(:,:) , mat2(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  n = size(A, 1) 
  allocate(T(n,n), X(n,n), mat(n,n), mat2(n,n))

  mat = zeroc
  mat2 = zeroc
  T = f * A
  
  do i = kmax_, 1, -1
    X = mat
    call gemm("n", "n", n, n, n, (2.0_wp,0.0_wp), T, n, mat, n, zeroc, B, n)
    mat = B - mat2 + cheb_coeff(i)*eyef_c(n)
    mat2 = X
  end do  

  call gemm("n", "n", n, n, n, onec, T, n, mat, n, zeroc, B, n)
  B = B - mat2 + cheb_coeff(0)*eyef_c(n)
end subroutine expm_cheb_clen_c

subroutine expm_cheb_clen_cr(f, A, B, kmax)
  complex(wp), intent(in)       :: f
  real(wp), intent(in)          :: A(:,:)
  complex(wp), intent(out)      :: B(:,:)
  integer, optional, intent(in) :: kmax
  !local
  integer                       :: i, n, kmax_
  complex(wp), allocatable      :: T(:,:), X(:,:), mat(:,:) , mat2(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  n = size(A, 1) 
  allocate(T(n,n), X(n,n), mat(n,n), mat2(n,n))

  mat = zeroc
  mat2 = zeroc
  T = f * A
  
  do i = kmax_, 1, -1
    X = mat
    call gemm("n", "n", n, n, n, (2.0_wp,0.0_wp), T, n, mat, n, zeroc, B, n)
    mat = B - mat2 + cheb_coeff(i)*eyef_c(n)
    mat2 = X
  end do  

  call gemm("n", "n", n, n, n, onec, T, n, mat, n, zeroc, B, n)
  B = B - mat2 + cheb_coeff(0)*eyef_c(n)
end subroutine expm_cheb_clen_cr


!******************************************************************************** 
!
! Chebyshev expansion to approximate action of the matrix exponetial
!
! Clenshaw algorithm used to evaluate Chebyshev polynomials arXiv:1507.03917v4
!
! INPUT:
!         f - prefactor in exponential
!         A - matrix to exponentiate
!         C - matrix to act on
!         kmax - truncation order of the Taylor series
!         tol - stopping criterion
!
! OUTPUT:
!         C - exp(f*A) C 
!
!******************************************************************************** 
subroutine expm_cheb_clen_act_r(f, A, C, kmax)
  real(wp), intent(in)          :: f
  real(wp), intent(in)          :: A(:,:)
  real(wp), intent(inout)       :: C(:,:)
  integer, optional, intent(in) :: kmax
  !local
  integer                       :: i, n, m, kmax_
  real(wp), allocatable         :: T(:,:), B(:,:), X(:,:), mat(:,:) , mat2(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  n = size(C, 1) 
  m = size(C, 2)
  allocate(T(n,n), B(n,m), X(n,m), mat(n,m), mat2(n,m))

  mat = zeror
  mat2 = zeror
  T = f * A
  
  do i = kmax_, 1, -1
    X = mat
    call gemm("n", "n", n, m, n, 2.0_wp, T, n, mat, n, zeror, B, n)
    mat = B - mat2 + cheb_coeff(i)*C
    mat2 = X
  end do  

  call gemm("n", "n", n, m, n, oner, T, n, mat, n, zeror, B, n)
  C = B - mat2 + cheb_coeff(0)*C
end subroutine expm_cheb_clen_act_r

subroutine expm_cheb_clen_act_c(f, A, C, kmax)
  complex(wp), intent(in)       :: f
  complex(wp), intent(in)       :: A(:,:)
  complex(wp), intent(inout)    :: C(:,:)
  integer, optional, intent(in) :: kmax
  !local
  integer                       :: i, n, m, kmax_
  complex(wp), allocatable      :: T(:,:), B(:,:), X(:,:), mat(:,:) , mat2(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  n = size(C, 1) 
  m = size(C, 2)
  allocate(T(n,n), B(n,m), X(n,m), mat(n,m), mat2(n,m))

  mat = zeroc
  mat2 = zeroc
  T = f * A
  
  do i = kmax_, 1, -1
    X = mat
    call gemm("n", "n", n, m, n, (2.0_wp, 0.0_wp), T, n, mat, n, zeroc, B, n)
    mat = B - mat2 + cheb_coeff(i)*C
    mat2 = X
  end do  

  call gemm("n", "n", n, m, n, onec, T, n, mat, n, zeroc, B, n)
  C = B - mat2 + cheb_coeff(0)*C
end subroutine expm_cheb_clen_act_c

subroutine expm_cheb_clen_act_cr(f, A, C, kmax)
  complex(wp), intent(in)       :: f
  real(wp), intent(in)          :: A(:,:)
  complex(wp), intent(inout)    :: C(:,:)
  integer, optional, intent(in) :: kmax
  !local
  integer                       :: i, n, m, kmax_
  complex(wp), allocatable      :: T(:,:), B(:,:), X(:,:), mat(:,:) , mat2(:,:)

  if (present(kmax)) then
    kmax_ = kmax
  else
    kmax_ = default_expm_order
  end if

  n = size(C, 1) 
  m = size(C, 2)
  allocate(T(n,n), B(n,m), X(n,m), mat(n,m), mat2(n,m))

  mat = zeroc
  mat2 = zeroc
  T = f * A
  
  do i = kmax_, 1, -1
    X = mat
    call gemm("n", "n", n, m, n, (2.0_wp, 0.0_wp), T, n, mat, n, zeroc, B, n)
    mat = B - mat2 + cheb_coeff(i)*C
    mat2 = X
  end do  

  call gemm("n", "n", n, m, n, onec, T, n, mat, n, zeroc, B, n)
  C = B - mat2 + cheb_coeff(0)*C
end subroutine expm_cheb_clen_act_cr

end module expm_cheb_mod
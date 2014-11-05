program zeroD

!Declaration of Variables
!###################################################################
implicit none

integer turbmod, casetype, i, j, ntmax, niter, n, m, x, RSM_mod, l

double precision residual, Smag, k0, e0, turbmod3
double precision ep, ka, epa, casa, cas, a2a, a2
double precision eta, sij_total, dt, t
double precision sij(3,3), aij(3,3), dsijdt(3,3), oij(3,3), delta(3,3)
double precision smax, omega, R, dsdt, pi, r12n
parameter (pi=3.14159d0)
double precision alpha_1, alpha_1s, alpha_3, alpha_3s, alpha_4, alpha_5
double precision cssg_1, cssg_1_star, cssg_3, cssg_3_star, cssg_4, cssg_5
double precision obar, snorm, sbar, r11, r22, r33, r12, r11n, r22n, r33n
double precision casn, kn, en, e_old, k_old, prodk, ce1, ce2, e_new
double precision k_new, kdiff, ediff, cmu, s, A
double precision u1u1_n, u2u2_n, u3u3_n, u1u2_n
double precision u1u1_old, u2u2_old, u3u3_old, u1u2_old
double precision u1u1_new, u2u2_new, u3u3_new, u1u2_new
double precision fthd, tthd, othd,  ABAR
double precision a11_old, a22_old, a33_old, a12_old
double precision u1u1_diff, u2u2_diff, u3u3_diff, u1u2_diff
double precision K,E,PROD,S11,S12,S13,S21,S22,S23,S31,S32,S33
double precision A11,A22,A33,A12,O21,SRCR11,SRCR22,SRCR33,SRCR12
double precision O12, cssg_2, aij_Sij, cas_exact, t1, t2, t3, t4, t5, t6
double precision ass, saa, term7, aos, aik_ajk, casdiff, t7, t2a, tt
double precision dcasdt, cas_new, cas_n, cas_old, s0, sbar0
character*6 Output 

!###################################################################
!End of the declaration of variables 


WRITE(*,*) '                     YOU ARE TESTING:'
WRITE(*,*) '**********************************************************'
WRITE(*,*) 'TURBULENCE MODEL:'
!k-epsilon = 1 RSM(SSG) = 2 Cas-k-eps = 3

WRITE(*,*) 'Models: 1 - k-epsilon; 2 - RSM; 3 - Cas-ke'
READ(*,*) turbmod

WRITE(*,*) 'IN CASE: '
WRITE(*,*) '1 - A Constant Shear'
WRITE(*,*) '2 - B Oscillating Shear'
WRITE(*,*) '3 - C Oscillating Strain'
WRITE(*,*) '4 - D Strain-Relax-Destrain'
WRITE(*,*) '5 - E Shear Flow with Rotation'
!casetype = 1 (A - constant shear), 2 (B - osc shear), 3(C - osc strain)
!4 (D - strain-rel-destain), 5 (E - constant rotation)
READ(*,*) casetype



!Code Numerics Definition
!###################################################################
ntmax = 3000   !number of time steps
residual = 1.d-4
niter = 10000000    !number of iteration for each time step 
dt = 1.d-3
!###################################################################
!End of code numerics defintion  

!Case specific parameters
!###################################################################

omega = 2*pi/1.d0
R = 0.75d0
smax = omega/R
!###################################################################
!End of case specific parameters 

!Turbulence Initiation
!###################################################################
k0 = 0.161448d0
eta = 3.3   !non-dimensional strain 

if (int(casetype).eq.1)then
	s0 = 1.65
	sbar0 = 0.5d0*abs(s0)**2
	else if (int(casetype).eq.2)then
		s0 =smax*cos(omega*0.0d0)
		sbar0 = 0.5d0*abs(s0)**2
	else if (int(casetype).eq.3)then
		s0 = smax*cos(omega*0.0d0)/(sin(omega*0.0d0)+5.d0/4.d0)
		sbar0 =	0.5d0*abs(s0)**2	
	else if (int(casetype).eq.4)then
		s0 = 0.0d0
		sbar0 = 0.5d0*abs(s0)**2				
	else if (int(casetype).eq.5)then
		s0 =0.0d0
		sbar0 = 0.5d0*abs(s0)**2
endif
e0 = (s0*k0)/eta
!###################################################################
!End of turbulence initiation set up

!Constants 
!###################################################################
ce1 = 1.44
ce2 = 1.92
a2 = 0.3
cmu = 0.09
fthd = 4.0d0/3.0d0
tthd = 2.0d0/3.0d0
othd = 1.0d0/3.0d0

!RSM is set up for SSG
cssg_1 = 1.7d0
cssg_1_star = 0.9d0
cssg_2 = 1.05d0
cssg_3 = 0.8d0
cssg_3_star = 0.65d0 
cssg_4 = 0.625d0
cssg_5 = 0.2d0

alpha_1 = 1.0d0 - cssg_1
alpha_1s = - 1.0d0 - cssg_1_star
alpha_3 = 2.0d0/3.0d0 - cssg_3/2.0d0
alpha_3s = cssg_3_star/2.0d0
alpha_4 = 2.0d0*(1.0d0 - cssg_4)
alpha_5 = 2.0d0*(1.0d0 - cssg_5)
!###################################################################
!End of Constants  

do i=1,3
  do j=1,3
    aij(i,j)=0.d0
    sij(i,j)=0.d0
    dsijdt(i,j)=0.d0
    oij(i,j)=0.d0
    delta(i,j)=0.d0
  enddo
  delta(i,i)=1.d0
enddo

!IMPLICIT NONE


turbmod3 = 0

if(turbmod.eq.1)then
	WRITE(*,*) 'Standard k-epsilon Model'
	else if (turbmod.eq.2)then
	WRITE(*,*) 'Reynolds Stress Model'
	else if (turbmod.eq.3)then
	WRITE(*,*) 'Cas implemented in the k-epsilon Model'
	else
	WRITE(*,*) 'Error!!'
endif


!######################  Write to File  ############################

if(turbmod.eq.1)then
OPEN(FILE='Output_keps.dat', unit=11)
REWIND(11)
WRITE(11,*) '        n,            x ,         t,                    s,     Pk,     k,    epsilon '
endif

if(turbmod.eq.2)then
OPEN(FILE='Output_RSM.dat', unit=11)
OPEN(FILE='Output_RSM_stresses.dat', unit=12)
REWIND(11)
REWIND(12)
WRITE(11,*) '        n,            x ,         t,                    s,     Pk,     k,    epsilon '
WRITE(12,*) 'a11,   a22,   a33,   a12,   u1u1,     u2u2,    u3u3,     u1u2 '
  !if (turbmod3.eq.1)then
!	  WRITE(13,*)'Cas'
 ! endif
endif

if(turbmod.eq.3)then
OPEN(FILE='Output_Cas.dat', unit=11)
REWIND(11)
WRITE(11,*) '        n,            x ,         t,                    s,     Pk,     k,    epsilon,  '
WRITE(12,*) 'a11,   a22,   a33,   a12,   u1u1,     u2u2,    u3u3,     u1u2,       cas'
endif

!#######################  End of Write to File #####################


if(casetype.eq.1)then

	WRITE(*,*) 'A Constant Shear Case'

	else if (casetype.eq.2)then
	WRITE(*,*) 'B Oscillating Shear Case'

	else if (casetype.eq.3)then
	WRITE(*,*) 'C Oscillating Strain Case'

	else if (casetype.eq.4)then
	WRITE(*,*) 'D Strain-Relax-Destrain Case'

	else if (casetype.eq.5)then
	WRITE(*,*) 'E Constant Rotation Case'

	else
	WRITE(*,*) 'Error!!'


endif

!WRITE(*,*) 'TIME STEP SIZE:'
!WRITE(*,*) dt
!WRITE(*,*) '**********************************************************'
!WRITE(*,*) '                      CALCULATION'
!WRITE(*,*) 
!Initialise k-epsilon
k_old = k0
e_old = e0
kn = k0

!Initialise RSM
u1u1_n = (2.0d0/3.0d0)*k0
u2u2_n = u1u1_n
u3u3_n = u1u1_n
u1u2_n = 0.0d0
u1u1_old = u1u1_n
u2u2_old = u2u2_n
u3u3_old = u3u3_n
u1u2_old = u1u2_n

!Initialise Cas-k-epsilon
cas_n = 0.3
cas_old = cas_n

n = 0
1001 continue
n = n + 1
if (n.ge.ntmax)then
	stop
endif
!do x=1,ntmax
  t = n * dt
  
!#############################################################################################################
!######################################### Forcing Functions Definiton #######################################
!#############################################################################################################
if(casetype.eq.1)then
!Constant Shear
	s = 1.65
	dsdt = 0.0
	sij(1,2) = 0.5d0*(s)
	sij(2,1) = sij(1,2)
	oij(1,2) = sij(1,2)
	oij(2,1) = -oij(1,2)
        dsijdt(1,2)=0.5d0*dsdt
        dsijdt(2,1)=dsijdt(1,2)

        snorm=abs(s) !my sbar
        sbar =0.5d0*snorm**2
        obar =0.5d0*snorm**2

!Oscillating Shear
	else if (casetype.eq.2)then
	if(n.eq.1)then
		s=smax*cos(omega*t)
  		dsdt=-omega*smax*sin(omega*t)
	else
		s=smax*sin(omega*t)
		dsdt=omega*smax*cos(omega*t)
	endif
	sij(1,2) = 0.5d0*(s+s)
	sij(2,1) = sij(1,2)
	oij(1,2) = sij(1,2)
	oij(2,1) = -oij(1,2)
        dsijdt(1,2)=0.5d0*dsdt
        dsijdt(2,1)=dsijdt(1,2)

        snorm=abs(s) !my sbar
        sbar =0.5d0*snorm**2
        obar =0.5d0*snorm**2
	!Declaring epsilon
	if(n.eq.0)then
		e0 = (sbar*k0)/eta
	endif

!Oscillating Strain
	else if (casetype.eq.3)then
	!s = t
	s = smax*cos(omega*t)/(sin(omega*t)+5.d0/4.d0)
	dsdt = -omega*smax*((sin(omega*t)/(sin(omega*t)+5.d0/4.d0))+(cos(omega*t)/(sin(omega*t)+5.d0/4.d0))**2)
	sij(1,1) = s*(2.d0/3.d0)
	sij(2,2) = s*(-1.d0/3.d0)
	sij(3,3) = sij(2,2)
	sij_total = sij(1,1)+sij(2,2)+sij(3,3)
	if(sij_total.ne.0.d0)then
		WRITE(*,*)'Error! - incompressibility is not achieved'
	endif
  	sij(1,1)=s*(1.d0-1.d0/3.d0)
  	sij(2,2)=s*(0.d0-1.d0/3.d0)
  	sij(3,3)=s*(0.d0-1.d0/3.d0)

 	dsijdt(1,1)=dsdt*(1.d0-1.d0/3.d0)
  	dsijdt(2,2)=dsdt*(0.d0-1.d0/3.d0)
  	dsijdt(3,3)=dsdt*(0.d0-1.d0/3.d0)

  	sbar = (sij(1,1)**2+sij(2,2)**2+sij(3,3)**2)
  	obar = sbar
  	snorm = sqrt(2.d0*sbar)
	!Declaring epsilon
	if(n.eq.0)then
		e0 = (sbar*k0)/eta
	endif


!Strain-relax-destrain
	else if (casetype.eq.4)then
		if ((t.ge.0.d0).and.(t.lt.0.65d0))then
			s = 0.0d0
		else if ((t.ge.0.65d0).and.(t.lt.1.05d0))then
			s = 9.65d0*t - 6.27d0
		else if ((t.ge.1.05d0).and.(t.lt.1.47d0))then
			s = -9.65d0*t + 14.19d0
		else if ((t.ge.1.47d0).and.(t.lt.1.89d0))then
			s = 0.0d0
		else if ((t.ge.1.89d0).and.(t.lt.2.05d0))then
			s = (-21.38*t + 40.41)
			!s = -25.8*t
		else if ((t.ge.2.05d0).and.(t.lt.2.50d0))then
			s = 7.6d0*t - 19.0d0
		endif
		WRITE(*,*) 's', s
		!READ(*,*)
	sij(1,1) = s*(2.d0/3.d0)
	sij(2,2) = s*(-1.d0/3.d0)
	sij(3,3) = sij(2,2)
	sij_total = sij(1,1)+sij(2,2)+sij(3,3)
	if(sij_total.ne.0.d0)then
		WRITE(*,*)'Error! - incompressibility is not achieved'
	endif
  	sij(1,1)=s*(1.d0-1.d0/3.d0)
  	sij(2,2)=s*(0.d0-1.d0/3.d0)
  	sij(3,3)=s*(0.d0-1.d0/3.d0)

 	dsijdt(1,1)=dsdt*(1.d0-1.d0/3.d0)
  	dsijdt(2,2)=dsdt*(0.d0-1.d0/3.d0)
  	dsijdt(3,3)=dsdt*(0.d0-1.d0/3.d0)

  	sbar = (sij(1,1)**2+sij(2,2)**2+sij(3,3)**2)
  	obar = sbar
  	snorm = sqrt(2.d0*sbar)
	!Declaring epsilon
	if(n.eq.0)then
		e0 = (sbar*k0)/eta
	endif


!Constant rotation 
	else if (casetype.eq.5)then
	s = 5
    
	else
	s = 0 

endif

!Printing gradients 

!WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++&
!	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
!WRITE(*,*)'1. FORCING FUNCTION'
!WRITE(*,*)
!WRITE(*,*) 'time step number:', n
!!WRITE(*,*) 'iteration:', x
!WRITE(*,*) 'time step:', t
!WRITE(*,*) 'S', s
!WRITE(*,*) 'dsij/dt'
!WRITE(*,*) dsijdt(1,1), dsijdt(1,2), dsijdt(1,3)
!WRITE(*,*) dsijdt(2,1), dsijdt(2,2), dsijdt(2,3)
!WRITE(*,*) dsijdt(3,1), dsijdt(3,2), dsijdt(3,3)
!WRITE(*,*) 'Sij                                          Oij'
!WRITE(*,*) sij(1,1), sij(1,2), sij(1,3), '|', oij(1,1), oij(1,2), oij(1,3)
!WRITE(*,*) sij(2,1), sij(2,2), sij(2,3), '|', oij(2,1), oij(2,2), oij(2,3)
!WRITE(*,*) sij(3,1), sij(3,2), sij(3,3), '|', oij(3,1), oij(3,2), oij(3,3)
!!WRITE(*,*) '|S|:', sbar
!WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


!WRITE(11,*)n, x, t, s
!==========================================  Turbulence Model  =================================================



x = 0
1002 continue
!1003 continue
if (x.ge.niter)then
	stop
endif

!@@@@@@@@@@ k-epsilon
if(turbmod.eq.1)then

prodk = 2*cmu*((k_old**2)/e_old)*(sbar**2.0d0)
k_new = kn + dt*(prodk - e_old)
e_new = en + dt*(prodk*ce1 - e_old*ce2)*(e_old/k_old)  
  
kdiff = k_new - k_old
ediff = e_new - e_old

!WRITE(*,*) k0, e0
!W!RITE(*,*) k_old, e_old
!WRITE(*,*) kdiff, ediff
!READ(*,*)
endif
!@@@@@@ END: k-epsilons

!@@@@@@@@@ RSM
if(turbmod.eq.2)then
	
prodk = u1u1_old**2 + u2u2_old**2 + u3u3_old**2 + 2.0d0*(u1u2_old**2)

a11_old = (u1u1_old/k_old) - tthd
a22_old = (u2u2_old/k_old) - tthd
a33_old = (u3u3_old/k_old) - tthd
a12_old = (u1u2_old/k_old)

A2 = a11_old**2 + a22_old**2 + a33_old**2 + 2*(a12_old**2)
ABAR = sqrt(A2)
K = k_old
S11 = sij(1,1)
S12 = sij(1,2)
S13 = sij(1,3)
S21 = sij(2,1)
S22 = sij(2,2)
S23 = sij(2,3)
S31 = sij(3,1)
S32 = sij(3,2)
S33 = sij(3,3)
A11 = a11_old
A22 = a22_old
A33 = a33_old
A12 = a12_old
O12 = -S12
E = e_old
PROD = prodk


SRCR11 = -K*(FTHD*S11+2.0*(A11*S11+A12*S12+A12*O12))-A11*(cssg_1*E+cssg_1_star*PROD)+cssg_2*E*(A11**2.0+A12**2.0-OTHD*A2)+&
		S11*K*(cssg_3-cssg_3_star*ABAR)+cssg_4*K*(2.0*(S11*A11+A12*S12)-TTHD*(A11*S11+A22*S22+A33*S33+2.0*A12*S12))+cssg_5*K*(2.0*A12*O12)-TTHD*E
SRCR22 = -K*(FTHD*S22+2.0*(A22*S22+A12*S12+A12*O21))-A22*(cssg_1*E+cssg_1_star*PROD)+cssg_2*E*(A22**2.0+A12**2.0-OTHD*A2)+&
		S22*K*(cssg_3-cssg_3_star*ABAR)+cssg_4*K*(2.0*(S22*A22+A12*S12)-TTHD*(A11*S11+A22*S22+A33*S33+2.0*A12*S12))+cssg_5*K*(2.0*A12*O21)-TTHD*E
SRCR33 = -K*(FTHD*S33+2.0*A33*S33)-A33*(cssg_1*E+cssg_1_star*PROD)+cssg_2*E*(A33**2.0-OTHD*A2)+S33*K*(cssg_3-cssg_3_star*ABAR)+&
		cssg_4*K*(2.0*S33*A33-TTHD*(A11*S11+A22*S22+A33*S33+2.0*A12*S12))-TTHD*E
SRCR12 = -K*(FTHD*S12+S12*(A11+A22)+A12*(S11+S22)+A11*O21+A22*O12)-A12*(cssg_1*E+cssg_1_star*PROD)+cssg_2*E*(A11*A12+A12*A22)+&
		S12*K*(cssg_3-cssg_3_star*ABAR)+cssg_4*K*(S12*(A11+A22)+A12*(S11+S22))+cssg_5*K*(A11*O21+A22*O12)

u1u1_new = u1u1_n + dt*SRCR11
u2u2_new = u2u2_n + dt*SRCR22
u3u3_new = u3u3_n + dt*SRCR33
u1u2_new = u1u2_n + dt*SRCR12

u1u1_diff = u1u1_new - u1u1_old
u2u2_diff = u2u2_new - u2u2_old
u3u3_diff = u3u3_new - u3u3_old
u1u2_diff = u1u2_new - u1u2_old

k_new = 0.5d0*(u1u1_new+u2u2_new+u3u3_new)
e_new = en + dt*(prodk*ce1 - e_old*ce2)*(e_old/k_old)  
ediff = e_new - e_old
!if(turbmod3.eq.1)then
aij_Sij = 0.0d0
	do i=1,3
		do j=1,3
			aij_Sij = aij_Sij + aij(i,j)*sij(i,j)	
		enddo
	enddo
cas_exact = aij_Sij / SBAR
WRITE(*,*)'Cas exact', cas_exact
WRITE(13,*) cas_exact
	!endif
	
WRITE(*,*) u1u1_diff, u2u2_diff, u3u3_diff, u1u2_diff, ediff
!READ(*,*)
endif



!@@@@@@ END: RSM

!@@@@@@@@@@@@@@ Cas Parameter @@@@@@@@@@@@@@@@

if(turbmod.eq.3)then
	
		WRITE(*,*)'Modelled value of cas'
		prodk = 2*cmu*((k_old**2)/e_old)*(sbar**2.0d0)
		k_new = kn + dt*(prodk - e_old)
		e_new = en + dt*(prodk*ce1 - e_old*ce2)
		
		!Modelled cas
		aik_ajk = 0.d0
		do i=1,3
	      do j=1,3!
			do l=1,3
			 ! aij(i,j) = aij(i,j)-xv*tt**2*(xcmu2*(sij(i,l)*oij(l,j)+sij(j,l)*oij(l,i)) &
			  !               -xcmu3*(sij(i,l)*sij(l,j)))
			  aik_ajk = aik_ajk + aij(i,l)*aij(l,j) 
			 enddo
		  enddo
		enddo
		
		a2=0.d0
		saa=0.d0 ! re-added this term \alpha_2 from quadratic part of SSG
		ass=0.d0
		aos=0.d0
		term7=0.d0
		do i=1,3
		  do j=1,3
		    a2=a2+aij(i,j)**2
			! original form:        
			term7=term7 - 1d0/snorm*dsijdt(i,j)*(cas-casa)*2.d0*sij(i,j)/snorm

			! new formulation for term 7: 1/S . dSij/dt . abs( Cas(model) - -Cas(from EVM))    
			!
		!	     term7=term7 - 1d0/snorm*dsijdt(i,j)*abs(cas + etal*0.09d0)&
		!	    *2.d0*sij(i,j)/snorm
			do l=1,3
			      saa=saa+sij(i,j)*aij(i,l)*aij(l,j)
			      ass=ass+sij(i,j)*aij(i,l)*sij(j,l)
			      aos=aos+sij(i,j)*aij(i,l)*oij(j,l)
		    enddo
		  enddo
		enddo


	prod = k*min(0.09d0*k/ep,cas/snorm)*snorm**2
	!    prod = k*cas*snorm

	t1 = alpha_1*cas/tt 
	t2 = alpha_1s*cas**2*snorm 
	t2a = - 1.05d0*saa*ep/(k*snorm) ! term reinstated 
	t3 = alpha_3*snorm 
	t4 = alpha_3s*sqrt(a2)*snorm ! generally I reduce alpha_3s by 0.7x 
	t5 = 1/snorm*alpha_4*ass
	t6 = 1/snorm*alpha_5*aos  
	t7 = + term7
	!
	dcasdt = t1 + t2 + t2a + t3 + t4 + t5 + t6 + t7
	cas = casn + dt*dcasdt
	
	k_new = kn + dt*(prodk - e_old)
	e_new = en + dt*(prodk*ce1 - e_old*ce2)*(e_old/k_old)  
  
	kdiff = k_new - k_old
	ediff = e_new - e_old
	casdiff = cas - dcasdt
	
endif

!@@@@@@@@@@@@@ END: Cas Parameter @@@@@@@@@@@@@

!WRITE(*,*)kdiff,ediff
!Printing turbulent quantities
!write(*,*) 
!WRITE(*,*)'2. TURBULENCE QUANTITIES'
!WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
!WRITE(*,*) 'Initial k and e,  k0:', k0,'e0:',e0
!WRITE(*,*) 'At previous (converged) time step, k:', kn, 'e:', en
!WRITE(*,*) 'k old:', k_old,'e_old:', e_old
!WRITE(*,*) 'Pk:', prodk
!WRITE(*,*) 'k_new:', k_new,'e_new:', e_new
!WRITE(*,*) 'residuals, k:', kdiff,'e:', ediff

!WRITE(*,*)

!if(n.eq.niter)then
!WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++&
!	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
!endif

if (turbmod.eq.1)then
	if ((abs(kdiff).le.residual).and.(abs(ediff).le.residual))then
	  WRITE(*,*)'Converged!'
	  kn = k_new
 	  en = e_new
	  k_old = k_new
	  e_old = e_new
	  WRITE(11,2001)n,x,t,s,prodk, k_new, e_new
	  if (k_new.ge.100000000)then
		  stop
	  endif
	!second star is formt - look for exponential form!
	  GOTO 1001
	else
	  WRITE(*,*)'Not Converged... Carry on calculating'
	  k_old = k_new
	  e_old = e_new
	  GOTO 1002
        endif
endif

if (turbmod.eq.2)then
	if ((abs(u1u1_diff).le.residual).and.(abs(u2u2_diff).le.residual).and.&
	   (abs(u3u3_diff).le.residual).and.(abs(u1u2_diff).le.residual).and.&
	   (abs(ediff).le.residual))then
       WRITE(*,*)'Converged!'
   	   k_new = 0.5d0*(u1u1_new + u2u2_new + u3u3_new)
   	   en = e_new
	   u1u1_n = u1u1_new
	   u2u2_n = u2u2_new
	   u3u3_n = u3u3_new
	   u1u2_n = u1u2_new
   	   WRITE(11,2001)n,x,t,s,prodk, k_new, e_new
	   WRITE(12,2002)a11_old, a22_old, a33_old, a12_old, u1u1_new, u2u2_new, u3u3_new, u1u2_new
   	   GOTO 1001
   	   else
   	   WRITE(*,*)'Not Converged... Carry on calculating'
   	   u1u1_old = u1u1_new
	   u2u2_old = u2u2_new
	   u3u3_old = u3u3_new
	   u1u2_old = u1u2_new
   	   e_old = e_new
   	   GOTO 1002
    endif		   
	   
endif

if (turbmod.eq.3)then
	if ((abs(kdiff).eq.residual).and.(abs(ediff).eq.residual).and.(abs(casdiff).eq.residual))then
		WRITE(*,*)'Converged!'
		k_old = k_new
		kn = k_new
		e_old = e_new
		en = e_new
		cas_n = cas_new
		WRITE(13,2001)n,x,t,s,prodk, k_new, e_new, cas
	 	WRITE(13,2002)a11_old, a22_old, a33_old, a12_old, u1u1_new, u2u2_new, u3u3_new, u1u2_new
		GOTO 1001
	    ELSE
		WRITE(*,*)'Not Converged... Carry on calculating'
		k_old = k_new
		e_old = e_new
		GOTO 1002	
	endif
endif	

!WRITE(11,2001)n,x,t,s,prodk, k_new, e_new
! Integer is I# where # is total including space and sign
! double precision Ew.d, where w is total including space and sign, d is total to right of decimal point 
2001 FORMAT(2I8,5E16.8)
2002 FORMAT(8E16.8)

end program zeroD

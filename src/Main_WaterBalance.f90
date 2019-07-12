! ===================================================================!
!   UBMOD  -- Water balance method of one-dimensional soil water     !
!             movement. Version 1.0.                                 !
!                                                                    !
!   Designed by Wei Mao, Yan Zhu and Jinzhong Yang.                  !
!                                                                    !
!   Cited: Mao W, Yang J, Zhu Y, et al. An efficient soil water      !
!          balance model based on hybrid numerical and statistical   !
!          methods[J]. Journal of hydrology, 2018, 559: 721-735,     !
!          doi: 10.1016/j.jhydrol.2018.02.074.                       !
!                                                                    !
!   Feel free to contact us if you have any question.                !
!       Email: weimao@whu.edu.cn, zyan0701@163.com                   !
!                                                                    !
!                            Last modified: Jul, 2019, by Wei Mao.   !
! ===================================================================!
! ====================================================================
!     Input files:
!     1. Essential files.
!         SELECTOR.IN     The basic input information.
!         uz.IN           The discrete information.
!     2. Optional files.
!         cropdat.dat     The simple crop model.
!         01.wea          Meteorological data.
! ====================================================================
! ====================================================================
!     Output files:
!         
! ====================================================================
!   storage Routing Method~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   *************************end of documentation.	
    PROGRAM WaterBalance
    USE parm
    IMPLICIT NONE
    REAL (KIND=KR) :: t1, t2, Zero
    Zero = 0.0_KR
    
    iof='/sim311/'	!the simulated project file fold

!---Open the interface files. The relpath is used here.
! ====================================================================
!   input files
    OPEN(33,file='Rh1D.in/'//trim(iof)//'/SELECTOR.IN', status='old')
    OPEN(32,file='Rh1D.in/'//trim(iof)//'/uz.IN',       status='old')
!   output files
    OPEN(90,file='Rh1D.out/'//trim(iof)//'/runtime.OUT',  status='unknown') ! Run time.
    OPEN(80,file='Rh1D.out/'//trim(iof)//'/thObs.dat',    status='unknown') ! Node data.
    OPEN(89,file='Rh1D.out/'//trim(iof)//'/balance1d.dat',status='unknown') ! Statistic boundary condition.
    OPEN(81,file='Rh1D.out/'//trim(iof)//'/th.dat',       status='unknown') 
    OPEN(99,file='Rh1D.out/'//trim(iof)//'/error.txt',    status='unknown') ! Error message.
! ====================================================================
      
!-----Begin of the program.
! ====================================================================
!     subroutine about input information.
!     call for basic information. 
    CALL SelectorIn
!     call for node information in 1D.
    CALL UzIn
! ====================================================================

!-----preparation of the calculation
! ====================================================================
!     CPU time.
    CALL CPU_time (t1)
!     The initial water amount in model.
    CALL Balance_Initial
!     Diffusion model.
    CALL Diffusion_Model
!     Call for reference Evaportranspiration and division of E&T.
    CALL Upper_Boundary
! ====================================================================

!-----Begin time loop.
100 CONTINUE
    
    CALL Tcontrol

! ====================================================================
!   Set upper boundary condition.
    CALL SetQ
  
! ====================================================================
!     Four main processes.
!     Firstly, divide infiltration and surface runoff.
!     There is no surface hydrological process by now.
          !if () then
          !    call SCS
          !elseif () then
          !    call Green-Ampt
          !else
          !    qair = Q_infiltration
          !endif

! ====================================================================
!     Secondly, advective movement driven by gravitational potential.
!       "Tipping-bucket" method.
    CALL redistribution

! ====================================================================
!     Thirdly, source/sink term.
!	  open the Files that stored E&T and the rain, and the writen ETa.
    IF(bup >= Zero) CALL SetET

! ====================================================================
!     Last, Diffusive soil water movement driven by matric potential.
    CALL unsatflow

! ====================================================================
!     Output control.
 !	Output the hydraulic head and soil moisture in 1D model.
    CALL Hthuz_out       
!	Call for the water balance in 1D and 3D model.
    CALL BalanceT   
!	call for new time and time step.
    WRITE(*,*)"t=",sngl(t)
            
    IF (bup == 1) THEN
        WRITE(150,'(A18,F10.3,2F10.6)')date,t,Epa,Tra
        CALL DateAdd(date,1,date)
        Epa=Zero
        Tra=Zero
    ENDIF

! ====================================================================
    IF(t-Tend >= Zero) THEN
        GOTO 200
    ELSE
        GOTO 100
    ENDIF
    
200 CALL CPU_time (t2)
    WRITE(90,*)'Real time [sec]',t2-t1
    CLOSE(90)
    CLOSE(80) 
    CLOSE(89) 
    CLOSE(81)
    CLOSE(88)
    CLOSE(98)
    CLOSE(99)
    CLOSE(110) 
    CLOSE(130)
    CLOSE(150) 
    
    STOP
    END PROGRAM WaterBalance

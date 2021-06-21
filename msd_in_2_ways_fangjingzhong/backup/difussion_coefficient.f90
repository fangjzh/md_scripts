PROGRAM diffusion_coefficient
!
! Purpose :
!    To calculate diffusion coefficient of atoms in periodic boundary condition
!
!    Record of revisions:
!      Date         Programer      Description of change
!      ====         =========      =====================
!   2013/05/24    Jingzhong Fang      Original code
!   2014/10/27    Jingzhong Fang      Fix bug of reading inputfile
!                                     formating inputfile

IMPLICIT NONE

! Data dictionary: declare constants
INTEGER, PARAMETER :: MAX_NUM=100               ! Max number of atoms
INTEGER, PARAMETER :: MAX_STEP=1000000          ! Max number of total steps of atomic coordinate data

! Atomic coordinate in angstrom
REAL, DIMENSION(MAX_NUM,MAX_STEP) :: x0,y0,z0   ! Real position of the atoms 
REAL, DIMENSION(MAX_NUM) :: xA,yA,zA            ! Atom's current position in periodic boundary condition 
REAL, DIMENSION(MAX_NUM) :: xB,yB,zB            ! Atom's previous position in periodic boundary condition
REAL :: dx,dy,dz                                 ! Value of xA-xB,yA-yB,zA-zB

! Data dictionary:
INTEGER :: n                                     ! Number of total steps(one set of data corresponding to one step)
INTEGER :: m                                     ! Number of atoms
INTEGER :: step_i                                ! Current step
INTEGER :: atom_i                                ! The numerical symbols of atom "i"
INTEGER :: indiv                                 ! Interval of steps 
INTEGER :: indiv_perfect                         ! The best indiv which has diff_coef_perfect

! Intermediate variable
INTEGER :: status                                ! The status of open file operation
INTEGER :: i
INTEGER :: j

! Data dictionary: Mean square any
REAL*8 :: MSD_single,MSD_sum                     ! Mean square displacement in angstrom squal
REAL*8 :: MSE,MSE_sum,MSE_minimum                ! Mean square error

! Time in femtosecond(1.0E-15s) 
REAL*8 :: delta_t                                ! Interval time between "current and previous position" of atoms

! diffusion coefficient in (angstrom squal/femtosecond)
REAL, DIMENSION(MAX_STEP) :: diff_coef_1,diff_coef_average
REAL*8 :: diff_coef_sum,diff_coef_perfect

! Box size in angstrom
REAL :: box_size_x,box_size_y,box_size_z          ! Box size          
REAL :: half_box_x,half_box_y,half_box_z          ! Half of box size

CHARACTER(len=20) :: file_input                  ! Name of input file(atom position) 
CHARACTER(len=20) :: file_output                 ! Name of output file 
CHARACTER(len=20) :: file_ctrl                   ! Name of control parametre file 
CHARACTER(len=72) :: words                       ! Any words

! Initialize filename
file_ctrl='prog_in.txt'                            ! File unit =  11                           
file_input='atomic.dat'                            ! File unit =  12
file_output='output_info.txt'                      ! File unit =  10

! Initialize data              
MSE_minimum=6     !!! to be adjusted by user
diff_coef_perfect=0.0
indiv_perfect=1

OPEN (UNIT=10, FILE=file_output, STATUS='REPLACE', ACTION='WRITE')

!!!!!!!!!!!!!!!!!!!!
!!!! Get the control parametre of this program start
!!!!!!!!!!!!!!!!!!!!
OPEN (UNIT=11, FILE=file_ctrl, STATUS='OLD', ACTION='READ', IOSTAT=status)
READ (11,99)words                                  ! "Control parametre file"
99 format(A72)
READ (11,99)words                                  ! Informations
WRITE (10,99)words                                 ! Informations                           
READ (11,99)words                                  ! "Box size in angstrom"
READ (11,*)box_size_x,box_size_y,box_size_z  

        half_box_x=box_size_x/2.0    
        half_box_y=box_size_y/2.0
        half_box_z=box_size_z/2.0

WRITE (10,100)box_size_x,box_size_y,box_size_z
100 format('Box size is: ',3F)
 
READ (11,99)words                                  ! "Number of atoms"
READ (11,*)m
WRITE (10,101)m
101 format('Number of atoms is: ',I)

READ (11,99)words                                  ! "Number of total steps"
READ (11,*)n
WRITE (10,102)n
102 format('Number of total steps is: ',I)

READ (11,99)words                                  ! "Interval time"
READ (11,*)delta_t
WRITE (10,103)delta_t
103 format('Delta_t is: ',F8.3,' fs')

CLOSE (UNIT=11, IOSTAT=status)
!!!!!!!!!!!!!!!!!!!!
!!!! Get the control parametre of this program end
!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!
!!!! Convert Atom's position in periodic boundary condition into real position start
!!!!!!!!!!!!!!!!!!!!

OPEN (UNIT=12, FILE=file_input, STATUS='OLD', ACTION='READ', IOSTAT=status)

read_step: DO i=1,n
              READ (12,*)words,step_i
 read_atom: DO j=1,m             
               READ (12,*)xA(j),yA(j),zA(j)         ! Read the coordinates of atoms 

   if_mark02: IF(i==1) THEN
                xB(j)=xA(j)
                yB(j)=yA(j)
                zB(j)=zA(j)
  
                x0(j,i)=xA(j)
                y0(j,i)=yA(j)
                z0(j,i)=zA(j)
              ELSE
                dx=xA(j)-xB(j)
                dy=yA(j)-yB(j)
                dz=zA(j)-zB(j)

                IF(dx.gt.half_box_x)  dx=dx-box_size_x
                IF(dy.gt.half_box_y)  dy=dy-box_size_y
                IF(dz.gt.half_box_z)  dz=dz-box_size_z
 
                IF(dx.lt.-half_box_x) dx=dx+box_size_x
                IF(dy.lt.-half_box_y) dy=dy+box_size_y
                IF(dz.lt.-half_box_z) dz=dz+box_size_z
     
                x0(j,i)=x0(j,i-1)+dx
                y0(j,i)=y0(j,i-1)+dy
                z0(j,i)=z0(j,i-1)+dz
     
                xB(j)=xA(j)
                yB(j)=yA(j)
                zB(j)=zA(j)
              END IF if_mark02
            END DO read_atom
           END DO read_step

CLOSE (UNIT=12, IOSTAT=status)
!!!!!!!!!!!!!!!!!!!!
!!!! Convert Atom's position in periodic boundary condition into real position end
!!!!!!!!!!!!!!!!!!!!


ergodic_step_length: DO indiv=5,n/10
                       diff_coef_sum=0.0
        ergodic_steps: DO i=1,n-indiv
                         MSD_sum=0.0d0                      ! Calculate the diffusion coefficient start
          ergodic_atoms: DO j=1,m
                           MSD_single=(x0(j,i+indiv)-x0(j,i))**2+(y0(j,i+indiv)-y0(j,i))**2+(z0(j,i+indiv)-z0(j,i))**2
                           MSD_sum=MSD_sum+MSD_single
                         END DO ergodic_atoms
                         diff_coef_1(i)=MSD_sum/(REAL(m)*6.0*REAL(indiv)*delta_t) 
                         diff_coef_sum=diff_coef_sum+diff_coef_1(i)
                       END DO ergodic_steps
                       diff_coef_average(indiv)=diff_coef_sum/REAL(n-indiv)
                                                            ! Calculate the diffusion coefficient end

                       MSE_sum=0.0d0                        ! Calculate the mean square error start
            do_mark01: DO i=1,n-indiv           
                         MSE_sum=MSE_sum+(1-diff_coef_1(i)/diff_coef_average(indiv))**2
                       END DO do_mark01
                       MSE=MSE_sum/REAL(n-indiv)                ! Calculate the mean square error end

            if_mark03: IF(MSE_minimum>MSE) THEN             ! Get the most stable diffusion coefficient start
                         MSE_minimum=MSE
                         diff_coef_perfect=diff_coef_average(indiv)
                         indiv_perfect =indiv
                       END IF if_mark03                     ! Get the most stable diffusion coefficient end
                     END DO ergodic_step_length
 
WRITE (10,104)diff_coef_perfect
104   format('Diffusion coefficient is: ',ES14.5)
WRITE (10,105)MSE_minimum
105   format('MSE: ',ES14.5)
WRITE (10,106)indiv_perfect
106   format('individual: ',I)

END PROGRAM diffusion_coefficient

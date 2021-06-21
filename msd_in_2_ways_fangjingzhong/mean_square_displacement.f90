PROGRAM mean_square_displacement
!
! Purpose :
!    To calculate mean square displacement of atoms in periodic boundary condition
! Note :
!    The method of this program is refering to the work of Tiefeng Wu @ Beihang University
!    Record of revisions:
!      Date         Programer      Description of change
!      ====         =========      =====================
!   2013/05/29    Jingzhong Fang      Original code
!   2014/10/27    Jingzhong Fang      Fix bug of reading inputfile while number of atome>1
!                                     Fix bug of msdsum(divided by m) while number of atome>1
!                                     formating inputfile
!                                     allocatable matrix.
IMPLICIT NONE

! Data dictionary: declare constants
!INTEGER, PARAMETER :: MAX_NUM=100               ! Max number of atoms
!INTEGER, PARAMETER :: MAX_STEP=1000000          ! Max number of total steps of atomic coordinate data

! Atomic coordinate in angstrom
!REAL, DIMENSION(MAX_NUM,MAX_STEP) :: x0,y0,z0   ! 
REAL,DIMENSION(:,:),allocatable :: x0,y0,z0      ! Real position of the atoms 
!REAL, DIMENSION(MAX_NUM) :: xA,yA,zA            !
REAL,DIMENSION(:),allocatable   :: xA,yA,zA      ! Atom's current position in periodic boundary condition 
!REAL, DIMENSION(MAX_NUM) :: xB,yB,zB            ! 
REAL,DIMENSION(:),allocatable   :: xB,yB,zB      ! Atom's previous position in periodic boundary condition
REAL :: dx,dy,dz                                 ! Value of xA-xB,yA-yB,zA-zB

! Data dictionary:
INTEGER :: n                                     ! Number of total steps(one set of data corresponding to one step)
INTEGER :: m                                     ! Number of atoms
INTEGER :: step_i                                ! Current step
INTEGER :: atom_i                                ! The numerical symbols of atom "i"
INTEGER :: indiv                                 ! Interval of steps 

! Intermediate variable
INTEGER :: status,istat                                ! The status of open file operation
INTEGER :: i
INTEGER :: j

! Data dictionary: Mean square any
!REAL, DIMENSION(MAX_STEP) :: time_1,MSD_1        ! 
REAL,DIMENSION(:),allocatable  :: time_1,MSD_1        ! MSD_1 (Average MSD) within time_1
REAL*8 :: MSD_single,MSD_sum                     ! Mean square displacement in angstrom squal
REAL*8 :: MSE,MSE_sum,MSE_minimum,MSD_sum_sum    ! Mean square error

! Statistics for linear fit of MSD data
REAL*8 :: sum_x=0.0,sum_x2=0.0,sum_y=0.0,sum_xy=0.0
REAL*8 :: x_average,y_average,y_int,slope
INTEGER :: msd_number

! diffusion coefficient in (angstrom squal/femtosecond)
REAL*8 :: diff_coef

! Time in femtosecond(1.0E-15s) 
REAL*8 :: delta_t                                ! Interval time between "current and previous position" of atoms

! Box size in angstrom
REAL :: box_size_x,box_size_y,box_size_z          ! Box size          
REAL :: half_box_x,half_box_y,half_box_z          ! Half of box size

CHARACTER(len=20) :: file_input                  ! Name of input file(atom position) 
CHARACTER(len=20) :: file_output                 ! Name of output file 
CHARACTER(len=20) :: file_ctrl                   ! Name of control parametre file 
CHARACTER(len=20) :: file_msd                    ! Name of MSD file
CHARACTER(len=72) :: words                       ! Any words

! Initialize filename
file_output='output_info.txt'                      ! File unit =  10
file_ctrl='prog_in.txt'                            ! File unit =  11                           
file_input='atomic.dat'                            ! File unit =  12
file_msd='msd_out.dat'                                 ! File unit =  13
!CALL getarg(1,file_input)
 

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
100 format('Box size is: ',3F18.7)
 
READ (11,99)words                                  ! "Number of atoms"
READ (11,*)m
WRITE (10,101)m
101 format('Number of atoms is: ',I8)

READ (11,99)words                                  ! "Number of total steps"
READ (11,*)n
WRITE (10,102)n
102 format('Number of total steps is: ',I8)

READ (11,99)words                                  ! "Interval time"
READ (11,*)delta_t
WRITE (10,103)delta_t
103 format('Delta_t is: ',F8.3,' fs')

CLOSE (UNIT=11, IOSTAT=status)

allocate (x0(m,n),y0(m,n),z0(m,n),stat=istat)
allocate (xA(m),yA(m),zA(m),stat=istat)
allocate (xB(m),yB(m),zB(m),stat=istat)
allocate (time_1(n),MSD_1(n),stat=istat)

!!!!!!!!!!!!!!!!!!!!
!!!! Get the control parametre of this program end
!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!
!!!! Convert Atom's position in periodic boundary condition into real position start
!!!!!!!!!!!!!!!!!!!!

OPEN (UNIT=12, FILE=file_input, STATUS='OLD', ACTION='READ', IOSTAT=status)
!READ (12,99)words
!READ (12,99)words
read_step: DO i=1,n
            READ (12,*)words,step_i
!write(*,*)'ni da ye',words,step_i
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

OPEN (UNIT=13, FILE=file_msd, STATUS='REPLACE', ACTION='WRITE')
WRITE(13,*)'time_fs    msd_angtrom_square'
                                                             ! Calculate the mean square displacement start
msd_number=0
ergodic_step_length: DO indiv=5,n/10
                       MSD_sum_sum=0.0
        ergodic_steps: DO i=1,n-indiv
                         MSD_sum=0.0d0                     
          ergodic_atoms: DO j=1,m
                           MSD_single=(x0(j,i+indiv)-x0(j,i))**2+(y0(j,i+indiv)-y0(j,i))**2+(z0(j,i+indiv)-z0(j,i))**2
                           MSD_sum=MSD_sum+MSD_single
                         END DO ergodic_atoms
                           MSD_sum_sum=MSD_sum_sum+MSD_sum
                       END DO ergodic_steps
                       MSD_1(indiv)=MSD_sum_sum/REAL(n-indiv)/m  !divided by number of atoms
                       time_1(indiv)=delta_t*indiv
                       WRITE (13,104)time_1(indiv),MSD_1(indiv)
104 FORMAT(F12.2,' ',F12.5)
                                                             ! Calculate the mean square displacement end
                                                             ! Statistics for linear fit of MSD start
                       sum_x=sum_x+time_1(indiv)
                       sum_x2=sum_x2+time_1(indiv)**2
                       sum_y=sum_y+MSD_1(indiv)
                       sum_xy=sum_xy+time_1(indiv)*MSD_1(indiv)
                       msd_number=msd_number+1                       
                                                            ! Statistics for linear fit of MSD end
                     END DO ergodic_step_length

CLOSE (UNIT=13, IOSTAT=status)
               
x_average=sum_x/REAL(msd_number)
y_average=sum_y/REAL(msd_number)
slope=(sum_xy-sum_x*y_average)/(sum_x2-sum_x*x_average)
y_int=y_average-slope*x_average

diff_coef=slope/6.0
WRITE(10,105) slope ,y_int,diff_coef
105 FORMAT (' ','Regression coefficients for the least-squares line:',&
         /,1X,'slope     (m)        =  ',ES14.5,&
         /,1X,'Intercept (b)        =  ',F12.3,&
         /,1X,'diffusion coefficient=  ',ES14.5,'  angstrom squal/femtosecond')

CLOSE (UNIT=10, IOSTAT=status)

 
END PROGRAM mean_square_displacement

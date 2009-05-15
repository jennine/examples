!> \file
!> $Id: AnalyticLaplaceExample.f90 20 2007-05-28 20:22:52Z cpb $
!> \author Chris Bradley
!> \brief This is an example program to solve an Analytic Laplace equation using openCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is openCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> Main program
PROGRAM ANALYTICLAPLACEEXAMPLE

  USE ANALYTIC_ANALYSIS_ROUTINES
  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CMISS
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE EQUATIONS_SET_ROUTINES
  USE FIELD_ROUTINES
  USE FIELD_IO_ROUTINES
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LAPLACE_EQUATIONS_ROUTINES
  USE LISTS
  USE MESH_ROUTINES
  USE MPI
  USE PROBLEM_CONSTANTS
  USE PROBLEM_ROUTINES
  USE REGION_ROUTINES
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(DP), PARAMETER :: ORIGIN(2)=(/-PI/2, -PI/2/)
  REAL(DP), PARAMETER :: HEIGHT=PI
  REAL(DP), PARAMETER :: WIDTH=PI
  REAL(DP), PARAMETER :: LENGTH=PI

  !Program types

  !Program variables

  INTEGER(INTG) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,INTERPOLATION_SPECIFICATIONS
  INTEGER(INTG) :: NUMBER_OF_DOMAINS

  INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES
  INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER
  INTEGER(INTG) :: MPI_IERROR

  INTEGER(INTG) :: EQUATIONS_SET_INDEX

  TYPE(BASIS_TYPE), POINTER :: BASIS
  TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
  TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
  TYPE(MESH_TYPE), POINTER :: MESH
  TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
  TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
  TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
  TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,GEOMETRIC_FIELD,DEPENDENT_FIELD
  TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
  TYPE(REGION_TYPE), POINTER :: REGION
  TYPE(SOLVER_TYPE), POINTER :: SOLVER
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS

  LOGICAL :: EXPORT_FIELD,IMPORT_FIELD
  TYPE(VARYING_STRING) :: FILE,METHOD
  CHARACTER *100 BUFFER

  REAL(SP) :: START_USER_TIME(1),STOP_USER_TIME(1),START_SYSTEM_TIME(1),STOP_SYSTEM_TIME(1)

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables

  INTEGER(INTG) :: ERR
  TYPE(VARYING_STRING) :: ERROR

  INTEGER(INTG) :: DIAG_LEVEL_LIST(5)
  CHARACTER(LEN=MAXSTRLEN) :: DIAG_ROUTINE_LIST(1),TIMING_ROUTINE_LIST(1)

#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise cmiss
  CALL CMISS_INITIALISE(ERR,ERROR,*999)

  !Set all diganostic levels on for testing
  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5
  !DIAG_ROUTINE_LIST(1)="FINISH_CREATE_FIELD"
  !CALL DIAGNOSTICS_SET_ON(ALL_DIAG_TYPE,DIAG_LEVEL_LIST,"LaplaceExample",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
  !CALL DIAGNOSTICS_SET_ON(ALL_DIAG_TYPE,DIAG_LEVEL_LIST,"",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
  !CALL DIAGNOSTICS_SET_ON(IN_DIAG_TYPE,DIAG_LEVEL_LIST,"LaplaceExample",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
  !CALL DIAGNOSTICS_SET_ON(IN_DIAG_TYPE,DIAG_LEVEL_LIST,"",DIAG_ROUTINE_LIST,ERR,ERROR,*999)

  TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !Put all output to a file.
  CALL OUTPUT_SET_ON("AnalyticLaplace",ERR,ERROR,*999)

  !Calculate the start times
  CALL CPU_TIMER(USER_CPU,START_USER_TIME,ERR,ERROR,*999)
  CALL CPU_TIMER(SYSTEM_CPU,START_SYSTEM_TIME,ERR,ERROR,*999)

  !Get the number of computational nodes
  NUMBER_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999
  !Get my computational node number
  MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999

!  !Read in the number of elements in the X & Y directions, and the number of partitions on the master node (number 0)
!  IF(MY_COMPUTATIONAL_NODE_NUMBER==0) THEN
!    WRITE(*,'("Enter the number of elements in the X direction :")')
!    READ(*,*) NUMBER_GLOBAL_X_ELEMENTS
!    WRITE(*,'("Enter the number of elements in the Y direction :")')
!    READ(*,*) NUMBER_GLOBAL_Y_ELEMENTS
!    WRITE(*,'("Enter the number of elements in the Z direction :")')
!    READ(*,*) NUMBER_GLOBAL_Z_ELEMENTS
!    WRITE(*,'("Enter the number of domains :")')
!    READ(*,*) NUMBER_OF_DOMAINS
!  ENDIF

  IF(MY_COMPUTATIONAL_NODE_NUMBER==0) THEN
    IF(COMMAND_ARGUMENT_COUNT()==5) THEN
      !GET THE PARAMETERS FROM THE COMMAND LINE ARGUMENT
      CALL GET_COMMAND_ARGUMENT(1,BUFFER) 
      READ(BUFFER,*) NUMBER_GLOBAL_X_ELEMENTS
      CALL GET_COMMAND_ARGUMENT(2,BUFFER)
      READ(BUFFER,*) NUMBER_GLOBAL_Y_ELEMENTS
      CALL GET_COMMAND_ARGUMENT(3,BUFFER)
      READ(BUFFER,*) NUMBER_GLOBAL_Z_ELEMENTS
      CALL GET_COMMAND_ARGUMENT(4,BUFFER)
      READ(BUFFER,*) NUMBER_OF_DOMAINS
      CALL GET_COMMAND_ARGUMENT(5,BUFFER)
      READ(BUFFER,*) INTERPOLATION_SPECIFICATIONS
    ELSE
!TODO more detailed error message
      !CALL FLAG_ERROR("Incorrect number of argements.",ERR,ERROR,*999)
      NUMBER_GLOBAL_X_ELEMENTS=3
      NUMBER_GLOBAL_Y_ELEMENTS=3
      NUMBER_GLOBAL_Z_ELEMENTS=0
      NUMBER_OF_DOMAINS=1
      INTERPOLATION_SPECIFICATIONS=1
    ENDIF
  ENDIF

  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(INTERPOLATION_SPECIFICATIONS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"COMPUTATIONAL ENVIRONMENT:",ERR,ERROR,*999)
  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Total number of computaional nodes = ",NUMBER_COMPUTATIONAL_NODES, &
    & ERR,ERROR,*999)
  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  My computational node number = ",MY_COMPUTATIONAL_NODE_NUMBER,ERR,ERROR,*999)

  !Start the creation of a new RC coordinate system
  NULLIFY(COORDINATE_SYSTEM)
  CALL COORDINATE_SYSTEM_CREATE_START(1,COORDINATE_SYSTEM,ERR,ERROR,*999)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,2,ERR,ERROR,*999)
  ELSE
    !Set the coordinate system to be 3D
    CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,3,ERR,ERROR,*999)
  ENDIF
  !Finish the creation of the coordinate system
  CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*999)

  !Start the creation of the region
  NULLIFY(REGION)
  CALL REGION_CREATE_START(1,REGION,ERR,ERROR,*999)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*999)
  !Finish the creation of the region
  CALL REGION_CREATE_FINISH(REGION,ERR,ERROR,*999)

  FILE="LaplaceExample"
  METHOD="Fortran"
  IMPORT_FIELD=.FALSE.
  IF(IMPORT_FIELD) THEN
     CALL FIELD_IO_FILEDS_IMPORT(FILE, METHOD, REGION, MESH, 1, DECOMPOSITION, 1,  &
      & DECOMPOSITION_CALCULATED_TYPE, FIELD_U_VARIABLE_TYPE, FIELD_ARITHMETIC_MEAN_SCALING, ERR, ERROR, *999)
  ELSE
    !Start the creation of a basis (default is trilinear lagrange)
    NULLIFY(BASIS)
    CALL BASIS_CREATE_START(1,BASIS,ERR,ERROR,*999)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the basis to be a bilinear Lagrange basis
      CALL BASIS_NUMBER_OF_XI_SET(BASIS,2,ERR,ERROR,*999)
      ! Using Cubic Hermite.
      CALL BASIS_INTERPOLATION_XI_SET(BASIS,(/INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS/),ERR,ERROR,*999)
    ELSE
      !Set the basis to be a trilinear Lagrange basis
      CALL BASIS_NUMBER_OF_XI_SET(BASIS,3,ERR,ERROR,*999)
      ! Using Cubic Hermite.
      CALL BASIS_INTERPOLATION_XI_SET(BASIS,(/INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS, &
          & INTERPOLATION_SPECIFICATIONS/),ERR,ERROR,*999)
    ENDIF
    !Finish the creation of the basis
    CALL BASIS_CREATE_FINISH(BASIS,ERR,ERROR,*999)

    !Start the creation of a generated mesh in the region
    NULLIFY(GENERATED_MESH)
    NULLIFY(MESH)
    CALL GENERATED_MESH_CREATE_START(1,REGION,GENERATED_MESH,ERR,ERROR,*999)
    !Set up a regular 100x100 mesh
    CALL GENERATED_MESH_TYPE_SET(GENERATED_MESH,1, &
        & ERR,ERROR,*999)
    CALL GENERATED_MESH_BASIS_SET(GENERATED_MESH,BASIS,ERR,ERROR,*999)

    !Define the mesh on the region
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      CALL GENERATED_MESH_EXTENT_SET(GENERATED_MESH,(/WIDTH,HEIGHT/),ERR,ERROR,*999)
      CALL GENERATED_MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/), &
        & ERR,ERROR,*999)
      CALL GENERATED_MESH_ORIGIN_SET(GENERATED_MESH,ORIGIN,ERR,ERROR,*999)
    ELSE
      CALL GENERATED_MESH_EXTENT_SET(GENERATED_MESH,(/WIDTH,HEIGHT,LENGTH/),ERR,ERROR,*999)
      CALL GENERATED_MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH,(/NUMBER_GLOBAL_X_ELEMENTS, &
        & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS/), ERR,ERROR,*999)
    ENDIF

    !Finish the creation of a generated mesh in the region
    CALL GENERATED_MESH_CREATE_FINISH(GENERATED_MESH,1,MESH,ERR,ERROR,*999)

    !Create a decomposition
    NULLIFY(DECOMPOSITION)
    CALL DECOMPOSITION_CREATE_START(1,MESH,DECOMPOSITION,ERR,ERROR,*999)
    !Set the decomposition to be a general decomposition with the specified number of domains
    CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)
    CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*999)
    CALL DECOMPOSITION_CREATE_FINISH(MESH,DECOMPOSITION,ERR,ERROR,*999)

    !Start to create a default (geometric) field on the region
    NULLIFY(GEOMETRIC_FIELD)
    CALL FIELD_CREATE_START(1,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)
    !Set the decomposition to use
    CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*999)
    !Set the domain to be used by the field components
    !NB these are needed now as the default mesh component number is 1
    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1,1,ERR,ERROR,*999)
    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,2,1,ERR,ERROR,*999)
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,3,1,ERR,ERROR,*999)
    ENDIF
    !Finish creating the field
    CALL FIELD_CREATE_FINISH(GEOMETRIC_FIELD,ERR,ERROR,*999)

    !Update the geometric field parameters
    CALL GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE(GEOMETRIC_FIELD,GENERATED_MESH,ERR,ERROR,*999)

  ENDIF

  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD)) GEOMETRIC_FIELD=>REGION%FIELDS%FIELDS(1)%PTR

  !Create the equations_set
  NULLIFY(EQUATIONS_SET)
  CALL EQUATIONS_SET_CREATE_START(1,REGION,GEOMETRIC_FIELD,EQUATIONS_SET,ERR,ERROR,*999)
  !Set the equations set to be a standard Laplace problem
  CALL EQUATIONS_SET_SPECIFICATION_SET(EQUATIONS_SET,EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_LAPLACE_EQUATION_TYPE, &
    & EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE,ERR,ERROR,*999)
  !Finish creating the equations set
  CALL EQUATIONS_SET_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)

  !Create the equations set analytic field variables
  NULLIFY(DEPENDENT_FIELD)
  CALL EQUATIONS_SET_DEPENDENT_CREATE_START(EQUATIONS_SET,2,DEPENDENT_FIELD,ERR,ERROR,*999)
  !Finish the equations set dependent field variables
  CALL EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)

  !Create the equations set analytic field variables
  NULLIFY(ANALYTIC_FIELD)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL EQUATIONS_SET_ANALYTIC_CREATE_START(EQUATIONS_SET,EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2,3,ANALYTIC_FIELD, &
      & ERR,ERROR,*999)
  ELSE
    CALL EQUATIONS_SET_ANALYTIC_CREATE_START(EQUATIONS_SET,EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2,3,ANALYTIC_FIELD, &
      & ERR,ERROR,*999)
  ENDIF
  !Finish the equations set dependent field variables
  CALL EQUATIONS_SET_ANALYTIC_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)

  !Create the equations set equations
  NULLIFY(EQUATIONS)
  CALL EQUATIONS_SET_EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
  !Set the equations matrices sparsity type
  CALL EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,EQUATIONS_SPARSE_MATRICES,ERR,ERROR,*999)
  !CALL EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,EQUATIONS_FULL_MATRICES,ERR,ERROR,*999)
  !Set the equations set output
  !CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_TIMING_OUTPUT,ERR,ERROR,*999)
  CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_MATRIX_OUTPUT,ERR,ERROR,*999)
  !CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_ELEMENT_MATRIX_OUTPUT,ERR,ERROR,*999)

  CALL EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)

  !Set up the boundary conditions as per the analytic solution
  CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC(EQUATIONS_SET,ERR,ERROR,*999)
  
  !Create the problem
  NULLIFY(PROBLEM)
  CALL PROBLEM_CREATE_START(1,PROBLEM,ERR,ERROR,*999)
  !Set the problem to be a standard Laplace problem
  CALL PROBLEM_SPECIFICATION_SET(PROBLEM,PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_LAPLACE_EQUATION_TYPE, &
    & PROBLEM_STANDARD_LAPLACE_SUBTYPE,ERR,ERROR,*999)
  !Finish creating the problem
  CALL PROBLEM_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !Create the problem control loop
  CALL PROBLEM_CONTROL_LOOP_CREATE_START(PROBLEM,ERR,ERROR,*999)
  !Finish creating the problem control
  CALL PROBLEM_CONTROL_LOOP_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !Start the creation of the problem solver
  NULLIFY(SOLVER)
  CALL PROBLEM_SOLVERS_CREATE_START(PROBLEM,ERR,ERROR,*999)
  CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,1,SOLVER,ERR,ERROR,*999)
  !CALL SOLVER_LINEAR_TYPE_SET(SOLVER,SOLVER_LINEAR_DIRECT_SOLVE_TYPE,ERR,ERROR,*999)
  !CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_NO_OUTPUT,ERR,ERROR,*999)
  !CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_TIMING_OUTPUT,ERR,ERROR,*999)
  !CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_SOLVER_OUTPUT,ERR,ERROR,*999)
  CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_MATRIX_OUTPUT,ERR,ERROR,*999)
  !Finish the creation of the problem solver
  CALL PROBLEM_SOLVERS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !Create the problem solver equations
  NULLIFY(SOLVER)
  NULLIFY(SOLVER_EQUATIONS)
  CALL PROBLEM_SOLVER_EQUATIONS_CREATE_START(PROBLEM,ERR,ERROR,*999)
  CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,1,SOLVER,ERR,ERROR,*999)
  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
  !CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER,SOLVER_FULL_MATRICES,ERR,ERROR,*999)
  !Add in the equations set
  CALL SOLVER_EQUATIONS_EQUATIONS_SET_ADD(SOLVER_EQUATIONS,EQUATIONS_SET,EQUATIONS_SET_INDEX,ERR,ERROR,*999)
  !Finish the problem solver equations
  CALL PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !Solve the problem
  CALL PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*999)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    METHOD="FORTRAN"
    !CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)
    !CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)
    !CALL ANALYTIC_ANALYSIS_EXPORT(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, FILE, METHOD, ERR,ERROR,*999)
  ENDIF
  CALL ANALYTIC_ANALYSIS_OUTPUT(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"",ERR,ERROR,*999)

  !Output timing summary
  !CALL TIMING_SUMMARY_OUTPUT(ERR,ERROR,*999)

  !Calculate the stop times and write out the elapsed user and system times
  CALL CPU_TIMER(USER_CPU,STOP_USER_TIME,ERR,ERROR,*999)
  CALL CPU_TIMER(SYSTEM_CPU,STOP_SYSTEM_TIME,ERR,ERROR,*999)

  CALL WRITE_STRING_TWO_VALUE(GENERAL_OUTPUT_TYPE,"User time = ",STOP_USER_TIME(1)-START_USER_TIME(1),", System time = ", &
    & STOP_SYSTEM_TIME(1)-START_SYSTEM_TIME(1),ERR,ERROR,*999)

  CALL OUTPUT_SET_OFF(ERR,ERROR,*999)

  CALL CMISS_FINALISE(ERR,ERROR,*999)

  WRITE(*,'(A)') "Program successfully completed."

  STOP
999 CALL CMISS_WRITE_ERROR(ERR,ERROR)
  STOP 1

END PROGRAM ANALYTICLAPLACEEXAMPLE

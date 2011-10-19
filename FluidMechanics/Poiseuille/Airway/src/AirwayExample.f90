!> \file
!> \author Jennine Mitchell
!> \brief This is an example program to solve Poiseulle Flow in a 1D model of
!> \ the airway
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
!> The Original Code is OpenCMISS
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

!> \example FluidMechanics/Poiseuille/Static/src/StaticExample.f90
!! Example program to solve a static Poiseuille equation using openCMISS calls.
!<

!> Main program
PROGRAM AIRWAYPOISEUILLEEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Program parameters
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7 !geometry
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9 !viscosity/density
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumber=11 !radius could carry gravity
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldViscosityUserNumber=15

  !Program variables
  INTEGER(CMISSIntg) :: INTERPOLATION_TYPE
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  REAL(CMISSDP) :: POSITION(3),PIPE_LENGTH,VISCOSITY
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT

  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber,FirstNodeDomain,LastNodeDomain

  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfComponents=1
  INTEGER(CmissIntg), PARAMETER :: MeshComponentNumber=1



 ! LOGICAL :: EXPORT_FIELD

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,IndependentField
  TYPE(CMISSFieldType) :: MaterialsField, MaterialsFieldViscosity
  TYPE(CMISSFieldsType) :: Fields
!  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSMeshElementsType) :: MeshElements,MeshElements2
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSFieldType) :: EquationsSetField

  !Generic CMISS variables

  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber

!  !Equations sets
!  TYPE(CMISSEquationsSetType) :: EquationsSetPoiseuille
!  !Equations
!  TYPE(CMISSEquationsType) :: EquationsPoiseuille

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG

  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Get input arguments
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS >= 2) THEN
    !If we have enough arguments then use the first four for setting up the problem. The subsequent arguments may be used to
    !pass flags to, say, PETSc.
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_X_ELEMENTS
    IF(NUMBER_GLOBAL_X_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION_TYPE
    IF(INTERPOLATION_TYPE<=0) CALL HANDLE_ERROR("Invalid Interpolation specification.")
  ELSE
    !If there are not enough arguments default the problem specification
    NUMBER_GLOBAL_X_ELEMENTS=2
    INTERPOLATION_TYPE=1
  ENDIF

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap all errors
  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  !Output to a file
  CALL CMISSOutputSetOn("Poiseuille_Airway",Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,3,Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 3D RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

   !Start the creation of a basis, Lagrange or Simplex options ? Linear Lagrange as defualt for constant field interpolation -or do I do something else
   CALL CMISSBasisTypeInitialise(Basis,Err)
   CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
   CALL CMISSBasisTypeSet(Basis,CMISSBasisLagrangeHermiteTPType,Err)
   CALL CMISSBasisNumberOfXiSet(Basis,1,Err)
  !Set the basis xi interpolation and number of Gauss points
   CALL CMISSBasisInterpolationXiSet(Basis,[CMISSBasisLinearLagrangeInterpolation],Err)
   CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis,[2],Err)
   !Finish the creation of the basis
   CALL CMISSBasisCreateFinish(Basis,Err)

  !Create a mesh
  !The mesh will consist of twelve nodes.
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSNodesCreateStart(Region,12,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err) 


  CALL CMISSMeshTypeInitialise(Mesh,Err)
  !create a 1D airway airway mesh
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,1,Mesh,Err)
  CALL CMISSMeshNumberOfElementsSet(Mesh,11,Err)
  CALL CMISSMeshNumberOfComponentsSet(Mesh,2,Err)!add a second mesh component
  CALL CMISSMeshElementsTypeInitialise(MeshElements,Err)
  !for the node based fields
  CALL CMISSMeshElementsCreateStart(Mesh,1,Basis,MeshElements,Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,1,[1,2],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,2,[2,3],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,3,[2,4],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,4,[3,5],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,5,[3,6],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,6,[4,7],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,7,[4,8],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,8,[8,9],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,9,[8,10],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,10,[10,11],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements,11,[10,12],Err)
  CALL CMISSMeshElementsBasisSet(MeshElements,11,Basis,Err)
  CALL CMISSMeshElementsCreateFinish(MeshElements,Err)
  !for the element based fields
  CALL CMISSMeshElementsTypeInitialise(MeshElements2,Err)
  CALL CMISSMeshElementsCreateStart(Mesh,2,Basis,MeshElements2,Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,1,[1,2],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,2,[2,3],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,3,[2,4],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,4,[3,5],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,5,[3,6],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,6,[4,7],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,7,[4,8],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,8,[8,9],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,9,[8,10],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,10,[10,11],Err)
  CALL CMISSMeshElementsNodesSet(MeshElements2,11,[10,12],Err)
  CALL CMISSMeshElementsBasisSet(MeshElements2,11,Basis,Err)
  CALL CMISSMeshElementsCreateFinish(MeshElements2,Err)
  CALL CMISSMeshCreateFinish(Mesh,Err)
!stop-good till here

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a geometric field using unit scaling on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the label
  CALL CMISSFieldVariableLabelSet(GeometricField,CMISSFieldUVariableType,"coordinates",Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the field Scaling
  CALL CMISSFieldScalingTypeSet(GeometricField,CMISSFieldUnitScaling,Err)
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)
 
  !Set Geometric x,y,z values, for 12 nodes
  !Node 1, note only 1 version(location),scalar value only, node 1,3 components (x,y,z)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,2,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,1,3,0.1_CMISSDP,Err)
  !Node 2  note only 1 version(location),scalar value only, node 2,3 components (x,y,z)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,1,-80.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,2,-0.1_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,2,3,0.0_CMISSDP,Err)
  !Node 3
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       &  1,1,3,1,-95.70815_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       &  1,1,3,2,21.06008_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,3,3,-12.71314_CMISSDP,Err)
  !Node 4
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,4,1,-103.1027_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,4,2,-21.41029_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,4,3,15.40255_CMISSDP,Err)
  !Node 5
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,5,1,-93.35438_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,5,2,41.83544_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,5,3,-16.7010_CMISSDP,Err)
  !Node 6
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,6,1,-121.4918_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,6,2,27.14352_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,6,3,-25.33416_CMISSDP,Err)
  !Node 7
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,7,1,-99.60918_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,7,2,-30.53896_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,7,3,27.75301_CMISSDP,Err)
  !Node 8
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,8,1,-128.4728_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,8,2,-29.00407_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,8,3,16.42376_CMISSDP,Err)
  !Node 9
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,9,1,-128.4295_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,9,2,-39.54292_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,9,3,13.58564_CMISSDP,Err)
  !Node 10
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,10,1,-149.4737_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,10,2,-27.74580_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,10,3,12.63901_CMISSDP,Err)
  !Node 11
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,11,1,-154.3792_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,11,2,-31.4978_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,11,3,10.6421_CMISSDP,Err)
  !Node 12
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,12,1,-153.8843_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,12,2,-23.13308_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,&
       & 1,1,12,3,11.56305_CMISSDP,Err)

  !Geometric field update
  CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)


  !Create an independent field with constant interpolation to hold the radius values 
  CALL CMISSFieldTypeInitialise(IndependentField,Err)
  CALL CMISSFieldCreateStart(IndependentFieldUserNumber,Region,IndependentField,Err)
  CALL CMISSFieldTypeSet(IndependentField,CMISSFieldIndependentType,Err)
  CALL CMISSFieldMeshDecompositionSet(IndependentField,Decomposition,Err) 
  CALL CMISSFieldVariableLabelSet(IndependentField,CMISSFieldUVariableType,"radius",Err)
  CALL CMISSFieldNumberOfVariablesSet(IndependentField,FieldIndependentNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(IndependentField,CMISSFieldUVariableType,1,Err) 
!  CALL CMISSFieldComponentMeshComponentSet(IndependentField,CMISSFieldIndependentType,1,MeshComponentNumber,Err)
  CALL CMISSFieldComponentInterpolationSet(IndependentField,CMISSFieldUVariableType,1,&
	& CMISSFieldElementBasedInterpolation,Err)
  ! CALL CMISSFieldScalingTypeSet(IndependentField,CMISSFieldUnitScaling,Err)
  CALL CMISSFieldCreateFinish(IndependentField,Err)
  !set the radius values
  !element 1
  CALL CMISSFieldParameterSetUpdateElement(IndependentField,CMISSFieldIndependentType,&
	& CMISSFieldValuesSetType,1,1,6.0_CMISSDP,Err) 
  !element 2
  CALL CMISSFieldParameterSetUpdateElement(IndependentField,CMISSFieldIndependentType,&
	& CMISSFieldValuesSetType,2,1,4.0_CMISSDP,Err)
  !element 3
  CALL CMISSFieldParameterSetUpdateElement(IndependentField,CMISSFieldIndependentType,&
	& CMISSFieldValuesSetType,3,1,5.0_CMISSDP,Err)
  !element 4
  CALL CMISSFieldParameterSetUpdateElement(IndependentField,CMISSFieldIndependentType,&
	& CMISSFieldValuesSetType,4,1,2.0_CMISSDP,Err)
  !element 5
  CALL CMISSFieldParameterSetUpdateElement(IndependentField,CMISSFieldIndependentType,&
	& CMISSFieldValuesSetType,5,1,1.5_CMISSDP,Err)
  !element 6
  CALL CMISSFieldParameterSetUpdateElement(IndependentField,CMISSFieldIndependentType,&
	& CMISSFieldValuesSetType,6,1,2.0_CMISSDP,Err)
  !element 7
  CALL CMISSFieldParameterSetUpdateElement(IndependentField,CMISSFieldIndependentType,&
	& CMISSFieldValuesSetType,7,1,3.0_CMISSDP,Err)
  !element 8
  CALL CMISSFieldParameterSetUpdateElement(IndependentField,CMISSFieldIndependentType,&
	& CMISSFieldValuesSetType,8,1,1.5_CMISSDP,Err)
  !element 9
  CALL CMISSFieldParameterSetUpdateElement(IndependentField,CMISSFieldIndependentType,&
	& CMISSFieldValuesSetType,9,1,2.0_CMISSDP,Err)
  !element 10
  CALL CMISSFieldParameterSetUpdateElement(IndependentField,CMISSFieldIndependentType,&
	& CMISSFieldValuesSetType,10,1,1.5_CMISSDP,Err)
  !element 11
  CALL CMISSFieldParameterSetUpdateElement(IndependentField,CMISSFieldIndependentType,&
	& CMISSFieldValuesSetType,11,1,1.0_CMISSDP,Err)

  
!Radius -cant finish ?
  CALL CMISSFieldParameterSetUpdateStart(IndependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(IndependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)

stop

  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"Test_Airway","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"Test_Airway","FORTRAN",Err)
!  CALL CMISSFieldIOElementsExport(IndependentField,"Radius","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)
!#################################################################################################################################################################
!notes not sure why I cant "finish" the element based field -it prints out ok to cmgui, should this have been an equations set ?


!##################################################################################################################################################################

  !Create the equations_set
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetFluidmechanicsClass, &
    & CMISSEquationsSetPoiseuilleEquationType,CMISSEquationsSetStaticPoiseuilleSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the equations set independent field variables, in this case radius
  CALL CMISSEquationsSetIndependentCreateStart(EquationsSet,IndependentFieldUserNumber,IndependentField,Err)
  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables pressure, flow
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field variables (length...need to set up differently)
!  CALL CMISSFieldTypeInitialise(MaterialsField,Err)
!  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
!  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)
!  PIPE_LENGTH=SQRT(DOT_PRODUCT(POSITION,POSITION))
!  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,PIPE_LENGTH,Err)

  !Create the equations set material field variables viscosity
  CALL CMISSFieldTypeInitialise(MaterialsFieldViscosity,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldViscosityUserNumber,MaterialsFieldViscosity,Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)
  VISCOSITY=1.86E-5_CMISSDP !viscosity of air
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldViscosity,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Viscosity,Err)


stop

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemFluidmechanicsClass,CMISSProblemPoiseuilleEquationType, &
    & CMISSProblemStaticPoiseuilleSubtype,Err)
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  !Set solver parameters
  !CALL CMISSSolverLinearTypeSet(Solver,CMISSSolverLinearIterativeSolveType,Err)
  CALL CMISSSolverLinearTypeSet(Solver,CMISSSolverLinearDirectSolveType,Err)
  !CALL CMISSSolverLinearIterativeRelativeToleranceSet(Solver,1.0E-8_CMISSDP,Err)
  !CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(Solver,1.0E-8_CMISSDP,Err)
  !CALL CMISSSolverLinearIterativeMaximumIterationsSet(Solver,10000,Err)
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Set up the boundary conditions
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the fixed boundary conditions at the first node and last nodes
  FirstNodeNumber=1
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSRegionNodesGet(Region,Nodes,Err)
  CALL CMISSNodesNumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,FirstNodeNumber,1, &
      & CMISSBoundaryConditionFixed,100.0_CMISSDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,LastNodeNumber,1, &
      & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  ENDIF
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)

  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"Poiseuille","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"Poiseuille","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

END PROGRAM AIRWAYPOISEUILLEEXAMPLE

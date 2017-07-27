! StructuredHexMesh_Class.f90
! 
! Copyright 2017 Joseph Schoonover <schoonover.numerics@gmail.com>
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE StructuredHexMesh_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE LinkedList_Class
USE KeyRing_Class
! src/interp/
USE Quadrature
USE Lagrange_Class
! src/geom/
USE Surface_Class
USE MappedGeometry_3D_Class
USE HexElement_Class
USE Face_Class
USE Node_Class



IMPLICIT NONE


 
   TYPE StructuredHexMesh 
      INTEGER                                 :: nXElem, nYElem, nZElem
      INTEGER                                 :: nFaces
      TYPE( HexElement ), ALLOCATABLE         :: elements(:,:,:) !
      TYPE( MappedGeometry_3D ), ALLOCATABLE  :: geom(:,:,:)
      TYPE( Node  ), ALLOCATABLE              :: nodes(:,:,:)  
      TYPE( Face ), ALLOCATABLE               :: faces(:)
      INTEGER                                 :: cornerMap(1:3,1:nHexNodes) 
      INTEGER                                 :: sideMap(1:nHexFaces) 
      INTEGER                                 :: faceMap(1:nQuadNodes,1:nHexFaces) 
      INTEGER                                 :: edgeFaceMap(1:2,1:nQuadEdges)

      CONTAINS
      
      PROCEDURE :: Initialize => Initialize_StructuredHexMesh
      PROCEDURE :: Trash      => Trash_StructuredHexMesh
       
      PROCEDURE :: ConstructFaces            => ConstructFaces_StructuredHexMesh
      PROCEDURE :: ConstructElementNeighbors => ConstructElementNeighbors_StructuredHexMesh
      PROCEDURE :: DetermineOrientation      => DetermineOrientation_StructuredHexMesh 
      PROCEDURE :: ScaleTheMesh              => ScaleTheMesh_StructuredHexMesh
      PROCEDURE :: LoadDefaultMesh           => LoadDefaultMesh_StructuredHexMesh
      
      PROCEDURE :: ConstructDoublyPeriodicFaces => ConstructDoublyPeriodicFaces_StructuredHexMesh
      PROCEDURE :: LoadDoublyPeriodicMesh       => LoadDoublyPeriodicMesh_StructuredHexMesh
       
      PROCEDURE :: ReadPeaceMeshFile    => ReadPeaceMeshFile_StructuredHexMesh
      PROCEDURE :: WritePeaceMeshFile   => WritePeaceMeshFile_StructuredHexMesh
      PROCEDURE :: ReadUCDMeshFile      => ReadUCDMeshFile_StructuredHexMesh
      PROCEDURE :: WriteTecplot         => WriteTecplot_Hexmesh
      PROCEDURE :: WriteMaterialTecplot => WriteMaterialTecplot_Hexmesh
       
    END TYPE StructuredHexMesh


 INTEGER, PRIVATE, PARAMETER    :: nDims = 3
 INTEGER, PRIVATE, PARAMETER    :: BoundaryFlagDefault = NO_NORMAL_FLOW

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup StructuredHexMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R Initialize 
! 
!> \fn Initialize_StructuredHexMesh 
!! Initializes the attributes and "convenience" arrays of the StructuredHexMesh data-structure.
!! 
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(StructuredHexMesh) :: this <BR>
!! <B>INTEGER</B>        :: nNodes, nElems, nFaces, N <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize( nNodes, nElems, nFaces, N ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myStructuredHexMesh <td> StructuredHexMesh <td> 
!!   <tr> <td> in <th> nNodes <td> INTEGER <td> The number of nodes in the mesh
!!   <tr> <td> in <th> nElems <td> INTEGER <td> The number of elements in the mesh
!!   <tr> <td> in <th> nFaces <td> INTEGER <td> The number of faces in the mesh 
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree for the geometry within each element
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_StructuredHexMesh( myStructuredHexMesh, nNodes, nElems, nFaces, N )

   IMPLICIT NONE
   CLASS(StructuredHexMesh), INTENT(out) :: myStructuredHexMesh
   INTEGER, INTENT(in)         :: nNodes, nElems, nFaces, N
   !LOCAL
   INTEGER :: i

      ! A hexahedron element (hex-element for short) has six faces. Each face has geometry that 
      ! requires the use of two computational coordinates. The third computational coordinate is 
      ! fixed. The sideMap gives the value of the remaining computational coordinate for each face
      ! of the hex-element.
      ! 
      ! The face ordering is  1 - South, 2 - East, 3 - North, 4 - West, 5 - Bottom, 6 - Top 
      myHexmesh % sideMap(1:nHexFaces) = (/ 0, N, N, 0, 0, N /)
      
      ! The eight corner nodes in an approximation that uses Gauss-Lobatto points, typically CG-type
      ! methods, have fixed computational coordinates. The corner node numbering  starts in the 
      ! southwest corner (in the computational grid) of the bottom face, proceeds counter clockwise
      ! around the face, and then repeats for the top face. This gives the local ID for the corner
      ! nodes as
      ! 
      ! Bottom, SouthWest = 1
      ! Bottom, SouthEast = 2
      ! Bottom, NorthEast = 3
      ! Bottom, NorthWest = 4
      ! Top, SouthWest = 5
      ! Top, SouthEast = 6
      ! Top, NorthEast = 7
      ! Top, NorthWest = 8
      !
      ! The computational coordinates for the corner nodes is given assuming a Gauss-Lobatto 
      ! computational mesh is used. Note that for a Gauss mesh, the corner nodes are not included.
      
      myHexmesh % cornerMap(1, 1:nHexNodes) = (/ 0, N, N,  0,  0, N, N,  0 /)
      myHexmesh % cornerMap(2, 1:nHexNodes) = (/ 0,  0, N, N,  0,  0, N, N /)
      myHexmesh % cornerMap(3, 1:nHexNodes) = (/ 0,  0,  0,  0, N, N, N, N /)
      
      ! Mesh construction usually begins with the specification of elements and the corner nodes, in
      ! addition to the element geometry. From the element-to-node connectivity, we need to construct
      ! the unique faces in the mesh and specify the abutting elements and their relative orientation.
      ! This procedure is aided by a "convenience array" that lists the local corner node IDs in the
      ! local counter-clockwise direction beginning in the local southwest corner of the face. When 
      ! two elements share a face, the global node IDs for each element can be found using this 
      ! convenience array (called "faceMap") and the relative orientation of the neighboring elements
      ! can be determined. The first index cycles over the nodes which make up the face in a 
      ! counterclockwise direction. The second index cycles over the faces in the element.
      
      myStructuredHexMesh % faceMap(1:nQuadNodes, south)  = (/ 1, 2, 6, 5 /)
      myStructuredHexMesh % faceMap(1:nQuadNodes, east)   = (/ 2, 3, 7, 6 /)
      myStructuredHexMesh % faceMap(1:nQuadNodes, north)  = (/ 4, 3, 7, 8 /)
      myStructuredHexMesh % faceMap(1:nQuadNodes, west)   = (/ 1, 4, 8, 5 /)
      myStructuredHexMesh % faceMap(1:nQuadNodes, bottom) = (/ 1, 2, 3, 4 /)
      myStructuredHexMesh % faceMap(1:nQuadNodes, top)    = (/ 5, 6, 7, 8 /)
      
      ! Each of the faces can be identified by their four corner nodes. The geometry of the faces
      ! is described using two computational coordinates between [-1,1]X[-1,1]. This 2-D computational
      ! grid has its own "southwest", "southeast", "northeast", and "northwest" identifications. The
      ! The corner nodes of the faces are labeled in the order mentioned in the previous sentence
      ! ( counterclockwise starting from southwest ). For quick referencing when producing tri-linear
      ! elements, a book-keeping array is useful for ensuring that we reference each edge of the 
      ! face in the order of increasing computational coordinate. This is identical to what is done
      ! in the StructuredHexMeshClass.f90 for the "edgeMap". Here, it is called the "edgeFaceMap".
      ! The first index references the starting(1) or ending(2) node. The second index references
      ! the edge of the face with the first being the southern edge and increasing the second index
      ! proceeds counter-clockwise.
      
      myHexmesh % edgeFaceMap(1, 1:nQuadEdges) = (/ 1, 2, 4, 1 /)
      myHexmesh % edgeFaceMap(2, 1:nQuadEdges) = (/ 2, 3, 3, 4 /)
      
      ! The number of nodes, the number of elements, and the number of faces are stored in this data
      ! structure for convenience. In another implementation (planned for the next version), the 
      ! number of elements, nodes, and faces is dynamic; that implementation 
      ! requires the use of dynamic storage, e.g. a linked-list like structure for elements, edges,
      ! and nodes.
      myStructuredHexMesh % nNodes = nNodes
      myStructuredHexMesh % nElems = nElems
      myStructuredHexMesh % nFaces = nFaces
      
      ALLOCATE( myHexmesh % elements(1:nElems) )
      ALLOCATE( myHexmesh % geom(1:nElems) )
      ALLOCATE( myHexmesh % nodes(1:nNodes) )
      ALLOCATE( myHexmesh % faces(1:nFaces) )
     
      DO i = 1, myHexmesh % nNodes
         CALL myHexmesh % nodes(i) % Initialize( )
      ENDDO 
      DO i = 1, myStructuredHexMesh % nElems
         CALL myStructuredHexMesh % elements(i) % Initialize( )
         CALL myStructuredHexMesh % geom(i) % Initialize( N )
      ENDDO
      DO i = 1, myStructuredHexMesh % nFaces
         CALL myStructuredHexMesh % faces(i) % Initialize( )
      ENDDO

 END SUBROUTINE Initialize_StructuredHexMesh
!
!> \addtogroup StructuredHexMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_StructuredHexMesh 
!! Frees memory associated with each of the attributes of the StructuredHexMesh data structure
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(StructuredHexMesh) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myStructuredHexMesh <td> StructuredHexMesh <td> On output, the memory associated with 
!!                         attributes of the StructuredHexMesh structure have been freed 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_StructuredHexMesh( myStructuredHexMesh )

   IMPLICIT NONE
   CLASS(StructuredHexMesh), INTENT(inout) :: myStructuredHexMesh
  ! LOCAL
   INTEGER :: i
      
      DO i = 1, myStructuredHexMesh % nElems
         CALL myStructuredHexMesh % geom(i) % Trash( )
      ENDDO

      DEALLOCATE( myStructuredHexMesh % nodes, &
                  myStructuredHexMesh % faces, &
                  myStructuredHexMesh % geom, &
                  myStructuredHexMesh % elements )

 END SUBROUTINE Trash_StructuredHexMesh
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup StructuredHexMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R ConstructFaces
! 
!> \fn ConstructFaces_StructuredHexMesh  
!! Uses the element-to-node connectivity to construct all of the unique faces in a mesh
!! 
!! Similar to the edge construction in a QuadMesh, we loop over the elements and construct a face
!! using the element-to-node connectivity and a convenience array that relates the local nodes
!! to each face. A face is identified by its four corner nodes. The four corner node IDs for each
!! elements face are gathered and compared against a list of already identified faces. If the
!! face is not in the list, then a new face is added to the list and the element that added
!! the face to the list is designated the "primary" element. If the face is on the list,
!! the orientation of the secondary element relative to the primary element is determined.
!! 
!! This routine depends on <BR>
!!     Module StructuredHexMesh_Class : \ref DetermineOrientation_StructuredHexMesh
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DataType) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ConstructFaces( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myStructuredHexMesh <td> StructuredHexMesh <td> 
!!                         On <B>input</B>, the element and node information has been filled in, <BR>
!!                         On <B>output</B>, the unique faces have been identified and the face 
!!                         information has been filled in. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ConstructFaces_StructuredHexMesh( myStructuredHexMesh )

   IMPLICIT NONE
   CLASS( StructuredHexMesh ), INTENT(inout) :: myStructuredHexMesh
   ! LOCAL
   TYPE( KeyRing ) :: KeyCabinet(1:myStructuredHexMesh % nNodes)
   INTEGER :: nEls, nNodes, iEl, nFaces, k, j  
   INTEGER :: locNodeIDs(1:nQuadNodes), globNodeIDs(1:nQuadNodes)
   INTEGER :: keyRingID, globFaceID
   INTEGER :: N
   LOGICAL :: keyExists

      nNodes = myStructuredHexMesh % nNodes
      nEls   = myStructuredHexMesh % nElems
      N      = myStructuredHexMesh % geom(1) % N
      nFaces = 0
      
      DO k = 1, nNodes
         CALL keyCabinet(k) % Build( )
      ENDDO

      DO iEl = 1, nEls ! Loop over the elements in the mesh

         DO k = 1, nHexFaces ! Loop over the faces of each element

            ! In this first step, we want to identify the unique global node ID's for this face
            ! To do this, we start by gathering the local (to the element) node ID's.
            DO j = 1, nQuadNodes   
               locNodeIDs(j) = myStructuredHexMesh % faceMap(j,k) ! starting local node for this Face
            ENDDO
            
            ! Now, we extract, for this element, the global node ID's for each local node ID
            DO j = 1, nQuadNodes 
               globNodeIDs(j) = myStructuredHexMesh % elements(iEl) % nodeIDs( locNodeIDs(j) )
            ENDDO
            ! Our key cabinet has many key-rings with each key-ring containing a set of notched keys.
            ! Each notched key corresponds to a unique face in the mesh. To enable fast searching
            ! for a unique face, we address our key-rings according to the minimum global node ID
            ! that resides on the current face. If another face shares this minimum global node ID,
            ! then it is possible that the face has already been generated. If not, our search is 
            ! limited only to the key-ring with the same key-ring ID. Here, we grab the key-ring ID.
            keyRingID = MINVAL( globNodeIDs )
            
            ! Now, we check to see if a notched key already exists with the same set of global ID's
            ! This is done by making a call to "FindDataForNotches". This routine searches through
            ! the current key ring for the key which has the same global ID's (though not 
            ! necessarily in the same order). If the face has been built already, the face ID is 
            ! returned in globFaceID and the logical, "keyExists", is set to TRUE. Otherwise, 
            ! globFaceID is set to zero and "keyExsists" is set to FALSE.
            CALL KeyCabinet(keyRingID) % FindDataForNotches( globNodeIDS, nQuadNodes, &
                                                             globFaceID, keyExists )
            ! Here is where the conditional processing begins. 
            ! 
            ! If this face has already been found then we do nothing
            !
            ! If this is a new face, we increment the number of faces and add to the key-ring.

            ! Add element to corner-node connectivity list
            IF( .NOT.(keyExists) )then ! this is a new face

               nFaces = nFaces + 1
               CALL KeyCabinet(keyRingID) % AddToList( nFaces, globNodeIDs, nQuadNodes )
               
            ENDIF
            

         ENDDO ! k, Loop over the faces of each element
        
      ENDDO! iEl, Loop over the elements in the mesh
 

      DO k = 1, nNodes
         CALL keyCabinet(k) % Trash( ) ! Trash the Facetable
         CALL keyCabinet(k) % Build( ) ! and rebuild a blank cabinet 
      ENDDO
      
      ! Re-allocate space for the mesh Faces

      DEALLOCATE( myStructuredHexMesh % Faces )

      ALLOCATE( myStructuredHexMesh % Faces( 1:nFaces ) )
      DO k = 1, nFaces
         CALL myStructuredHexMesh % faces(k) % Initialize( )
      ENDDO
      nFaces = 0

      DO iEl = 1, nEls ! Loop over the elements in the mesh

         DO k = 1, nHexFaces ! Loop over the faces of each element

            ! In this first step, we want to identify the unique global node ID's for this face
            ! To do this, we start by gathering the local (to the element) node ID's.
            DO j = 1, nQuadNodes   
               locNodeIDs(j) = myStructuredHexMesh % faceMap(j,k) ! starting local node for this Face
            ENDDO
            
            ! Now, we extract, for this element, the global node ID's for each local node ID
            DO j = 1, nQuadNodes 
               globNodeIDs(j) = myStructuredHexMesh % elements(iEl) % nodeIDs(locNodeIDs(j))
            ENDDO
            
            ! Our key cabinet has many key-rings with each key-ring containing a set of notched keys.
            ! Each notched key corresponds to a unique face in the mesh. To enable fast searching
            ! for a unique face, we address our key-rings according to the minimum global node ID
            ! that resides on the current face. If another face shares this minimum global node ID,
            ! then it is possible that the face has already been generated. If not, our search is 
            ! limited only to the key-ring with the same key-ring ID. Here, we grab the key-ring ID.
            keyRingID = MINVAL( globNodeIDs )
            
            ! Now, we check to see if a notched key already exists with the same set of global ID's
            ! This is done by making a call to "FindDataForNotches". This routine searches through
            ! the current key ring for the key which has the same global ID's (though not 
            ! necessarily in the same order). If the face has been built already, the face ID is 
            ! returned in globFaceID and the logical, "keyExists", is set to TRUE. Otherwise, 
            ! globFaceID is set to zero and "keyExsists" is set to FALSE.
            CALL KeyCabinet(keyRingID) % FindDataForNotches( globNodeIDS, nQuadNodes, &
                                                             globFaceID, keyExists )

            ! Here is where the conditional processing begins. 
            ! 
            ! If this face has already been found then we need to determine the face orientation
            ! of the secondary element relative to the primary element.
            !
            ! If this is a new face, we set the primary element information, and default the
            ! secondary element information.

            ! Add element to corner-node connectivity list
            IF( keyExists )then ! this face has already been found

               ! Since this face exists, we need to compare the relative face orientation of the 
               ! secondary element to the primary element. .

               CALL myStructuredHexMesh % DetermineOrientation( globFaceID, globNodeIDs )

               myStructuredHexMesh % faces( globFaceID ) % elementIDs(2)   = iEl
               myStructuredHexMesh % faces( globFaceID ) % elementSides(2) = k
               
            ELSE ! This is a new face

               ! First, we store the key-ring information
               nFaces = nFaces + 1
               CALL KeyCabinet(keyRingID) % AddToList( nFaces, globNodeIDs, nQuadNodes )

               ! Now, we set the primary element information
               myStructuredHexMesh % faces( nFaces ) % faceID          = nFaces
               myStructuredHexMesh % faces( nFaces ) % nodeIDs         = globNodeIDs
               myStructuredHexMesh % faces( nFaces ) % elementIDs(1)   = iEl
               myStructuredHexMesh % faces( nFaces ) % elementSides(1) = k
               ! Now we default the secondary element information and the swap flag
               myStructuredHexMesh % faces( nFaces ) % elementIDs(2)   = BoundaryFlagDefault
               myStructuredHexMesh % faces( nFaces ) % elementSides(2) = k
               myStructuredHexMesh % faces( nFaces ) % iStart          = 0
               myStructuredHexMesh % faces( nFaces ) % iInc            = 1
               myStructuredHexMesh % faces( nFaces ) % jStart          = 0
               myStructuredHexMesh % faces( nFaces ) % jInc            = 1
               myStructuredHexMesh % faces( nFaces ) % swapDimensions  = 0

            ENDIF

         ENDDO ! k, Loop over the faces of each element
        
      ENDDO! iEl, Loop over the elements in the mesh
      DO k = 1, nNodes
         CALL keyCabinet(k) % Trash( ) ! Trash the Facetable
      ENDDO
      myStructuredHexMesh % nFaces = nFaces
      

 END SUBROUTINE ConstructFaces_StructuredHexMesh
! 
SUBROUTINE ConstructDoublyPeriodicFaces_StructuredHexMesh( myStructuredHexMesh, nXElem, nYElem, nZElem )

! Assumes structured mesh
   IMPLICIT NONE
   CLASS( StructuredHexMesh ), INTENT(inout) :: myStructuredHexMesh
   INTEGER, INTENT(in) :: nXElem, nYElem, nZElem
   ! LOCAL
   INTEGER ::  e1, e2, s1, s2, nFaces, i, j, k, l, iFace

      nFaces = (nZElem+1)*nXElem*nYElem + (nXElem+1)*nYElem*nZElem + (nYElem+1)*nXElem*nZElem
      
      ! Re-allocate space for the mesh Faces
      DEALLOCATE( myStructuredHexMesh % Faces )
      ALLOCATE( myStructuredHexMesh % Faces( 1:nFaces ) )
      DO k = 1, nFaces
         CALL myStructuredHexMesh % faces(k) % Initialize( )
      ENDDO
      myStructuredHexMesh % nFaces = nFaces
      iFace = 0
      
      DO k = 1, nZElem
         DO j = 1, nYElem
            DO i = 1, nXElem
               
                  e1 = i + nXElem*( j-1 + nYElem*(k-1) ) ! Primary element ID
                  ! Element e1's southern boundary
                  s1 = SOUTH 
                  IF( j==1 )THEN
                     iFace = iFace + 1
                     e2 = i + nXElem*( nYElem-1 + (nYElem)*(k-1) )! Enforce periodicity with the "northern" most element
                     s2 = NORTH
                     myStructuredHexMesh % faces(iFace) % faceID = iFace
                     DO l = 1, 4
                        myStructuredHexMesh % faces(iFace) % nodeIDs(l) = myStructuredHexMesh % elements(e1) % nodeIDs( myStructuredHexMesh % faceMap(l,SOUTH) )
                     ENDDO
                     myStructuredHexMesh % faces(iFace) % elementIDs(1) = e1
                     myStructuredHexMesh % faces(iFace) % elementIDs(2) = e2
                     myStructuredHexMesh % faces(iFace) % elementSides(1) = s1
                     myStructuredHexMesh % faces(iFace) % elementSides(2) = s2
                     myStructuredHexMesh % faces(iFace) % iStart          = 0
                     myStructuredHexMesh % faces(iFace) % iInc            = 1
                     myStructuredHexMesh % faces(iFace) % jStart          = 0
                     myStructuredHexMesh % faces(iFace) % jInc            = 1
                     myStructuredHexMesh % faces(iFace) % swapDimensions  = 0
                  ELSE
                     iFace = iFace + 1
                     e2 = i + nXElem*( j-2 + (nYElem)*(k-1) )
                     s2 = NORTH
                     myStructuredHexMesh % faces(iFace) % faceID = iFace
                     DO l = 1, 4
                        myStructuredHexMesh % faces(iFace) % nodeIDs(l) = myStructuredHexMesh % elements(e1) % nodeIDs( myStructuredHexMesh % faceMap(l,SOUTH) )
                     ENDDO
                     myStructuredHexMesh % faces(iFace) % elementIDs(1) = e1
                     myStructuredHexMesh % faces(iFace) % elementIDs(2) = e2
                     myStructuredHexMesh % faces(iFace) % elementSides(1) = s1
                     myStructuredHexMesh % faces(iFace) % elementSides(2) = s2
                     myStructuredHexMesh % faces(iFace) % iStart          = 0
                     myStructuredHexMesh % faces(iFace) % iInc            = 1
                     myStructuredHexMesh % faces(iFace) % jStart          = 0
                     myStructuredHexMesh % faces(iFace) % jInc            = 1
                     myStructuredHexMesh % faces(iFace) % swapDimensions  = 0
                  ENDIF
                  
                  ! Element e1's western boundary
                  s1 = WEST
                  IF( i==1 )THEN
                     iFace = iFace + 1
                     e2 = nXElem + nXElem*( j-1 + (nYElem)*(k-1) ) ! Enforce periodicity with the "eastern" most element
                     s2 = EAST
                     myStructuredHexMesh % faces(iFace) % faceID = iFace
                     DO l = 1, 4
                        myStructuredHexMesh % faces(iFace) % nodeIDs(l) = myStructuredHexMesh % elements(e1) % nodeIDs( myStructuredHexMesh % faceMap(l,WEST) )
                     ENDDO
                     myStructuredHexMesh % faces(iFace) % elementIDs(1) = e1
                     myStructuredHexMesh % faces(iFace) % elementIDs(2) = e2
                     myStructuredHexMesh % faces(iFace) % elementSides(1) = s1
                     myStructuredHexMesh % faces(iFace) % elementSides(2) = s2
                     myStructuredHexMesh % faces(iFace) % iStart          = 0
                     myStructuredHexMesh % faces(iFace) % iInc            = 1
                     myStructuredHexMesh % faces(iFace) % jStart          = 0
                     myStructuredHexMesh % faces(iFace) % jInc            = 1
                     myStructuredHexMesh % faces(iFace) % swapDimensions  = 0
                  ELSE
                     iFace = iFace + 1
                     e2 = i-1 + nXElem*( j-1 + (nYElem)*(k-1) )
                     s2 = EAST
                     myStructuredHexMesh % faces(iFace) % faceID = iFace
                     DO l = 1, 4
                        myStructuredHexMesh % faces(iFace) % nodeIDs(l) = myStructuredHexMesh % elements(e1) % nodeIDs( myStructuredHexMesh % faceMap(l,WEST) )
                     ENDDO
                     myStructuredHexMesh % faces(iFace) % elementIDs(1) = e1
                     myStructuredHexMesh % faces(iFace) % elementIDs(2) = e2
                     myStructuredHexMesh % faces(iFace) % elementSides(1) = s1
                     myStructuredHexMesh % faces(iFace) % elementSides(2) = s2
                     myStructuredHexMesh % faces(iFace) % iStart          = 0
                     myStructuredHexMesh % faces(iFace) % iInc            = 1
                     myStructuredHexMesh % faces(iFace) % jStart          = 0
                     myStructuredHexMesh % faces(iFace) % jInc            = 1
                     myStructuredHexMesh % faces(iFace) % swapDimensions  = 0
                  ENDIF
                  
                  ! Element e1's bottom boundary
                  s1 = BOTTOM
                  IF( k==1 )THEN
                     iFace = iFace + 1
                     e2 = NO_NORMAL_FLOW
                     s2 = TOP
                     myStructuredHexMesh % faces(iFace) % faceID = iFace
                     DO l = 1, 4
                        myStructuredHexMesh % faces(iFace) % nodeIDs(l) = myStructuredHexMesh % elements(e1) % nodeIDs( myStructuredHexMesh % faceMap(l,BOTTOM) )
                     ENDDO
                     myStructuredHexMesh % faces(iFace) % elementIDs(1) = e1
                     myStructuredHexMesh % faces(iFace) % elementIDs(2) = e2
                     myStructuredHexMesh % faces(iFace) % elementSides(1) = s1
                     myStructuredHexMesh % faces(iFace) % elementSides(2) = s2
                     myStructuredHexMesh % faces(iFace) % iStart          = 0
                     myStructuredHexMesh % faces(iFace) % iInc            = 1
                     myStructuredHexMesh % faces(iFace) % jStart          = 0
                     myStructuredHexMesh % faces(iFace) % jInc            = 1
                     myStructuredHexMesh % faces(iFace) % swapDimensions  = 0
                  ELSE
                     iFace = iFace + 1
                     e2 = i + nXElem*( j-1 + (nYElem)*(k-2) )
                     s2 = TOP
                     myStructuredHexMesh % faces(iFace) % faceID = iFace
                     DO l = 1, 4
                        myStructuredHexMesh % faces(iFace) % nodeIDs(l) = myStructuredHexMesh % elements(e1) % nodeIDs( myStructuredHexMesh % faceMap(l,BOTTOM) )
                     ENDDO
                     myStructuredHexMesh % faces(iFace) % elementIDs(1) = e1
                     myStructuredHexMesh % faces(iFace) % elementIDs(2) = e2
                     myStructuredHexMesh % faces(iFace) % elementSides(1) = s1
                     myStructuredHexMesh % faces(iFace) % elementSides(2) = s2
                     myStructuredHexMesh % faces(iFace) % iStart          = 0
                     myStructuredHexMesh % faces(iFace) % iInc            = 1
                     myStructuredHexMesh % faces(iFace) % jStart          = 0
                     myStructuredHexMesh % faces(iFace) % jInc            = 1
                     myStructuredHexMesh % faces(iFace) % swapDimensions  = 0
                  ENDIF
                  
            ENDDO ! i
            
            e1 = nXElem + nXElem*( j-1 + (nYElem)*(k-1) )
            s1 = EAST
            iFace = iFace + 1
            e2 = 1 + nXElem*( j-1 + (nYElem)*(k-1) ) ! Enforce periodicity with the "western" most element
            s2 = WEST
            myStructuredHexMesh % faces(iFace) % faceID = iFace
            DO l = 1, 4
               myStructuredHexMesh % faces(iFace) % nodeIDs(l) = myStructuredHexMesh % elements(e1) % nodeIDs( myStructuredHexMesh % faceMap(l,EAST) )
            ENDDO
            myStructuredHexMesh % faces(iFace) % elementIDs(1) = e1
            myStructuredHexMesh % faces(iFace) % elementIDs(2) = e2
            myStructuredHexMesh % faces(iFace) % elementSides(1) = s1
            myStructuredHexMesh % faces(iFace) % elementSides(2) = s2
            myStructuredHexMesh % faces(iFace) % iStart          = 0
            myStructuredHexMesh % faces(iFace) % iInc            = 1
            myStructuredHexMesh % faces(iFace) % jStart          = 0
            myStructuredHexMesh % faces(iFace) % jInc            = 1
            myStructuredHexMesh % faces(iFace) % swapDimensions  = 0
            
         ENDDO ! j
         
         DO i = 1, nXElem
            iFace = iFace + 1
            e1 = i + nXElem*( nYElem-1 + (nYElem)*(k-1) )
            s1 = NORTH
            e2 = i + nXElem*( nYElem*(k-1) )
            s2 = SOUTH
            myStructuredHexMesh % faces(iFace) % faceID = iFace
            DO l = 1, 4
               myStructuredHexMesh % faces(iFace) % nodeIDs(l) = myStructuredHexMesh % elements(e1) % nodeIDs( myStructuredHexMesh % faceMap(l,NORTH) )
            ENDDO
            myStructuredHexMesh % faces(iFace) % elementIDs(1) = e1
            myStructuredHexMesh % faces(iFace) % elementIDs(2) = e2
            myStructuredHexMesh % faces(iFace) % elementSides(1) = s1
            myStructuredHexMesh % faces(iFace) % elementSides(2) = s2
            myStructuredHexMesh % faces(iFace) % iStart          = 0
            myStructuredHexMesh % faces(iFace) % iInc            = 1
            myStructuredHexMesh % faces(iFace) % jStart          = 0
            myStructuredHexMesh % faces(iFace) % jInc            = 1
            myStructuredHexMesh % faces(iFace) % swapDimensions  = 0
         ENDDO

      ENDDO ! k
         
      DO j = 1, nYElem
         DO i = 1, nXElem
            e1 = i + nXElem*( j-1 + (nYElem)*(nZElem-1) ) ! Primary element ID
            iFace = iFace + 1
            e2 = PRESCRIBED
            s1 = TOP
            s2 = s1
            myStructuredHexMesh % faces(iFace) % faceID = iFace
            DO l = 1, 4
               myStructuredHexMesh % faces(iFace) % nodeIDs(l) = myStructuredHexMesh % elements(e1) % nodeIDs( myStructuredHexMesh % faceMap(l,TOP) )
            ENDDO
            myStructuredHexMesh % faces(iFace) % elementIDs(1) = e1
            myStructuredHexMesh % faces(iFace) % elementIDs(2) = e2
            myStructuredHexMesh % faces(iFace) % elementSides(1) = s1
            myStructuredHexMesh % faces(iFace) % elementSides(2) = s2
            myStructuredHexMesh % faces(iFace) % iStart          = 0
            myStructuredHexMesh % faces(iFace) % iInc            = 1
            myStructuredHexMesh % faces(iFace) % jStart          = 0
            myStructuredHexMesh % faces(iFace) % jInc            = 1
            myStructuredHexMesh % faces(iFace) % swapDimensions  = 0
         ENDDO
      ENDDO
   
      PRINT*, 'iFace : ',iFace

 END SUBROUTINE ConstructDoublyPeriodicFaces_StructuredHexMesh
!
!> \addtogroup StructuredHexMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R DetermineOrientation
! 
!> \fn DetermineOrientation_StructuredHexMesh
!! Given a face-ID and an array of global node ID's, this routine determines the relative orientation
!! of the face to a face represented by the incoming global node ID's.
!! 
!!  This support routine takes as input the mesh, a face ID number, and a list of node IDs.
!!  This routine assumes that the primary element information for the given face ID has been
!!  filled upon issuing a call to this routine. If the primary element shares this face with
!!  another element, this routine is called and the secondary element node ID's are passed in.
!!  The order of the secondary node ID's relative to the order of the primary node ID's determines
!!  the orientation of the secondary element relative to the primary element. We need to know
!!  if the roles of the computational coordinates are flipped, and (simultaneously) if the 
!!  incrementing of the computational coordinates are reversed (or not).
!!  This routine determines the orientation by determining how many times the ordered secondary
!!  node ID's need to be shifted in order to match the ordered primary node ID's.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(StructuredHexMesh) :: this <BR>
!! <B>INTEGER</B>       :: faceID <BR>
!! <B>INTEGER</B>       :: secondaryNodes(1:4) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % DetermineOrientation( faceID, secondaryNodes ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myStructuredHexMesh <td> StructuredHexMesh <td> 
!!                         On <B>input</B>, StructuredHexMesh with myStructuredHexMesh % faces(faceID) primary element 
!!                         information filled in, <BR>
!!                         On <B>output</B>, secondary element information is filled in for this face
!!   <tr> <td> in <th> faceID <td> INTEGER <td> ID for the face in question
!!   <tr> <td> in <th> secondaryNodes(1:4) <td> INTEGER <td> Global node ID's associated with the
!!                                                           secondary element.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE DetermineOrientation_StructuredHexMesh( myStructuredHexMesh, faceID, secondaryNodes ) 

   IMPLICIT NONE
   CLASS( StructuredHexMesh ), INTENT(inout) :: myStructuredHexMesh
   INTEGER, INTENT(in)             :: faceID
   INTEGER, INTENT(in)             :: secondaryNodes(1:nQuadNodes)
   ! Local
   INTEGER :: primaryNodes(1:nQuadNodes)
   INTEGER :: nShifts, i, N
   LOGICAL :: theyMatch 
   
      primaryNodes = myStructuredHexMesh % faces( faceID ) % nodeIDs
      N = myStructuredHexMesh % geom(1) % N
      nShifts = 0
      theyMatch = .FALSE.
      
      DO i = 1, nQuadNodes

         ! First, we compare the primary and secondary nodes. This routine returns a zero 
         ! if the arrays match, and a one if they do not match.
         theyMatch = CompareArray( primaryNodes, secondaryNodes, nQuadNodes )
         
         IF( theyMatch )THEN
            EXIT
         ELSE
            nShifts = nShifts + 1
            CALL ForwardShift( primaryNodes, nQuadNodes )
         ENDIF
         
      ENDDO
      
      IF( theyMatch )THEN
      
         SELECT CASE ( nShifts )
         
            CASE (0)
               myStructuredHexMesh % faces( faceID ) % iStart          = 0
               myStructuredHexMesh % faces( faceID ) % iInc            = 1
               myStructuredHexMesh % faces( faceID ) % jStart          = 0
               myStructuredHexMesh % faces( faceID ) % jInc            = 1
               myStructuredHexMesh % faces( faceID ) % swapDimensions  = 0
            CASE (1)
               myStructuredHexMesh % faces( faceID ) % iStart          = 0
               myStructuredHexMesh % faces( faceID ) % iInc            = 1
               myStructuredHexMesh % faces( faceID ) % jStart          = N
               myStructuredHexMesh % faces( faceID ) % jInc            = -1
               myStructuredHexMesh % faces( faceID ) % swapDimensions  = 1
            CASE (2)
               myStructuredHexMesh % faces( faceID ) % iStart          = N
               myStructuredHexMesh % faces( faceID ) % iInc            = -1
               myStructuredHexMesh % faces( faceID ) % jStart          = N
               myStructuredHexMesh % faces( faceID ) % jInc            = -1
               myStructuredHexMesh % faces( faceID ) % swapDimensions  = 0
            CASE (3)
               myStructuredHexMesh % faces( faceID ) % iStart          = N
               myStructuredHexMesh % faces( faceID ) % iInc            = -1
               myStructuredHexMesh % faces( faceID ) % jStart          = 0
               myStructuredHexMesh % faces( faceID ) % jInc            = 1
               myStructuredHexMesh % faces( faceID ) % swapDimensions  = 1
            CASE DEFAULT
               PRINT*, 'Module StructuredHexMeshClass : S/R DetermineOrientation : Did not catch a match. Revise this subroutine. '
               PRINT*, 'Also, this is a reminder to add an exception-handler. Stopping. '
               STOP
               
         END SELECT
      ELSE
      
         PRINT*, 'Module StructuredHexMeshClass : S/R DetermineOrientation : Did not catch a match. Revise this subroutine. '
         PRINT*, 'Also, this is a reminder to add an exception-handler. Stopping. '
         STOP
         
      ENDIF
   

 END SUBROUTINE DetermineOrientation_StructuredHexMesh
!
!> \addtogroup StructuredHexMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R ConstructElementNeighbors 
! 
!> \fn sConstructElementNeighbors_StructuredHexMesh 
!! Uses the edge-to-element connectivity to construct the element neighbors attribute.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(myStructuredHexMesh) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ConstructElementNeighbors( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myStructuredHexMesh <td> StructuredHexMesh <td>
!!                         On <B>input</B>, the edge-to-element information has been constructed, <BR>
!!                         On <B>output</B>, the element neighbors information has been filled in. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ConstructElementNeighbors_StructuredHexMesh( myStructuredHexMesh )

   IMPLICIT NONE
   CLASS( StructuredHexMesh ), INTENT(inout) :: myStructuredHexMesh
   ! LOCAL
   INTEGER :: e(1:2), s(1:2), nFaces, iFace
  
      nFaces = myStructuredHexMesh % nFaces
      DO iFace = 1, nFaces
         e = myStructuredHexMesh % faces( iFace ) % elementIDs
         s = myStructuredHexMesh % faces( iFace ) % elementSides
         IF( e(2) > 0 )THEN
            myStructuredHexMesh % elements(e(1)) % neighbors(s(1))      = e(2)
            myStructuredHexMesh % elements(e(2)) % neighbors(ABS(s(2))) = e(1)
         ENDIF
      ENDDO

     
 END SUBROUTINE ConstructElementNeighbors_StructuredHexMesh
!
!> \addtogroup StructuredHexMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R ScaleTheMesh_StructuredHexMesh 
! 
!> \fn ScaleTheMesh  
!! Scales the element geometry and corner node positions by the provided x-scale, y-scale, and z-scale. 
!! 
!! This routine depend on <BR>
!!   Module \ref MappedGeometry_2D_Class, S/R ScaleGeometry_MappedGeometry_2D <BR>
!!   Module \ref Node_Class, S/R ScalePosition_Node 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(StructuredHexMesh)    :: this <BR>
!! <B>TYPE</B>(Lagrange) :: interp <BR>
!! <B>REAL</B>(prec)        :: xScale, yScale, zScale <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ScaleGeometry( interp, xScale, yScale, zScale ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myStructuredHexMesh <td> StructuredHexMesh <td> 
!!   <tr> <td> in <th> interp <td> Lagrange <td> 
!!   <tr> <td> in <th> xScale <td> REAL(prec) <td>
!!   <tr> <td> in <th> yScale <td> REAL(prec) <td>  
!!   <tr> <td> in <th> zScale <td> REAL(prec) <td>  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ScaleTheMesh_StructuredHexMesh( myStructuredHexMesh, interp, xScale, yScale, zScale  )

   IMPLICIT NONE
   CLASS( StructuredHexMesh ), INTENT(inout) :: myStructuredHexMesh
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda ), INTENT(in) :: interp
#else
   TYPE( Lagrange ), INTENT(in)    :: interp
#endif
   REAL(prec), INTENT(in)          :: xScale, yScale, zScale
   ! LOCAL
   INTEGER :: nElems, nNodes, i
   
      nElems = myStructuredHexMesh % nElems
      nNodes = myStructuredHexMesh % nNodes

      DO i = 1, nElems
         CALL myStructuredHexMesh % geom(i) % ScaleGeometry( interp, xScale, yScale, zScale )
      ENDDO
      
      DO i = 1, nNodes
         CALL myStructuredHexMesh % nodes(i) % ScalePosition( xScale, yScale, zScale )
      ENDDO

 END SUBROUTINE ScaleTheMesh_StructuredHexMesh
!
!> \addtogroup StructuredHexMesh_Class
!! @{ 
! ================================================================================================ !
! S/R LoadDefaultMesh
! 
!> \fn LoadDefaultMesh_StructuredHexMesh  
!! Constructs a "structured" spectral element mesh with nXelem-by-nYelem elements. 
!! 
!! This routine builds a mesh with nXelem elements in the x-direction,  nYelem elements in the
!! y-direction, and nZelem elements in the z-direction. The mesh nodes are between [0,1]x[0,1]x[0,1].
!!  After construction, the user can call "ScaleTheMesh" to change the physical extents of the domain.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(StructuredHexMesh)    :: this <BR>
!! <B>TYPE</B>(Lagrange) :: interp <BR>
!! <B>INTEGER</B>           :: nXElem, nYElem, nZElem <BR>
!!         .... <BR>
!!     <B>CALL</B> this % RoutineName( interp, nXElem, nYElem, nZElem ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myStructuredHexMesh <td> StructuredHexMesh <td> On output, contains the "structured" mesh
!!                                                       in the unstructured spectral element
!!                                                       mesh format.
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D interpolant that stores the computational
!!                                                  quadrature mesh for each element.
!!   <tr> <td> in <th> nXElem <td> INTEGER <td> The number of desired elements in the x-direction
!!   <tr> <td> in <th> nYElem <td> INTEGER <td> The number of desired elements in the y-direction
!!   <tr> <td> in <th> nZElem <td> INTEGER <td> The number of desired elements in the z-direction
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE LoadDefaultMesh_StructuredHexMesh( myStructuredHexMesh, interp, nXelem, nYelem, nZelem  )

   IMPLICIT NONE
   CLASS( StructuredHexMesh ), INTENT(inout)  :: myStructuredHexMesh
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda ), INTENT(in)    :: interp
#else
   TYPE( Lagrange ), INTENT(in)     :: interp
#endif
   INTEGER, INTENT(in)              :: nXelem, nYelem, nZelem
   ! LOCAL
   TYPE( Surface ) :: boundSurfs(1:nHexFaces)

   REAL(prec) :: x, y, z, dxElem, dyElem, dzElem
   REAL(prec) :: x1(1:nDims), x2(1:nDims), x3(1:nDims), x4(1:nDims)
   REAL(prec) :: c1(1:nDims), c2(1:nDims)
   REAL(prec), ALLOCATABLE :: xc(:,:,:), s(:)

   INTEGER :: nNodes, nElems, nFaces, gPolyDeg
   INTEGER :: nodes(1:nHexNodes)
   INTEGER :: n1, n2, n3, n4
   INTEGER :: iNode, iEl, iSide, iX, iY, iZ, i, j

      
      dxElem = ONE/nXElem
      dyElem = ONE/nYElem
      dzElem = ONE/nZElem
      
      ! ** "Hard-wired" values for a structured mesh with no holes ** !
      nNodes   = (nXElem+1)*(nYElem+1)*(nZElem+1)
      nElems   = (nXElem)*(nYElem)*(nZElem)
      nFaces   = (nXElem)*(nYElem)*(nZElem+1) + (nXElem)*(nZElem)*(nYElem+1) + (nYElem)*(nZElem)*(nXElem+1)
      gPolyDeg = 1
      ! ************************************************************************* !

      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'gPolyDeg  : ', gPolyDeg

      ! Generate the chebyshev points of order gPolyDeg
      ! These are the points used to define the parametric
      ! curves for the element boundaries

      ALLOCATE( s(0:gPolyDeg), xc(0:gPolyDeg,0:gPolyDeg,1:3) )
      s(0) = -ONE
      s(1) = ONE
      
      ! ---- Initialize the mesh (empty) ---- !
      CALL myStructuredHexMesh % Initialize( nNodes, nElems, nFaces, interp % N ) 
      
      
      ! ---- Read in the corner nodes ---- !
      DO iZ = 1, nZElem + 1
         z = dZElem*(REAL(iZ-1,prec))
         DO iY = 1, nYElem+1
            y = dYElem*(REAL(iY-1,prec))
            DO iX = 1, nXElem+1
               iNode = iX + (iY-1)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)
               x = dXElem*(REAL(iX-1,prec))
               myStructuredHexMesh % nodes(iNode) % x = x
               myStructuredHexMesh % nodes(iNode) % y = y
               myStructuredHexMesh % nodes(iNode) % z = z
               myStructuredHexMesh % nodes(iNode) % nodeID = iNode
            ENDDO
         ENDDO
      ENDDO
  
      ! Do the element information
      xc = ZERO
      ! Do the initial build for the parametric surfaces
      DO iSide = 1, nHexFaces
         CALL boundSurfs(iSide) % Build( xc, s, gPolyDeg, nDims ) 
      ENDDO
   
      DO iZ = 1, nZElem
         DO iY = 1, nYElem
            DO iX = 1, nXElem

               iEl = iX + (iY-1)*(nXElem) + (iZ-1)*(nXElem)*(nYElem)
               ! Calculate the global node IDs for this element.
               nodes(1) = iX + (iY-1)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)    ! Southwest
               nodes(2) = iX + 1 + (iY-1)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)! SouthEast
               nodes(3) = iX + 1 + (iY)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)  ! NorthEast
               nodes(4) = iX + (iY)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)      ! NorthWest
               nodes(5) = iX + (iY-1)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)      ! Southwest
               nodes(6) = iX + 1 + (iY-1)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)  ! SouthEast
               nodes(7) = iX + 1 + (iY)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)    ! NorthEast
               nodes(8) = iX + (iY)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)        ! NorthWest

               myStructuredHexMesh % elements(iEl) % nodeIDs = nodes
         
               DO iSide = 1, nHexFaces ! Loop over the sides of the quads

                  ! To build the current face, we construct a plane that passes through
                  ! the four corner nodes. Here, we grab the global node ID's for the four
                  ! corner nodes.
                  n1 = nodes( myStructuredHexMesh % faceMap(1,iSide) )
                  n2 = nodes( myStructuredHexMesh % faceMap(2,iSide) )
                  n3 = nodes( myStructuredHexMesh % faceMap(3,iSide) ) 
                  n4 = nodes( myStructuredHexMesh % faceMap(4,iSide) ) 

                  x1(1) = myStructuredHexMesh % nodes(n1) % x
                  x1(2) = myStructuredHexMesh % nodes(n1) % y
                  x1(3) = myStructuredHexMesh % nodes(n1) % z
                  
                  x2(1) = myStructuredHexMesh % nodes(n2) % x
                  x2(2) = myStructuredHexMesh % nodes(n2) % y
                  x2(3) = myStructuredHexMesh % nodes(n2) % z
                  
                  x3(1) = myStructuredHexMesh % nodes(n3) % x
                  x3(2) = myStructuredHexMesh % nodes(n3) % y
                  x3(3) = myStructuredHexMesh % nodes(n3) % z
                  
                  x4(1) = myStructuredHexMesh % nodes(n4) % x
                  x4(2) = myStructuredHexMesh % nodes(n4) % y
                  x4(3) = myStructuredHexMesh % nodes(n4) % z
                    
                  DO j = 0, gPolyDeg
                     DO i = 0, gPolyDeg
                     ! Transfinite inerpolation with linear blending is used to construct the face
                        c1 = ( HALF*(x2-x1)*(ONE+s(i)) + x1 ) 
                        c2 = ( HALF*(x3-x4)*(ONE+s(i)) + x4 )
                        xc(i,j,1:nDims) = HALF*(c2-c1)*(ONE+s(j)) + c1
                     ENDDO
                  ENDDO
                  CALL boundSurfs(iSide) % Reset( xc ) 

               ENDDO
               CALL myStructuredHexMesh % geom(iEl) % GenerateMesh( interp, boundSurfs )
               CALL myStructuredHexMesh % geom(iEl) % GenerateMetrics( interp )

            ENDDO
         ENDDO
      ENDDO ! iEl, cycle over the elements

      CALL myStructuredHexMesh % ConstructFaces( )
      nFaces = myStructuredHexMesh % nFaces
      PRINT*, 'nFaces    : ', nFaces

      CALL myStructuredHexMesh % ConstructElementNeighbors( )
      
      ! Clear up memory
      DEALLOCATE( s, xc )

      DO iSide = 1, nHexFaces
         CALL boundSurfs(iSide) % Trash( )
      ENDDO
  
 END SUBROUTINE LoadDefaultMesh_StructuredHexMesh
!
 SUBROUTINE LoadDoublyPeriodicMesh_StructuredHexMesh( myStructuredHexMesh, interp, nXelem, nYelem, nZelem  )

   IMPLICIT NONE
   CLASS( StructuredHexMesh ), INTENT(inout)  :: myStructuredHexMesh
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda ), INTENT(in)    :: interp
#else
   TYPE( Lagrange ), INTENT(in)     :: interp
#endif
   INTEGER, INTENT(in)              :: nXelem, nYelem, nZelem
   ! LOCAL
   TYPE( Surface ) :: boundSurfs(1:nHexFaces)

   REAL(prec) :: x, y, z, dxElem, dyElem, dzElem
   REAL(prec) :: x1(1:nDims), x2(1:nDims), x3(1:nDims), x4(1:nDims)
   REAL(prec) :: c1(1:nDims), c2(1:nDims)
   REAL(prec), ALLOCATABLE :: xc(:,:,:), s(:)

   INTEGER :: nNodes, nElems, nFaces, gPolyDeg
   INTEGER :: nodes(1:nHexNodes)
   INTEGER :: n1, n2, n3, n4
   INTEGER :: iNode, iEl, iSide, iX, iY, iZ, i, j

      
      dxElem = ONE/nXElem
      dyElem = ONE/nYElem
      dzElem = ONE/nZElem
      
      ! ** "Hard-wired" values for a structured mesh with no holes ** !
      nNodes   = (nXElem+1)*(nYElem+1)*(nZElem+1)
      nElems   = (nXElem)*(nYElem)*(nZElem)
      nFaces   = (nXElem)*(nYElem)*(nZElem+1) + (nXElem)*(nZElem)*(nYElem+1) + (nYElem)*(nZElem)*(nXElem+1)
      gPolyDeg = 1
      ! ************************************************************************* !

      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'gPolyDeg  : ', gPolyDeg

      ! Generate the chebyshev points of order gPolyDeg
      ! These are the points used to define the parametric
      ! curves for the element boundaries

      ALLOCATE( s(0:gPolyDeg), xc(0:gPolyDeg,0:gPolyDeg,1:3) )
      s(0) = -ONE
      s(1) = ONE
      
      ! ---- Initialize the mesh (empty) ---- !
      CALL myStructuredHexMesh % Initialize( nNodes, nElems, nFaces, interp % N ) 
      
      
      ! ---- Read in the corner nodes ---- !
      DO iZ = 1, nZElem + 1
         z = dZElem*(REAL(iZ-1,prec))
         DO iY = 1, nYElem+1
            y = dYElem*(REAL(iY-1,prec))
            DO iX = 1, nXElem+1
               iNode = iX + (iY-1)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)
               x = dXElem*(REAL(iX-1,prec))
               myStructuredHexMesh % nodes(iNode) % x = x
               myStructuredHexMesh % nodes(iNode) % y = y
               myStructuredHexMesh % nodes(iNode) % z = z
               myStructuredHexMesh % nodes(iNode) % nodeID = iNode
            ENDDO
         ENDDO
      ENDDO
  
      ! Do the element information
      xc = ZERO
      ! Do the initial build for the parametric surfaces
      DO iSide = 1, nHexFaces
         CALL boundSurfs(iSide) % Build( xc, s, gPolyDeg, nDims ) 
      ENDDO
   
      DO iZ = 1, nZElem
         DO iY = 1, nYElem
            DO iX = 1, nXElem

               iEl = iX + (iY-1)*(nXElem) + (iZ-1)*(nXElem)*(nYElem)
               ! Calculate the global node IDs for this element.
               nodes(1) = iX + (iY-1)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)    ! Southwest
               nodes(2) = iX + 1 + (iY-1)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)! SouthEast
               nodes(3) = iX + 1 + (iY)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)  ! NorthEast
               nodes(4) = iX + (iY)*(nXElem+1) + (iZ-1)*(nXElem+1)*(nYElem+1)      ! NorthWest
               nodes(5) = iX + (iY-1)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)      ! Southwest
               nodes(6) = iX + 1 + (iY-1)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)  ! SouthEast
               nodes(7) = iX + 1 + (iY)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)    ! NorthEast
               nodes(8) = iX + (iY)*(nXElem+1) + (iZ)*(nXElem+1)*(nYElem+1)        ! NorthWest

               myStructuredHexMesh % elements(iEl) % nodeIDs = nodes
         
               DO iSide = 1, nHexFaces ! Loop over the sides of the quads

                  ! To build the current face, we construct a plane that passes through
                  ! the four corner nodes. Here, we grab the global node ID's for the four
                  ! corner nodes.
                  n1 = nodes( myStructuredHexMesh % faceMap(1,iSide) )
                  n2 = nodes( myStructuredHexMesh % faceMap(2,iSide) )
                  n3 = nodes( myStructuredHexMesh % faceMap(3,iSide) ) 
                  n4 = nodes( myStructuredHexMesh % faceMap(4,iSide) ) 

                  x1(1) = myStructuredHexMesh % nodes(n1) % x
                  x1(2) = myStructuredHexMesh % nodes(n1) % y
                  x1(3) = myStructuredHexMesh % nodes(n1) % z
                  
                  x2(1) = myStructuredHexMesh % nodes(n2) % x
                  x2(2) = myStructuredHexMesh % nodes(n2) % y
                  x2(3) = myStructuredHexMesh % nodes(n2) % z
                  
                  x3(1) = myStructuredHexMesh % nodes(n3) % x
                  x3(2) = myStructuredHexMesh % nodes(n3) % y
                  x3(3) = myStructuredHexMesh % nodes(n3) % z
                  
                  x4(1) = myStructuredHexMesh % nodes(n4) % x
                  x4(2) = myStructuredHexMesh % nodes(n4) % y
                  x4(3) = myStructuredHexMesh % nodes(n4) % z
                    
                  DO j = 0, gPolyDeg
                     DO i = 0, gPolyDeg
                     ! Transfinite inerpolation with linear blending is used to construct the face
                        c1 = ( HALF*(x2-x1)*(ONE+s(i)) + x1 ) 
                        c2 = ( HALF*(x3-x4)*(ONE+s(i)) + x4 )
                        xc(i,j,1:nDims) = HALF*(c2-c1)*(ONE+s(j)) + c1
                     ENDDO
                  ENDDO
                  CALL boundSurfs(iSide) % Reset( xc ) 

               ENDDO
               CALL myStructuredHexMesh % geom(iEl) % GenerateMesh( interp, boundSurfs )
               CALL myStructuredHexMesh % geom(iEl) % GenerateMetrics( interp )

            ENDDO
         ENDDO
      ENDDO ! iEl, cycle over the elements

      CALL myStructuredHexMesh % ConstructDoublyPeriodicFaces( nXElem, nYElem, nZElem )
      nFaces = myStructuredHexMesh % nFaces
      PRINT*, 'nFaces    : ', nFaces

      CALL myStructuredHexMesh % ConstructElementNeighbors( )
      
      ! Clear up memory
      DEALLOCATE( s, xc )

      DO iSide = 1, nHexFaces
         CALL boundSurfs(iSide) % Trash( )
      ENDDO
  
 END SUBROUTINE LoadDoublyPeriodicMesh_StructuredHexMesh
!
!
!==================================================================================================!
!--------------------------------- Mesh File I/O Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup StructuredHexMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R ReadPeaceMeshFile 
! 
!> \fn ReadPeaceMeshFile_StructuredHexMesh 
!! Reads a PeaceMesh file and constructs the StructuredHexMesh data structure.
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(StructuredHexMesh)    :: this <BR>
!! <B>CHARACTER</B>        :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ReadPeaceMeshFile( filename ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myStructuredHexMesh <td> StructuredHexMesh <td>
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the PeaceMesh file
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ReadPeaceMeshFile_StructuredHexMesh( myStructuredHexMesh, filename )

   IMPLICIT NONE
   CLASS( StructuredHexMesh ), INTENT(out)   :: myStructuredHexMesh
   CHARACTER(*), INTENT(in)        :: filename
   ! LOCAL
   INTEGER :: nNodes, nElems, nFaces, N
   INTEGER :: iFace, iNode, iEl
   INTEGER :: fUnit, k, i, j, l, row, col


      PRINT*, 'Mesh File : '//TRIM( filename )//'.pc.mesh'
      
      ! Get a new file unit
      OPEN( UNIT    = NEWUNIT(fUnit), &
            FILE    = TRIM( filename )//'.pc.mesh', &
            FORM    = 'UNFORMATTED',&
            STATUS  = 'OLD', &
            ACCESS  = 'DIRECT', &
            CONVERT = 'BIG_ENDIAN', &
            RECL    = SIZEOF(nNodes) ) ! How to do variable record length

      ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
      k = 1
      READ( fUnit, rec=k )nNodes
      k = k+1
      READ( fUnit, rec=k )nElems
      k = k+1
      READ( fUnit, rec=k )nFaces
      k = k+1
      READ( fUnit, rec=k )N
      k = k+1

      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'nFaces    : ', nFaces
      PRINT*, 'N         : ', N
      ! ---- Initialize the quadrature mesh (empty) ---- !
      CALL myStructuredHexMesh % Initialize( nNodes, nElems, nFaces, N ) 
      
      ! ---- Read in the element connectivity ---- !
      DO iEl = 1, nElems
         READ( fUnit, rec=k ) myStructuredHexMesh % elements(iEl) % elementID 
         k = k+1
         DO i = 1, nHexNodes
            READ( fUnit, rec=k ) myStructuredHexMesh % elements(iEl) % nodeIDs(i)
            k = k+1
         ENDDO
      ENDDO 
      
      ! ---- Read in the face information ---- !

      DO iFace = 1, nFaces
         READ( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % faceID
         k = k+1
         READ( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % boundaryID
         k = k+1
         DO i = 1, 4
            READ( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % nodeIDs(i)
            k = k+1
         ENDDO
         DO i = 1, 2
            READ( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % elementIDs(i)
            k = k+1
            READ( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % elementSides(i)
            k = k+1
         ENDDO
         READ( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % iStart
         k = k+1
         READ( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % iInc
         k = k+1
         READ( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % jStart
         k = k+1
         READ( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % jInc
         k = k+1
         READ( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % swapDimensions
         k = k+1
      ENDDO
      
      CLOSE( fUnit )
      ! Get a new file unit
      OPEN( UNIT    = NEWUNIT(fUnit), &
            FILE    = TRIM( filename )//'.pc.geom', &
            FORM    = 'UNFORMATTED',&
            STATUS  = 'OLD', &
            ACCESS  = 'DIRECT', &
            CONVERT = 'BIG_ENDIAN', &
            RECL    = prec ) ! How to do variable record length
      
      ! ---- Read in the corner nodes ---- !
      k = 1
      DO iNode = 1, nNodes  ! Loop over the nodes in the file
         READ( fUnit, rec=k ) myStructuredHexMesh % nodes(iNode) % x
         k = k+1
         READ( fUnit, rec=k ) myStructuredHexMesh % nodes(iNode) % y
         k = k+1
         READ( fUnit, rec=k ) myStructuredHexMesh % nodes(iNode) % z
         k = k+1
      ENDDO

      ! ---- Read in the element information ---- !
      DO iEl = 1, nElems
         DO l = 0, N
            DO j = 0, N
               DO i = 0, N
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % x(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % y(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % z(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dxds(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dxdp(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dxdq(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dyds(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dydp(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dydq(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dzds(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dzdp(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dzdq(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % J(i,j,l)
                  k = k+1
                  DO col = 1, 3
                     DO row = 1, 3
                        READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % Ja(i,j,l,row,col)
                        k = k+1
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO l = 1, nHexFaces
            DO j = 0, N
               DO i = 0, N
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % xBound(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % yBound(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % zBound(i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % nHat(1,i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % nHat(2,i,j,l)
                  k = k+1
                  READ( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % nHat(3,i,j,l)
                  k = k+1
               ENDDO
            ENDDO
         ENDDO
      ENDDO 

      CLOSE( fUnit )
     
      CALL myStructuredHexMesh % ConstructElementNeighbors( )
  
 END SUBROUTINE ReadPeaceMeshFile_StructuredHexMesh
!
!> \addtogroup StructuredHexMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R WritePeaceMeshFile 
! 
!> \fn WritePeaceMeshFile_StructuredHexMesh 
!! Writes a PeaceMesh file using the StructuredHexMesh data structure.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(StructuredHexMesh)    :: this <BR>
!! <B>CHARACTER</B>        :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % WritePeaceMeshFile( filename ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myStructuredHexMesh <td> StructuredHexMesh <td>
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the PeaceMesh file
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
SUBROUTINE WritePeaceMeshFile_StructuredHexMesh( myStructuredHexMesh, filename )

   IMPLICIT NONE
   CLASS( StructuredHexMesh ), INTENT(in)   :: myStructuredHexMesh
   CHARACTER(*), INTENT(in)       :: filename
   ! LOCAL
   INTEGER :: nNodes, nElems, nFaces, N
   INTEGER :: iFace, iNode, iEl
   INTEGER :: fUnit, k, i, j, l, row, col


      PRINT*, 'Mesh File : '//TRIM( filename )//'.pc.mesh'
      nNodes = 1
      ! Get a new file unit
      OPEN( UNIT    = NEWUNIT(fUnit), &
            FILE    = TRIM( filename )//'.pc.mesh', &
            FORM    = 'UNFORMATTED',&
            STATUS  = 'REPLACE', &
            ACCESS  = 'DIRECT', &
            CONVERT = 'BIG_ENDIAN', &
            RECL    = SIZEOF(nNodes) ) ! How to do variable record length

      
      ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
      nNodes = myStructuredHexMesh % nNodes
      nElems = myStructuredHexMesh % nElems
      nFaces = myStructuredHexMesh % nFaces
      N      = myStructuredHexMesh % geom(1) % N
      k = 1
      WRITE( fUnit, rec=k )nNodes
      k = k+1
      WRITE( fUnit, rec=k )nElems
      k = k+1
      WRITE( fUnit, rec=k )nFaces
      k = k+1
      WRITE( fUnit, rec=k )N
      k = k+1

      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'nFaces    : ', nFaces
      PRINT*, 'N         : ', N

      ! ---- Read in the element connectivity ---- !
      DO iEl = 1, nElems
         WRITE( fUnit, rec=k ) myStructuredHexMesh % elements(iEl) % elementID 
         k = k+1
         DO i = 1, nHexNodes
            WRITE( fUnit, rec=k ) myStructuredHexMesh % elements(iEl) % nodeIDs(i)
            k = k+1
         ENDDO
      ENDDO 
      
      ! ---- Read in the face information ---- !

      DO iFace = 1, nFaces
         WRITE( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % faceID
         k = k+1
         WRITE( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % boundaryID
         k = k+1
         DO i = 1, 4
            WRITE( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % nodeIDs(i)
            k = k+1
         ENDDO
         DO i = 1, 2
            WRITE( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % elementIDs(i)
            k = k+1
            WRITE( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % elementSides(i)
            k = k+1
         ENDDO
         WRITE( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % iStart
         k = k+1
         WRITE( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % iInc
         k = k+1
         WRITE( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % jStart
         k = k+1
         WRITE( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % jInc
         k = k+1
         WRITE( fUnit, rec=k ) myStructuredHexMesh % faces(iFace) % swapDimensions
         k = k+1
      ENDDO
      
      CLOSE( fUnit )
      ! Get a new file unit
      OPEN( UNIT    = NEWUNIT(fUnit), &
            FILE    = TRIM( filename )//'.pc.geom', &
            FORM    = 'UNFORMATTED',&
            STATUS  = 'REPLACE', &
            ACCESS  = 'DIRECT', &
            CONVERT = 'BIG_ENDIAN', &
            RECL    = prec ) ! How to do variable record length
      
      ! ---- Read in the corner nodes ---- !
      k = 1
      DO iNode = 1, nNodes  ! Loop over the nodes in the file
         WRITE( fUnit, rec=k ) myStructuredHexMesh % nodes(iNode) % x
         k = k+1
         WRITE( fUnit, rec=k ) myStructuredHexMesh % nodes(iNode) % y
         k = k+1
         WRITE( fUnit, rec=k ) myStructuredHexMesh % nodes(iNode) % z
         k = k+1
      ENDDO

      ! ---- Read in the element information ---- !
      DO iEl = 1, nElems
         DO l = 0, N
            DO j = 0, N
               DO i = 0, N
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % x(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % y(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % z(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dxds(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dxdp(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dxdq(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dyds(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dydp(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dydq(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dzds(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dzdp(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % dzdq(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % J(i,j,l)
                  k = k+1
                  DO col = 1, 3
                     DO row = 1, 3
                        WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % Ja(i,j,l,row,col)
                        k = k+1
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO l = 1, nHexFaces
            DO j = 0, N
               DO i = 0, N
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % xBound(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % yBound(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % zBound(i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % nHat(1,i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % nHat(2,i,j,l)
                  k = k+1
                  WRITE( fUnit, rec=k ) myStructuredHexMesh % geom(iEl) % nHat(3,i,j,l)
                  k = k+1
               ENDDO
            ENDDO
         ENDDO
      ENDDO 
      
      CLOSE( fUnit )
     
 END SUBROUTINE WritePeaceMeshFile_StructuredHexMesh
!
!> \addtogroup StructuredHexMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R ReadUCDMeshFile 
! 
!> \fn ReadUCDMeshFile_StructuredHexMesh 
!! Reads in a ucd unstructured mesh file. It is assumed that the nodes and elements are numbered
!! between 1 to nNodes and 1 to nElems, respectively. To ensure this in Trellis or CuBit, type
!! "compress all" before exporting ucd file.
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(StructuredHexMesh)    :: this <BR>
!! <B>CHARACTER</B>        :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ReadUCDMeshFile( filename ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myStructuredHexMesh <td> StructuredHexMesh <td>
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the PeaceMesh file
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ReadUCDMeshFile_StructuredHexMesh( myStructuredHexMesh, interp, filename )

   IMPLICIT NONE
   CLASS( StructuredHexMesh ), INTENT(out)   :: myStructuredHexMesh
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda ), INTENT(in)    :: interp
#else
   TYPE( Lagrange ), INTENT(in)     :: interp
#endif
   CHARACTER(*), INTENT(in)        :: filename
   ! LOCAL
   CHARACTER(60) :: longDummy
   CHARACTER(3)  :: shortDummy
   INTEGER :: nNodes, nElems, nFaces
   INTEGER :: iFace, iNode, iEl, iSide
   INTEGER :: fUnit, k, i, j, l, row, col, n1, n2, n3, n4
   REAL(prec) :: xc(0:1,0:1,1:3), s(0:1)
   REAL(prec) :: x1(1:3), x2(1:3), x3(1:3), x4(1:3), c1(1:3), c2(1:3)
   TYPE( Surface ) :: boundSurfs(1:nHexFaces)


      PRINT*, 'Mesh File : '//TRIM( filename )
      
      ! Get a new file unit
      OPEN( UNIT    = NEWUNIT(fUnit), &
            FILE    = TRIM( filename ), &
            FORM    = 'FORMATTED',&
            STATUS  = 'OLD', &
            ACCESS  = 'SEQUENTIAL' ) 
 
      READ( fUnit, * ) longDummy
      WRITE( *, * ) longDummy
      
      READ( fUnit, * ) longDummy
      WRITE( *, * ) longDummy
      
      READ( fUnit, * ) longDummy
      WRITE( *, * ) longDummy
      
      READ( fUnit, * ) nNodes, nElems, i, j, l ! i,j, and l are just fillers for now.
      
      
      ! ---- Gather the number of nodes, number of elements, and number of edges ---- !
      ! Default nFaces = 1 for initial build
      nFaces = 1

      PRINT*, 'nNodes    : ', nNodes
      PRINT*, 'nElems    : ', nElems
      PRINT*, 'N         : ', interp % N
      ! ---- Initialize the quadrature mesh (empty) ---- !
      CALL myStructuredHexMesh % Initialize( nNodes, nElems, nFaces, interp % N ) 
      
      ! Read in the corner node positions
      DO iNode = 1, nNodes
         READ( fUnit, * ) myStructuredHexMesh % nodes(iNode) % nodeID, &
                          myStructuredHexMesh % nodes(iNode) % x, &
                          myStructuredHexMesh % nodes(iNode) % y, &
                          myStructuredHexMesh % nodes(iNode) % z
      ENDDO
      
      ! ---- Read in the element connectivity ---- !
      DO iEl = 1, nElems
         myStructuredHexMesh % elements(iEl) % elementID = iEl
         READ( fUnit, * ) i, j, shortDummy, &
                          myStructuredHexMesh % elements(iEl) % nodeIDs(1:8)
      ENDDO 
      CLOSE( fUnit )
      
      CALL myStructuredHexMesh % ConstructFaces( )
      CALL myStructuredHexMesh % ConstructElementNeighbors( )
      
     ! Do the element information
      s(0) = -1.0_prec
      s(1) = 1.0_prec
      xc   = 0.0_prec
      ! Do the initial build for the parametric surfaces
      DO iSide = 1, nHexFaces
         CALL boundSurfs(iSide) % Build( xc, s, 1, nDims ) 
      ENDDO
      
      
      ! For now, we assume tri-linear mappings so that the geometry can be constructed
      ! using only the corner node locations. These corner nodes will be used to generate 
      ! the bounding surfaces for each element, and transfinite interpolation with linear
      ! blending is used to construct the internal mesh geometry.
      DO iEl = 1, nElems
      
          DO iSide = 1, nHexFaces ! Loop over the sides of the quads

             ! To build the current face, we construct a plane that passes through
             ! the four corner nodes. Here, we grab the global node ID's for the four
             ! corner nodes.
             n1 = myStructuredHexMesh % elements(iEl) % nodeIDs( myStructuredHexMesh % faceMap(1,iSide) )
             n2 = myStructuredHexMesh % elements(iEl) % nodeIDs( myStructuredHexMesh % faceMap(2,iSide) )
             n3 = myStructuredHexMesh % elements(iEl) % nodeIDs( myStructuredHexMesh % faceMap(3,iSide) ) 
             n4 = myStructuredHexMesh % elements(iEl) % nodeIDs( myStructuredHexMesh % faceMap(4,iSide) ) 

             x1(1) = myStructuredHexMesh % nodes(n1) % x
             x1(2) = myStructuredHexMesh % nodes(n1) % y
             x1(3) = myStructuredHexMesh % nodes(n1) % z
               
             x2(1) = myStructuredHexMesh % nodes(n2) % x
             x2(2) = myStructuredHexMesh % nodes(n2) % y
             x2(3) = myStructuredHexMesh % nodes(n2) % z
               
             x3(1) = myStructuredHexMesh % nodes(n3) % x
             x3(2) = myStructuredHexMesh % nodes(n3) % y
             x3(3) = myStructuredHexMesh % nodes(n3) % z
               
             x4(1) = myStructuredHexMesh % nodes(n4) % x
             x4(2) = myStructuredHexMesh % nodes(n4) % y
             x4(3) = myStructuredHexMesh % nodes(n4) % z
                 
             DO j = 0, 1
                DO i = 0, 1
                   ! Transfinite inerpolation with linear blending is used to construct the face
                   c1 = ( HALF*(x2-x1)*(ONE+s(i)) + x1 ) 
                   c2 = ( HALF*(x3-x4)*(ONE+s(i)) + x4 )
                   xc(i,j,1:nDims) = HALF*(c2-c1)*(ONE+s(j)) + c1
                ENDDO
             ENDDO
             CALL boundSurfs(iSide) % Reset( xc ) 

          ENDDO
          CALL myStructuredHexMesh % geom(iEl) % GenerateMesh( interp, boundSurfs )
          CALL myStructuredHexMesh % geom(iEl) % GenerateMetrics( interp )
      
      ENDDO
      
      
      DO iSide = 1, nHexFaces
         CALL boundSurfs(iSide) % Trash( ) 
      ENDDO
  
 END SUBROUTINE ReadUCDMeshFile_StructuredHexMesh
!
!> \addtogroup StructuredHexMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R WriteTecplot
! 
!> \fn WriteTecplot_StructuredHexMesh  
!! Writes a tecplot file of the mesh geometry. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(StructuredHexMesh) :: this <BR>
!! <B>CHARACTER</B>      :: filename <BR>
!!         .... <BR>
!!     ! To write a file with a specified name <BR>
!!     <B>CALL</B> this % WriteTecplot( filename ) <BR>
!!     ! Or, the file "mesh.tec" will be used with the following calling sequence <BR>
!!     <B>CALL</B> this % WriteTecplot( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myStructuredHexMesh <td> StructuredHexMesh <td>
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the tecplot file
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE WriteTecplot_Hexmesh( myStructuredHexMesh, filename )

  IMPLICIT NONE
  CLASS(StructuredHexMesh), INTENT(inout)     :: myStructuredHexMesh
  CHARACTER(*), INTENT(in), OPTIONAL :: filename  
  ! Local
  INTEGER :: iS, iP, iQ, N, iEl, fUnit, eID
  CHARACTER(7) :: zoneID

      N = myStructuredHexMesh % geom(1) % N

    IF( PRESENT(filename) )THEN
    
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= TRIM(filename)//'.tec', &
             FORM='formatted')
    ELSE
    
       OPEN( UNIT=NEWUNIT(fUnit), &
             FILE= 'mesh.tec', &
             FORM='formatted')

    ENDIF
    
    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "Jacobian", "dxds", "dxdp", "dxdq", "dyds", "dydp", "dydq", "dzds", "dzdp", "dzdq" '


    DO iEl = 1, myStructuredHexMesh % nElems

       eID = myStructuredHexMesh % elements(iEl) % elementID
       WRITE(zoneID,'(I7.7)') eID
       WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',N+1,', J=',N+1,', K=',N+1,',F=POINT'

       DO iQ = 0, N
          DO iP = 0, N
             DO iS = 0,N

                WRITE(fUnit,*)  myStructuredHexMesh % geom(iEl) % x(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % y(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % z(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % J(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % dxds(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % dxdp(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % dxdq(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % dyds(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % dydp(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % dydq(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % dzds(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % dzdp(iS,iP,iQ), &
                               myStructuredHexMesh % geom(iEl) % dzdq(iS,iP,iQ)
             ENDDO
          ENDDO
      ENDDO

    ENDDO
    
    CLOSE(UNIT=fUnit)

    RETURN

 END SUBROUTINE WriteTecplot_Hexmesh
!
!> \addtogroup StructuredHexMesh_Class 
!! @{ 
! ================================================================================================ !
! S/R WriteMaterialTecplot
! 
!> \fn WriteMaterialTecplot_StructuredHexMesh  
!! Writes a tecplot file of the mesh geometry and a given "material" field. 
!! 
!! When performing domain decomposition (e.g. with the DecomposeStructuredHexMesh.f90 program), the 
!! decomposition can be visualized by assigning each element a number that corresponds to the 
!! process ID it has been assigned to. This routine takes in an array of "material ID's" that
!! identify each element.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(StructuredHexMesh) :: this <BR>
!! <B>REAL</B>(prec)     :: pID(1:this % nElems) <BR>
!! <B>CHARACTER</B>      :: filename <BR>
!!         .... <BR>
!!     ! To write a file with a specified name <BR>
!!     <B>CALL</B> this % WriteMaterialTecplot( pID, filename ) <BR>
!!     ! Or, the file "mesh.tec" will be used with the following calling sequence <BR>
!!     <B>CALL</B> this % WriteMaterialTecplot( pID ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myStructuredHexMesh <td> StructuredHexMesh <td>
!!   <tr> <td> in <th> materialIDs(1:myStructuredHexMesh % nElems) <td> REAL(prec) <td> 
!!                     Array of values that are assigned to each element.
!!   <tr> <td> in <th> filename <td> CHARACTER <td> Base-name of the MaterialTecplot file
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE WriteMaterialTecplot_Hexmesh( myStructuredHexMesh, materialIDs, filename )

   IMPLICIT NONE
   CLASS(StructuredHexMesh), INTENT(inout)      :: myStructuredHexMesh
   REAL(prec), INTENT(in)             :: materialIDs(1:myStructuredHexMesh % nElems)
   CHARACTER(*), INTENT(in), OPTIONAL :: filename  
   ! Local
   INTEGER :: iS, iP, iQ, N, iEl,fUnit
   CHARACTER(7) :: zoneID

      N = myStructuredHexMesh % geom(1) % N

      IF( PRESENT(filename) )THEN
         OPEN( UNIT=NEWUNIT(fUnit), &
               FILE= TRIM(filename)//'.tec', &
               FORM='formatted')
      ELSE
         OPEN( UNIT=NEWUNIT(fUnit), &
               FILE= 'mesh.tec', &
               FORM='formatted')
      ENDIF
           
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "materialIDs" '
   
      DO iEl = 1, myStructuredHexMesh % nElems
   
         WRITE(zoneID,'(I7.7)') myStructuredHexMesh % elements(iEl) % elementID
         WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',N+1,', J=',N+1,', K=',N+1,',F=POINT'
   
         DO iQ = 0, N
            DO iP = 0, N
               DO iS = 0, N
                  WRITE(fUnit,*) myStructuredHexMesh % geom(iEl) % x(iS,iP,iQ), &
                                 myStructuredHexMesh % geom(iEl) % y(iS,iP,iQ), &
                                 myStructuredHexMesh % geom(iEl) % z(iS,iP,iQ), &
                                 materialIDs(iEl)
               ENDDO
            ENDDO
         ENDDO
   
      ENDDO
       
      CLOSE(UNIT=fUnit)
      RETURN

 END SUBROUTINE WriteMaterialTecplot_Hexmesh
!
END MODULE StructuredHexMesh_Class

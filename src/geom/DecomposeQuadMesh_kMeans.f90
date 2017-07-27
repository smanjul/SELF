PROGRAM DecomposeQuadMesh_kMeans

 ! src/common/
 USE ModelPrecision
 USE ConstantsDictionary
 USE CommonRoutines
 USE Timing
 ! src/nodal/
 USE NodalStorage_Class
 ! src/geom/
 USE QuadMesh_Class
 USE Params_Class
 ! src/boundary/
 USE BoundaryCommunicator_Class
 
 ! To make use of the METIS data types and intrinsics
 !USE metis
  
 IMPLICIT NONE
 
 INTEGER, PARAMETER           :: nTrials = 50 ! This is the number of k-means trials that are run
                                               ! with random initial cluster centers
 INTEGER, PARAMETER           :: iterMax = 500 ! Maximum number of iterates for determining the cluster centers
 REAL(prec), PARAMETER        :: tolerance = 10.0_prec**(-5)
 
 
 TYPE( NodalStorage )         :: nodal
 TYPE( QuadMesh )             :: mesh, procMesh
 TYPE( RunParams )            :: params
 TYPE( MultiTimers )          :: timers
 TYPE( BoundaryCommunicator ) :: bcom
 INTEGER                      :: procID, rnInt, fUnit, iTrial, iter
 INTEGER                      :: nElems, nProc, pID
 INTEGER                      :: iEl, jEl, iSide, iNode, nAdj, nID, elID, nBe, nMPI
 INTEGER                      :: n1, n2, ng1, ng2, m1, m2, mg1, mg2, e2, iEdge, jEdge, localID
 INTEGER, ALLOCATABLE         :: globalToLocal(:,:), nElPerProc(:), nodeLogic(:,:), nNodePerProc(:)
 INTEGER, ALLOCATABLE         :: globalToLocalNode(:,:), elementToCluster(:,:)
 REAL(prec), ALLOCATABLE      :: materials(:), elementCenter(:,:), centroidDistances(:,:,:)
 REAL(prec), ALLOCATABLE      :: centroidCenters(:,:), newCC(:,:), cD(:), nElemInCluster(:,:)
 REAL(prec), ALLOCATABLE      :: centroidDistMag(:)
 REAL(prec)                   :: dc, rn
 CHARACTER(4)                 :: pIDChar
 
      ! Read in the parameters
      CALL params % Build( )
      CALL timers % Build( )
      
      CALL timers % AddTimer( 'Spectral Partition Time', 1 )
      ! Build an interpolant
      CALL nodal % Build( N = params % polyDeg, &
                          nPlot = params % nPlot, &
                          quadrature = GAUSS,  &
                          approxForm = DG )

      ! Build the Geometry
      IF( params % SpecMeshFile == NADA )THEN
         CALL mesh % LoadDefaultMesh( nodal % interp, params % nXElem, params % nYElem )
      ELSE
         CALL mesh % ReadSpecMeshFile( nodal % interp, params % SpecMeshFile )
      ENDIF

      nElems = mesh % nElems
      nProc  = params % nProc
      ALLOCATE( materials(1:nElems), &
                globalToLocal(1:nElems,1:2), &
                nElPerProc(0:nProc-1), &
                nodeLogic(1:mesh % nNodes, 0:nProc-1), &
                nNodePerProc(0:nProc-1), &
                globalToLocalNode(1:mesh % nNodes, 0:nProc-1), &
                elementCenter(1:2,1:mesh % nElems), &
                elementToCluster(1:mesh % nElems,1:nTrials), &
                centroidDistances(1:mesh % nElems, 1:nProc,1:nTrials), &
                centroidDistMag(1:nTrials), &
                centroidCenters(1:2,1:nProc), &
                newCC(1:2,1:nProc), &
                cD(1:nProc), &
                nElemInCluster(1:nProc,1:nTrials) )
                
      materials         = 0.0_prec
      globalToLocal     = 0
      nElPerProc        = 0
      nodeLogic         = 0
      nNodePerProc      = 0
      elementCenter     = 0.0_prec
      elementToCluster  = 0
      centroidDistances = 0.0_prec
      centroidCenters   = 0.0_prec
      newCC             = 0.0_prec
      cD                = 0.0_prec
      nElemInCluster    = 0.0_prec
      centroidDistMag   = 0.0_prec
   
      CALL timers % StartTimer( 1 )
      
      ! Initialize element centers list
      ! The center physical location of the elements is used in the k-means clustering algorithm
      ! to group nearby elements
      DO iEl = 1, mesh % nElems
         elementCenter(1:2,iEl) = mesh % geom(iEl) % CalculateLocation( nodal % interp, &
                                                                        (/ 0.0_prec, 0.0_prec /) )
      ENDDO
                                                                        
      ! Loop over the random centroid distributions
      ! k-means must be initialized with a set of centroids. Here, we conduct "nTrials" experiments
      ! where the initial centroids are randomly selected from the element list.
      OPEN( UNIT   = NewUnit(fUnit), &
            FILE   = 'kMeans_Changes.txt', &
            STATUS = 'REPLACE', &
            ACTION = 'WRITE', &
            FORM   = 'FORMATTED' ) 
            
      DO iTrial = 1, nTrials
         
         ! Set the centroid centers here.
         
         ! >> In each loop
         ! >>    Loop over centroids until centroids do not change
         DO pID = 1, nProc
            CALL RANDOM_NUMBER( rn )
            centroidCenters(1,pID) = rn
            CALL RANDOM_NUMBER( rn )
            centroidCenters(2,pID) = rn        
         ENDDO
        
         DO iter = 1, itermax
         
            ! Loop over centroids and calculate the distances between the elements and the centroid
            ! centers.
            DO pID = 1, nProc
               DO iEl = 1, mesh % nElems
                  ! >>    >> Sort elements into cluster lists based on distance from each centroid
                  centroidDistances(iEl, pID, iTrial) = sqrt( ( elementCenter(1,iEl) - centroidCenters(1,pID) )**2 + &
                                                              ( elementCenter(2,iEl) - centroidCenters(2,pID) )**2 )
               ENDDO
            ENDDO
            
            ! Find the cluster center that is closest to each element and 
            ! recalculate the cluster centers
            newCC = 0.0_prec
            nElemInCluster(1:nProc,iTrial) = 0.0_prec
            DO iEl = 1, mesh % nElems
            
               cD  = centroidDistances(iEl,1:nProc, iTrial)
               pID = MINLOC( cD, 1 ) ! Find the cluster index that corresponds to the minimum distance
                                  ! between element and cluster center
               elementToCluster(iEl,iTrial) = pID
               
               newCC(1:2,pID) = newCC(1:2,pID) + elementCenter(1:2,iEl)
               nElemInCluster(pID,iTrial) = nElemInCluster(pID,iTrial) + 1.0_prec
               
            ENDDO
            
            ! Divide by the number of elements in each cluster to obtain the centroid
            ! Also calculate the magnitude of the change in the cluster centers to determine
            ! if the method has converged for this trial
            dc = 0.0_prec
            DO pID = 1, nProc
               newCC(1:2,pID) = newCC(1:2,pID)/nElemInCluster(pID,iTrial)
               dc = dc +  sqrt( ( newCC(1,pID) - centroidCenters(1,pID) )**2 + &
                                ( newCC(2,pID) - centroidCenters(2,pID) )**2 )
               centroidCenters(1:2,pID) = newCC(1:2,pID) ! Update the centroid centers
            ENDDO
            
            WRITE( fUnit, * ) iTrial, dc
            
            IF( dc <= tolerance )THEN
               EXIT
            ENDIF
            
         ENDDO
         
      ENDDO
      CLOSE( fUnit )
      
      DO iTrial = 1, nTrials
      
         centroidDistMag(iTrial) = 0.0_prec
         
         DO pID = 1, nProc
            DO iEl = 1, mesh % nElems
               centroidDistMag(iTrial) = centroidDistMag(iTrial) + centroidDistances(iEl,pID,iTrial)
            ENDDO
         ENDDO
         
      ENDDO
      
      ! Pick the most ideal trial
     ! iTrial = 1 ! ** FOR NOW ** !
      iTrial = MINLOC( centroidDistMag, 1 )
      
      ! Write element-to-cluster list to file - this file can be used in developing a multigrid
      ! preconditioner.
      OPEN( UNIT   = NewUnit(fUnit), &
            FILE   = TRIM( params % PeaceMeshFile )//'.clusters', &
            FORM   = 'FORMATTED', &
            ACTION = 'WRITE', &
            STATUS = 'REPLACE' )
            
      WRITE( fUnit, * ) nProc
      DO iEl = 1, mesh % nElems
         WRITE( fUnit, * ) elementToCluster(iEl,iTrial)
      ENDDO
      
      CLOSE( fUnit )
            
      ! Now, we need to build each mesh and and the associated connecivity
      DO iEl = 1, nElems

         !procID               = partitions(iEl-1)      
         procID               = elementToCluster(iEl,iTrial) - 1
         nElPerProc(procID)   = nElPerProc(procID) + 1 ! Increment the local element ID for this process
         globalToLocal(iEl,1) = nElPerProc(procID)     ! Store the local element ID on this process
         globalToLocal(iEl,2) = procID                 ! Store the process ID
         
         DO iNode = 1, 4
            nID = mesh % elements(iEl) % nodeIDs(iNode) 
            nodeLogic(nID,procID) = 1
         ENDDO
         
      ENDDO
      DO nID = 1, mesh % nNodes
         nNodePerProc = nNodePerProc + nodeLogic(nID,:)
         globalToLocalNode(nID,:) = nNodePerProc
      ENDDO
      
      DO procID = 0, nProc-1
         PRINT*, 'Process ID :',procID, ', nElems :', nElPerProc(procID)
         PRINT*, 'Process ID :',procID, ', nNodes :', nNodePerProc(procID)
      ENDDO
      
      CALL timers % StopTimer( 1 )
      
      ! Now we generate the local mesh for each process
      DO procID = 0, nProc-1
      
         CALL procMesh % Initialize(  nNodePerProc(procID), &
                                      nElPerProc(procID), &
                                      1, params % polyDeg )
                                      
         !PRINT*, 'Setting up local nodes'
         ! The nodes on the local node mesh are assigned using the global mesh and the "nodeLogic"
         ! array
         DO nID = 1, mesh % nNodes
         
            ! We only assign node information if it is "owned" by this process
            IF( nodeLogic(nID,procID) == 1 )THEN 
               ! Obtain the local node ID from the "global-to-local" array
               localID = globalToLocalNode( nID, procID )
               ! Assign the local node attributes
               procMesh % nodes(localID) % nodeID   = nID
               procMesh % nodes(localID) % nodeType = mesh % nodes(nID) % nodeType
               procMesh % nodes(localID) % x        = mesh % nodes(nID) % x
               procMesh % nodes(localID) % y        = mesh % nodes(nID) % y
               procMesh % nodes(localID) % z        = mesh % nodes(nID) % z
            ENDIF
            
         ENDDO

         ! Set the element-node connectivity
         DO iEl = 1, mesh % nElems
            ! We only assign element information if it is "owned" by this process
            IF( globalToLocal(iEl,2) == procID )THEN
               elID = globalToLocal(iEl,1) ! Obtain the local element ID from the global-to-local array
               
               DO iNode = 1, 4
                  ! First we grab the global node ID for this element's corner node
                  nID = mesh % elements(iEl) % nodeIDs(iNode)
                  ! Then we convert this to a local node ID for this process
                  localID = globalToLocalNode( nID, procID )
                  ! And now we assign the local Node ID to the local element's corner node ID's
                  procMesh % elements(elID) % nodeIDs(iNode)  = localID
               ENDDO
               procMesh % elements(elID) % elementID = iEl
                  
               ! Set the geometry
               procMesh % geom(elID) % nHat   = mesh % geom(iEl) % nHat
               procMesh % geom(elID) % xBound = mesh % geom(iEl) % xBound
               procMesh % geom(elID) % yBound = mesh % geom(iEl) % yBound
               procMesh % geom(elID) % x      = mesh % geom(iEl) % x
               procMesh % geom(elID) % y      = mesh % geom(iEl) % y
               procMesh % geom(elID) % J      = mesh % geom(iEl) % J
               procMesh % geom(elID) % dxds   = mesh % geom(iEl) % dxds
               procMesh % geom(elID) % dxdp   = mesh % geom(iEl) % dxdp
               procMesh % geom(elID) % dyds   = mesh % geom(iEl) % dyds
               procMesh % geom(elID) % dydp   = mesh % geom(iEl) % dydp
            ENDIF
            
         ENDDO
         
         ! Set the edges
         ! First we count the local edges and fill in the edge connectivity
         
         CALL procMesh % ConstructEdges( )
         PRINT*, 'Process ID :',procID, ', nEdges :', procMesh % nEdges
         
         ! Assign the start and inc
         DO iEdge = 1, procMesh % nEdges
            n2 = procMesh % edges(iEdge) % elementSides(2)
            IF(n2 < 0)THEN
               procMesh % edges(iEdge) % start = params % polyDeg -1
               procMesh % edges(iEdge) % inc   = -1
            ELSE
               procMesh % edges(iEdge) % start = 1
               procMesh % edges(iEdge) % inc   = 1
            ENDIF
            
        ENDDO
         
         ! Now we need to cycle through the edges and find the boundary edges
         ! If we find a boundary edge (unassigned secondary element), then we need to search
         ! through the global mesh to find the matching edge to determine the appropriate 
         ! boundary condition. The boundary condition could either be a physical boundary condition,
         ! numerical boundary condition, or may need to indicate an MPI_SEND/RECV pair.
         
         ! Count the number of boundary edges
         nBe = 0
         DO iEdge = 1, procMesh % nEdges
            e2 = procMesh % edges(iEdge) % ElementIDs(2)
            IF( e2 <= 0 )THEN
               nBe = nBe + 1
            ENDIF
         ENDDO
         PRINT*, 'Process ID :',procID, ', nBEdge :', nBe
         !PRINT*, 'Setting up boundary communications'
         
         ! And initialize the boundary communicator
         CALL bCom % Initialize( nBe )
         ! Assign the boundary edge ID to the
         nBe = 0
         nMPI = 0
         DO iEdge = 1, procMesh % nEdges
         
            e2 = procMesh % edges(iEdge) % ElementIDs(2)
            IF( e2 <= 0 )THEN ! This edge is a boundary of this process's edge
            
               nBe = nBe + 1
               bCom % boundaryIDs(nBe) = iEdge
               ! Initially fill in the boundary condition as provided, and assign the external
               ! process ID to the current process ID (implying no MPI communication)
               bCom % extElemIDs(nBe) = e2
               bCom % extProcIDs(nBe) = procID 
               procMesh % edges(iEdge) % boundaryID = nBe
               
               ! Now we search through the global mesh for an identical edge
               ! Recall that the edge is identified by its terminating nodes.
               ! Here, we must use the global node ID's
               n1 = procMesh % edges(iEdge) % nodeIDs(1) ! local node ID
               n2 = procMesh % edges(iEdge) % nodeIDs(2) ! local node ID
               ng1 = min( procMesh % nodes(n1) % nodeID, procMesh % nodes(n2) % nodeID)  ! Global node ID
               ng2 = max( procMesh % nodes(n1) % nodeID, procMesh % nodes(n2) % nodeID)  ! Global node ID
               
               DO jEdge = 1, mesh % nEdges ! Cycle through the global mesh and search for an identical edge
                  ! Obtain the terminating global node IDs from the global mesh
                  m1 = mesh % edges(jEdge) % nodeIDs(1)
                  m2 = mesh % edges(jEdge) % nodeIDs(2)
                  mg1 = min( m1, m2 )
                  mg2 = max( m1, m2 )
              
                  IF( ng1 == mg1 .AND. ng2 == mg2 )THEN ! we have found an identical edge
                     
                     ! In this case, we have found an identical edge 
                     ! We need to copy over the secondary element information
                     iEl = mesh % edges(jEdge) % elementIDs(2)
                     IF( iEl > 0 )THEN ! If iEl > 0, then this edge is internal to the global mesh
                                       ! implying that an MPI communication is needed
                        nMPI = nMPI + 1
                        ! For this case, an MPI_SEND/RECV pair is needed and all secondary element
                        ! information must be copied from the global mesh
                        elID = globalToLocal(iEl, 1) ! obtain the process local element ID
                        pID  = globalToLocal(iEl, 2) ! obtain the process ID
                        
                        ! If the process ID of the global-secondary element is owned by the current 
                        ! process, then we need to assign the global primary element as the local
                        ! secondary element
                        IF( pID == procID )THEN
                           iEl = mesh % edges(jEdge) % elementIDs(1)
                           elID = globalToLocal(iEl, 1) ! obtain the process local element ID
                           pID  = globalToLocal(iEl, 2) ! obtain the process ID
                        ENDIF
                        
                     !   print*,iEl, elID, pID
                        procMesh % edges(iEdge) % edgeID = jEdge
                        procMesh % edges(iEdge) % elementIDs(2) = -elID
                        procMesh % edges(iEdge) % start = mesh % edges(jEdge) % start
                        procMesh % edges(iEdge) % inc = mesh % edges(jEdge) % inc
                        
                        ! Update the boundary communicator information
                        bCom % extElemIDs(nBe) = elID
                        bCom % extProcIDs(nBe) = pID
                        
                     ELSE
                        ! In this case, a physical or numerical boundary is found and we only 
                        ! need to assign the secondary element ID
                        procMesh % edges(iEdge) % elementIDs(2) = iEl
                     ENDIF
                     
                     ! In either case, we have identified the edge and we can exit this inner loop
                     EXIT
                     
                  ENDIF
               
               ENDDO 
               
            ENDIF ! If this is a boundary edge
            
         ENDDO ! loop over the local edges
         PRINT*, 'Process ID :',procID, ', nMPI   :', nMPI
         
         ! Now we need to write a peace-mesh file and and communicator file
         WRITE( pIDChar, '(I4.4)' ) procID
         CALL procmesh % WriteTecplot( 'mesh.'//pIDChar )
         CALL procmesh % WritePeaceMeshFile( TRIM(params % PeaceMeshFile)//'.'//pIDChar )
         CALL bCom % WritePickup( 'ExtComm.'//pIDChar )
         
         CALL procMesh % Trash( )
         CALL bCom % Trash( )
         
      ENDDO
      
      ! For visualizing the decomposition
      materials = REAL( elementToCluster(1:mesh % nElems,iTrial) - 1, prec )
      CALL mesh % WriteMaterialTecplot( materials )
      
      CALL timers % Write_MultiTimers( )

      ! Clean up memory !
      DEALLOCATE( materials, &
                  globalToLocal, &
                  nElPerProc, &
                  nodeLogic, &
                  nNodePerProc, &
                  globalToLocalNode, &
                  elementCenter, &
                  elementToCluster, &
                  centroidDistances, &
                  centroidCenters, &
                  newCC, &
                  cD, &
                  nElemInCluster )
      
   !   CALL mesh % WriteTecplot( )
      CALL mesh % Trash( )
      CALL nodal % Trash( )
      CALL timers % Trash( )

END PROGRAM DecomposeQuadMesh_kMeans

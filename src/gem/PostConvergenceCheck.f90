
!-------------------------------------------------------------------------------------------------------------
!> \file    PostConvergenceCheck.f90
!> \brief   Post-convergence verification: check all non-stable phases
!!          (both solution and pure condensed) for negative driving forces
!!          that were missed by convergence shortcuts.
!!
!!          If a missing phase is detected, force-add it to the assemblage
!!          (bypassing normal gatekeeping) so the re-convergence loop can
!!          find the true global minimum.
!-------------------------------------------------------------------------------------------------------------

subroutine PostConvergenceCheck(lConv)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    logical, intent(inout) :: lConv
    integer :: i, j, k, iFirst, iLast, iMaxDrivingForce, iBestSoln, iWeakestCon
    real(8) :: dMaxDrivingForce, dBestDrivingForce, dWeakestConMoles
    logical :: lPhaseChange, lTraceOn
    real(8), allocatable :: dBestMolFraction(:)
    character(len=32) :: cTraceEnv

    ! Env-var-gated diagnostic tracing (TC_TRACE_POSTCONV=1)
    lTraceOn = .FALSE.
    call GET_ENVIRONMENT_VARIABLE('TC_TRACE_POSTCONV', cTraceEnv)
    if (TRIM(cTraceEnv) == '1') lTraceOn = .TRUE.

    ! Note: the snapshot is taken from inside CheckConvergence only when the
    ! full feasibility battery has been satisfied.  We do NOT snapshot at the
    ! top of PostConvergenceCheck because the prior round may have force-added
    ! a phase that is not yet mass-balanced.
    if (lTraceOn) then
        write(*,'(A)') '=== PostConvergenceCheck enter ==='
        write(*,'(A,I0,A,I0,A,I0)') '  nConPhases=', nConPhases, ' nSolnPhases=', nSolnPhases, &
            ' nSolnPhasesSys=', nSolnPhasesSys
        write(*,'(A)') '  Current assemblage:'
        do i = 1, nConPhases
            write(*,'(A,I4,A,I6,A,ES12.4)') '    slot ', i, ' species=', iAssemblage(i), &
                ' moles=', dMolesPhase(i)
        end do
        do i = nElements - nSolnPhases + 1, nElements
            write(*,'(A,I4,A,I6,A,ES12.4)') '    slot ', i, ' solnIdx=', -iAssemblage(i), &
                ' moles=', dMolesPhase(i)
            if (-iAssemblage(i) > 0) write(*,'(A,A)') '      name=', TRIM(cSolnPhaseName(-iAssemblage(i)))
        end do
    end if

    ! =====================================================================
    ! Part 1: Check all non-stable PURE CONDENSED phases
    ! =====================================================================
    call CompDrivingForce(iMaxDrivingForce, dMaxDrivingForce)
    if (dMaxDrivingForce < dTolerance(4)) then
        ! A pure condensed phase has negative driving force -- should be added.
        ! Force-add it directly.  If the assemblage is already full, remove the
        ! weakest pure condensed phase first.
        if (nConPhases + nSolnPhases >= nElements) then
            if (nConPhases <= 0) then
                lConv = .FALSE.
                return
            end if

            iWeakestCon = 1
            dWeakestConMoles = dMolesPhase(1)
            do i = 2, nConPhases
                if (dMolesPhase(i) < dWeakestConMoles) then
                    dWeakestConMoles = dMolesPhase(i)
                    iWeakestCon = i
                end if
            end do

            do i = iWeakestCon, nConPhases - 1
                iAssemblage(i) = iAssemblage(i + 1)
                dMolesPhase(i) = dMolesPhase(i + 1)
            end do
            dMolesPhase(nConPhases) = 0D0
            iAssemblage(nConPhases) = 0
            nConPhases = nConPhases - 1
        end if

        nConPhases = nConPhases + 1
        iAssemblage(nConPhases) = iMaxDrivingForce
        dMolesPhase(nConPhases) = dTolerance(9)  ! Min moles for new phase
        iterLast = iterGlobal
        lConv = .FALSE.
        return
    end if

    ! =====================================================================
    ! Part 2: Check all non-stable SOLUTION phases
    ! =====================================================================
    dBestDrivingForce = 0D0
    iBestSoln = 0
    iRetrySolnPhase = 0
    if (allocated(dRetryMolFraction)) dRetryMolFraction = 0D0
    if (allocated(dBestMolFraction)) deallocate(dBestMolFraction)
    allocate(dBestMolFraction(nSpeciesPhase(nSolnPhasesSys)))
    dBestMolFraction = 0D0

    do i = 1, nSolnPhasesSys
        if (lSolnPhases(i)) cycle

        ! Compute driving force via CompMolFraction (which calls Subminimization
        ! for non-ideal phases)
        if (cSolnPhaseType(i) == 'IDMX') then
            call CompMolFraction(i)
            if (lTraceOn) write(*,'(A,I3,A,A,A,ES14.6)') '  [DF] i=', i, ' name=', &
                TRIM(cSolnPhaseName(i)), ' IDMX df=', dDrivingForceSoln(i)
        else
            ! First try the phase's current composition estimate.  This often
            ! sits closer to a liquid basin than the extremum scans alone.
            call CompMolFraction(i)
            if (lTraceOn) write(*,'(A,I3,A,A,A,ES14.6)') '  [DF] i=', i, ' name=', &
                TRIM(cSolnPhaseName(i)), ' curstart df=', dDrivingForceSoln(i)
            if (dDrivingForceSoln(i) < dBestDrivingForce) then
                dBestDrivingForce = dDrivingForceSoln(i)
                iBestSoln = i
                iFirst = nSpeciesPhase(i-1) + 1
                iLast  = nSpeciesPhase(i)
                dBestMolFraction(iFirst:iLast) = dMolFraction(iFirst:iLast)
            end if

            call CheckMiscibilityGap(i, lPhaseChange)
            if (lTraceOn) write(*,'(A,I3,A,A,A,ES14.6,A,L1)') '  [DF] i=', i, ' name=', &
                TRIM(cSolnPhaseName(i)), ' MGmultistart df=', dDrivingForceSoln(i), ' changed=', lPhaseChange
            if (lPhaseChange) then
                if (dDrivingForceSoln(i) < dBestDrivingForce) then
                    dBestDrivingForce = dDrivingForceSoln(i)
                    iBestSoln = i
                    iFirst = nSpeciesPhase(i-1) + 1
                    iLast  = nSpeciesPhase(i)
                    dBestMolFraction(iFirst:iLast) = dMolFraction(iFirst:iLast)
                end if
            end if
        end if

        if (dDrivingForceSoln(i) < dBestDrivingForce) then
            dBestDrivingForce = dDrivingForceSoln(i)
            iBestSoln = i
            iFirst = nSpeciesPhase(i-1) + 1
            iLast  = nSpeciesPhase(i)
            dBestMolFraction(iFirst:iLast) = dMolFraction(iFirst:iLast)
        end if
    end do

    if (iBestSoln > 0 .AND. dBestDrivingForce < dTolerance(4)) then
        if (lTraceOn) write(*,'(A,I0,A,A,A,ES14.6)') '  [DECISION] force-add soln idx=', iBestSoln, &
            ' name=', TRIM(cSolnPhaseName(iBestSoln)), ' df=', dBestDrivingForce
        iRetrySolnPhase = iBestSoln
        if (.NOT. allocated(dRetryMolFraction)) allocate(dRetryMolFraction(nSpeciesPhase(nSolnPhasesSys)))
        dRetryMolFraction = 0D0
        dRetryMolFraction = dBestMolFraction
        ! Force-add the solution phase with most negative driving force
        if (nConPhases + nSolnPhases >= nElements) then
            if (nConPhases <= 0) then
                if (allocated(dBestMolFraction)) deallocate(dBestMolFraction)
                return
            end if

            ! If the assemblage is already full, evict the weakest pure condensed
            ! phase to make room for the missing solution phase.
            iWeakestCon = 1
            dWeakestConMoles = dMolesPhase(1)
            do i = 2, nConPhases
                if (dMolesPhase(i) < dWeakestConMoles) then
                    dWeakestConMoles = dMolesPhase(i)
                    iWeakestCon = i
                end if
            end do

            do i = iWeakestCon, nConPhases - 1
                iAssemblage(i) = iAssemblage(i + 1)
                dMolesPhase(i) = dMolesPhase(i + 1)
            end do
            dMolesPhase(nConPhases) = 0D0
            iAssemblage(nConPhases) = 0
            nConPhases = nConPhases - 1
        end if

        nSolnPhases = nSolnPhases + 1
        k = nElements - nSolnPhases + 1
        iAssemblage(k) = -iBestSoln

        ! Set initial moles and species moles
        dMolesPhase(k) = dTolerance(9)
        iFirst = nSpeciesPhase(iBestSoln-1) + 1
        iLast  = nSpeciesPhase(iBestSoln)
        dMolFraction(iFirst:iLast) = dBestMolFraction(iFirst:iLast)
        do j = iFirst, iLast
            dMolesSpecies(j) = dMolFraction(j) * dMolesPhase(k)
        end do

        lSolnPhases(iBestSoln) = .TRUE.
        dDrivingForceSoln(iBestSoln) = 0D0
        iterLast = iterGlobal

        if (allocated(dBestMolFraction)) deallocate(dBestMolFraction)
        lConv = .FALSE.
        return
    end if

    if (allocated(dBestMolFraction)) deallocate(dBestMolFraction)

end subroutine PostConvergenceCheck

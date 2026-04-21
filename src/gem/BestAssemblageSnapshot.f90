!-------------------------------------------------------------------------------------------------------------
!> \file    BestAssemblageSnapshot.f90
!> \brief   Best-Gibbs assemblage snapshot support: records the lowest-Gibbs
!!          feasible (converged) assemblage seen at any point during the GEM
!!          solve.  If the final converged state is a local minimum that is
!!          higher in Gibbs energy than an earlier transit, the snapshot can
!!          be restored to reach the true global minimum.
!-------------------------------------------------------------------------------------------------------------

subroutine UpdateBestAssemblageSnapshot

    !-------------------------------------------------------------------------
    ! If the current converged assemblage has a lower integral Gibbs energy
    ! than any previously stored snapshot, copy the full state into the
    ! best-snapshot buffers.
    !-------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: i
    real(8) :: dGibbsSys
    character(len=32) :: cTraceEnv
    logical :: lTraceOn

    lTraceOn = .FALSE.
    call GET_ENVIRONMENT_VARIABLE('TC_TRACE_POSTCONV', cTraceEnv)
    if (TRIM(cTraceEnv) == '1') lTraceOn = .TRUE.

    ! Compute integral Gibbs energy of the system from element potentials and
    ! element totals.  This matches the convention used by CheckConvergence's
    ! running minimum (dMinGibbs).
    dGibbsSys = 0D0
    do i = 1, nElements
        dGibbsSys = dGibbsSys + dElementPotential(i) * dMolesElement(i) &
                              * dTemperature * dIdealConstant
    end do

    ! Reject if assemblage contains dummy species or fewer than 1 phase
    if (nConPhases + nSolnPhases <= 0) return

    if (.NOT. lBestAssemblageValid .OR. dGibbsSys < dBestSnapshotGibbs) then
        dBestSnapshotGibbs        = dGibbsSys
        nBestConPhases            = nConPhases
        nBestSolnPhases           = nSolnPhases
        iBestAssemblage           = iAssemblage
        dBestMolesPhaseSnap       = dMolesPhase
        dBestElementPotentialSnap = dElementPotential
        dBestMolFractionSnap      = dMolFraction
        dBestMolesSpeciesSnap     = dMolesSpecies
        lBestSolnPhasesSnap       = lSolnPhases
        lBestAssemblageValid      = .TRUE.
        if (lTraceOn) then
            write(*,'(A,ES16.8,A,I0,A,I0)') '[SNAPSHOT] saved G=', dGibbsSys, &
                ' nCon=', nConPhases, ' nSoln=', nSolnPhases
            do i = 1, nConPhases
                write(*,'(A,I0,A,I0,A,ES12.4)') '  con slot ', i, ' sp=', iAssemblage(i), ' mol=', dMolesPhase(i)
            end do
            do i = nElements - nSolnPhases + 1, nElements
                write(*,'(A,I0,A,I0,A,ES12.4)') '  soln slot ', i, ' idx=', -iAssemblage(i), ' mol=', dMolesPhase(i)
            end do
        end if
    else
        if (lTraceOn) write(*,'(A,ES16.8,A,ES16.8)') '[SNAPSHOT] skip: cur=', dGibbsSys, ' best=', dBestSnapshotGibbs
    end if

end subroutine UpdateBestAssemblageSnapshot


subroutine RestoreBestAssemblageSnapshot

    !-------------------------------------------------------------------------
    ! Restore the state saved by UpdateBestAssemblageSnapshot into the
    ! primary solver arrays.  Used when the final converged state has a
    ! higher Gibbs energy than an earlier transit.
    !-------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    if (.NOT. lBestAssemblageValid) return

    nConPhases        = nBestConPhases
    nSolnPhases       = nBestSolnPhases
    iAssemblage       = iBestAssemblage
    dMolesPhase       = dBestMolesPhaseSnap
    dElementPotential = dBestElementPotentialSnap
    dMolFraction      = dBestMolFractionSnap
    dMolesSpecies     = dBestMolesSpeciesSnap
    lSolnPhases       = lBestSolnPhasesSnap

end subroutine RestoreBestAssemblageSnapshot

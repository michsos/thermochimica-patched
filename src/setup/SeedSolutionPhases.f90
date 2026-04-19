
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    SeedSolutionPhases.f90
    !> \brief   Seed the Leveling assemblage with solution phases that have negative driving forces.
    !
    !> \details After the LevelingSolver produces an initial assemblage (typically all pure condensed
    !!          phases), this routine checks each solution phase for stability using the Subminimization
    !!          algorithm, which accounts for non-ideal mixing energies.  If a solution phase has a
    !!          negative driving force (indicating it is more stable than the current Gibbs plane
    !!          predicts), it is swapped into the assemblage in place of the weakest pure condensed phase.
    !!
    !!          This is critical for SUBQ phases (like Slag-liq) whose stability derives primarily from
    !!          short-range-order mixing contributions that the linear Leveling algorithm cannot capture.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine SeedSolutionPhases

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, iFirst, iLast, nVar, iBestPhase, iWeakestSlot
    real(8) :: dBestDrivingForce, dWeakestMoles, dDF
    logical :: lPhasePass

    ! Only proceed if the leveling assemblage is all pure condensed phases
    ! and there are solution phases in the system to seed.
    if (nSolnPhasesSys == 0) return
    if (INFOThermo /= 0) return

    ! Check if there are any solution phases already in the assemblage
    ! (from PostLevelingSolver).  If so, skip.
    do i = 1, nElements
        if (iPhase(iAssemblage(i)) > 0) return  ! Solution phase already present
    end do

    ! We need dMolFraction, dMolesSpecies, etc. for Subminimization.
    ! These are allocated in InitGEMSolver, but we're before that.
    ! Allocate temporary arrays if needed.
    if (.NOT. allocated(dMolFraction)) allocate(dMolFraction(nSpecies))
    if (.NOT. allocated(dMolesSpecies)) allocate(dMolesSpecies(nSpecies))
    if (.NOT. allocated(dPartialExcessGibbs)) allocate(dPartialExcessGibbs(nSpecies))
    if (.NOT. allocated(dPartialExcessGibbsLast)) allocate(dPartialExcessGibbsLast(nSpecies))
    if (.NOT. allocated(dSumMolFractionSoln)) then
        k = MAX(1, nSolnPhasesSys)
        allocate(dSumMolFractionSoln(k))
    end if
    if (.NOT. allocated(dDrivingForceSoln)) then
        k = MAX(1, nSolnPhasesSys)
        allocate(dDrivingForceSoln(k))
    end if
    if (.NOT. allocated(dEffStoichSolnPhase)) then
        k = MAX(1, nSolnPhasesSys)
        allocate(dEffStoichSolnPhase(k, nElements))
    end if
    if (.NOT. allocated(lSolnPhases)) then
        k = MAX(1, nSolnPhasesSys)
        allocate(lSolnPhases(k))
    end if
    if (.NOT. allocated(lMiscibility)) then
        k = MAX(1, nSolnPhasesSys)
        allocate(lMiscibility(k))
    end if

    ! Initialize arrays
    dMolFraction = 0D0
    dMolesSpecies = 0D0
    dPartialExcessGibbs = 0D0
    dPartialExcessGibbsLast = 0D0
    dSumMolFractionSoln = 0D0
    dDrivingForceSoln = 0D0
    dEffStoichSolnPhase = 0D0
    lSolnPhases = .FALSE.
    lMiscibility = .FALSE.

    ! Initialize mole fractions for all solution phases
    ! (uniform distribution as initial estimate)
    do k = 1, nSolnPhasesSys
        nVar = nSpeciesPhase(k) - nSpeciesPhase(k - 1)
        if (nVar > 0) then
            do i = nSpeciesPhase(k - 1) + 1, nSpeciesPhase(k)
                dMolFraction(i) = 1D0 / DBLE(nVar)
            end do
        end if
    end do

    ! Find the solution phase with the most negative driving force
    iBestPhase = 0
    dBestDrivingForce = 0D0

    LOOP_SolnPhases: do i = 1, nSolnPhasesSys
        ! Skip if this phase has too few species
        nVar = nSpeciesPhase(i) - nSpeciesPhase(i - 1)
        if (nVar < 1) cycle

        ! Use CheckMiscibilityGap which does multi-start Subminimization
        ! from different composition extremes
        call CheckMiscibilityGap(i, lPhasePass)

        ! Get the driving force
        dDF = dDrivingForceSoln(i)

        if (lPhasePass .AND. (dDF < dBestDrivingForce)) then
            dBestDrivingForce = dDF
            iBestPhase = i
        end if
    end do LOOP_SolnPhases

    ! If we found a solution phase with negative driving force, swap it in
    if (iBestPhase > 0) then
        ! Find the weakest (smallest moles) pure condensed phase to replace
        iWeakestSlot = 1
        dWeakestMoles = dMolesPhase(1)
        do j = 2, nElements
            if (dMolesPhase(j) < dWeakestMoles) then
                dWeakestMoles = dMolesPhase(j)
                iWeakestSlot = j
            end if
        end do

        ! Swap: remove the weakest pure condensed phase, add the solution phase
        ! Shift the assemblage: move the last pure condensed into the gap
        if (iWeakestSlot < nElements) then
            iAssemblage(iWeakestSlot) = iAssemblage(nElements)
            dMolesPhase(iWeakestSlot) = dMolesPhase(nElements)
        end if

        ! Put the solution phase in the last slot (standard convention)
        iAssemblage(nElements) = -iBestPhase
        dMolesPhase(nElements) = dWeakestMoles  ! Give it the removed phase's moles

        ! Update phase counts
        nConPhases = nElements - 1
    end if

    return

end subroutine SeedSolutionPhases

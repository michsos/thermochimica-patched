
subroutine RetryCalculationFirstPhase
    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleReinit
    USE ModuleGEMSolver

    implicit none

    logical :: oldReinit, lHaveRetryMolFraction, lStableRetryFamily, lStableSlag
    integer :: iPhaseRetry, iFirst, iLast, iFallbackPhase, iStableSlot, iStablePhase

    oldReinit = lReinitRequested
    if ((INFOThermo == 12) .AND. (nSolnPhasesSys > 0)) then
        iPhaseRetry = iRetrySolnPhase
        if ((iPhaseRetry <= 0) .OR. (iPhaseRetry > nSolnPhasesSys)) iPhaseRetry = 1

        ! Avoid falling back to gas_ideal as the retry seed. Gas is useful for
        ! thermodynamic bookkeeping, but a gas-only reinit state is a poor
        ! escape path for condensed local minima.
        if ((cSolnPhaseType(iPhaseRetry) == 'IDMX') .AND. (TRIM(cSolnPhaseName(iPhaseRetry)) == 'gas_ideal') .AND. &
            (nSolnPhases > 0)) then
            iFallbackPhase = 0
            do iFallbackPhase = 1, nSolnPhasesSys
                if (iFallbackPhase == iPhaseRetry) cycle
                if ((cSolnPhaseType(iFallbackPhase) == 'IDMX') .AND. &
                    (TRIM(cSolnPhaseName(iFallbackPhase)) == 'gas_ideal')) cycle
                if (ANY(iAssemblage == -iFallbackPhase)) cycle
                iPhaseRetry = iFallbackPhase
                exit
            end do
        end if

        ! If the retry phase family is already stable and the current
        ! assemblage has no liquid, prefer re-seeding with a non-stable liquid
        ! variant. This helps escape alkali-silicate local minima that keep
        ! retrying the same solid solution family instead of exploring the
        ! missing Slag-liq basin.
        lStableRetryFamily = .FALSE.
        lStableSlag = .FALSE.
        do iStableSlot = nElements - nSolnPhases + 1, nElements
            if (iStableSlot < 1 .OR. iStableSlot > nElements) cycle
            iStablePhase = -iAssemblage(iStableSlot)
            if (iStablePhase <= 0) cycle
            if (TRIM(cSolnPhaseName(iStablePhase)) == TRIM(cSolnPhaseName(iPhaseRetry))) lStableRetryFamily = .TRUE.
            if (TRIM(cSolnPhaseName(iStablePhase)) == 'Slag-liq') lStableSlag = .TRUE.
        end do

        if (lStableRetryFamily .AND. (.NOT. lStableSlag) .AND. (TRIM(cSolnPhaseName(iPhaseRetry)) /= 'Slag-liq')) then
            do iFallbackPhase = 1, nSolnPhasesSys
                if (TRIM(cSolnPhaseName(iFallbackPhase)) /= 'Slag-liq') cycle
                if (ANY(iAssemblage == -iFallbackPhase)) cycle
                iPhaseRetry = iFallbackPhase
                exit
            end do
        end if

        iFirst = nSpeciesPhase(iPhaseRetry - 1) + 1
        iLast  = nSpeciesPhase(iPhaseRetry)
        lHaveRetryMolFraction = .FALSE.
        if (allocated(dRetryMolFraction)) then
            if ((iFirst >= 1) .AND. (iLast <= SIZE(dRetryMolFraction))) then
                lHaveRetryMolFraction = SUM(dRetryMolFraction(iFirst:iLast)) > 0D0
            end if
        end if

        if (iLast >= iFirst) then
            INFOThermo = 0
            lReinitRequested = .TRUE.
            call PostProcessThermo
            nSolnPhases = 1
            call SaveReinitData
            if (lReinitAvailable) then
                ! If post-convergence already identified a better assemblage that
                ! includes the retry phase, preserve that full assemblage for the
                ! recursive solve instead of collapsing back to a single uniform
                ! solution-phase seed.
                if (ANY(iAssemblage == -iPhaseRetry)) then
                    iAssemblage_Old = iAssemblage
                    dMolesPhase_Old = dMolesPhase
                    dMolFraction_Old = dMolFraction
                    if (lHaveRetryMolFraction) then
                        dMolFraction_Old(iFirst:iLast) = dRetryMolFraction(iFirst:iLast)
                    end if
                else
                    iAssemblage_Old = 0
                    iAssemblage_Old(nElements) = -iPhaseRetry
                    dMolFraction_Old = 0D0
                    if (lHaveRetryMolFraction) then
                        dMolFraction_Old(iFirst:iLast) = dRetryMolFraction(iFirst:iLast)
                    else if (SUM(dMolFraction(iFirst:iLast)) > 0D0) then
                        dMolFraction_Old(iFirst:iLast) = dMolFraction(iFirst:iLast)
                    else
                        dMolFraction_Old(iFirst:iLast) = 1D0 / DFLOAT(iLast - iFirst + 1)
                    end if
                    dMolesPhase_Old = 0D0
                    dMolesPhase_Old(nElements) = 1D0
                    ! When collapsing to a single-phase retry seed, discard the
                    ! previous local-minimum element-potential state. Reusing it
                    ! can pin the recursive solve to the same basin even with a
                    ! different solution-phase seed.
                    dElementPotential_Old = 0D0
                end if
                dChemicalPotential_Old = 0D0
                call ResetThermo
                lRetryAttempted = .TRUE.
                call Thermochimica
                iterGlobal = iterGlobal + iterGlobalMax
                ! Reset flag to original
                if (INFOThermo /= 0) INFOThermo = 12
            else
                INFOThermo = 12
                dElementMass = dElementMass * dNormalizeInput / dMassScale
            end if
            lReinitRequested = oldReinit
        end if
    end if

end subroutine RetryCalculationFirstPhase

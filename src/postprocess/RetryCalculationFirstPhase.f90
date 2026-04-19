
subroutine RetryCalculationFirstPhase
    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleReinit
    USE ModuleGEMSolver

    implicit none

    logical :: oldReinit, lHaveRetryMolFraction
    integer :: iPhaseRetry, iFirst, iLast

    oldReinit = lReinitRequested
    if ((INFOThermo == 12) .AND. (nSolnPhasesSys > 0)) then
        iPhaseRetry = iRetrySolnPhase
        if ((iPhaseRetry <= 0) .OR. (iPhaseRetry > nSolnPhasesSys)) iPhaseRetry = 1

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

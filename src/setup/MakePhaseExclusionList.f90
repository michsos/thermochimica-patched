subroutine MakePhaseExclusionList

    USE ModuleThermoIO
    USE ModuleParseCS

    implicit none

    integer :: i, j
    logical :: lFound

    ! If there is an "excluded except" list, then add all other phases to exclusion list
    if (nPhasesExcludedExcept > 0) then
        ! Solution phases
        loop_checkExclusionSolution: do i = 1, nSolnPhasesSysCS
            ! Check if phase is on exception list (exact or prefix match)
            do j = 1, nPhasesExcludedExcept
                if (cSolnPhaseNameCS(i) == cPhasesExcludedExcept(j)) cycle loop_checkExclusionSolution
                if (LEN_TRIM(cPhasesExcludedExcept(j)) > 0) then
                    if (INDEX(TRIM(cSolnPhaseNameCS(i)), TRIM(cPhasesExcludedExcept(j))) == 1) &
                        cycle loop_checkExclusionSolution
                end if
            end do
            ! Check if phase is already on exclusion list
            do j = 1, nPhasesExcluded
                if (cSolnPhaseNameCS(i) == cPhasesExcluded(j)) cycle loop_checkExclusionSolution
            end do
            ! If not, add to exclusion list
            nPhasesExcluded = nPhasesExcluded + 1
            cPhasesExcluded(nPhasesExcluded) = cSolnPhaseNameCS(i)
        end do loop_checkExclusionSolution

        ! Pure condensed phases
        loop_checkExclusionPureCondensed: do i = nSpeciesPhaseCS(nSolnPhasesSysCS) + 1, nSpeciesCS
            ! Check if dummy - don't exclude dummies automatically
            if (iPhaseCS(i) == -1) cycle loop_checkExclusionPureCondensed
            ! Check if phase is on exception list (exact or prefix match)
            do j = 1, nPhasesExcludedExcept
                if (cSpeciesNameCS(i) == cPhasesExcludedExcept(j)) cycle loop_checkExclusionPureCondensed
                if (LEN_TRIM(cPhasesExcludedExcept(j)) > 0) then
                    if (INDEX(TRIM(cSpeciesNameCS(i)), TRIM(cPhasesExcludedExcept(j))) == 1) &
                        cycle loop_checkExclusionPureCondensed
                end if
            end do
            ! Check if phase is already on exclusion list
            do j = 1, nPhasesExcluded
                if (cSpeciesNameCS(i) == cPhasesExcluded(j)) cycle loop_checkExclusionPureCondensed
            end do
            ! If not, add to exclusion list
            nPhasesExcluded = nPhasesExcluded + 1
            cPhasesExcluded(nPhasesExcluded) = cSpeciesNameCS(i)
        end do loop_checkExclusionPureCondensed
    end if


end subroutine MakePhaseExclusionList

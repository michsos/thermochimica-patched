subroutine ResetDormantPhaseCandidateSummary

    USE ModuleThermoIO, ONLY: iStrongestDormantPureConPhaseIndex, dStrongestDormantPureConDrivingForce, &
        iStrongestDormantSolnPhaseIndex, dStrongestDormantSolnDrivingForce, iDormantPhaseAddOrder

    implicit none

    iStrongestDormantPureConPhaseIndex = 0
    dStrongestDormantPureConDrivingForce = 0D0
    iStrongestDormantSolnPhaseIndex = 0
    dStrongestDormantSolnDrivingForce = 0D0
    iDormantPhaseAddOrder = 0

    return
end subroutine ResetDormantPhaseCandidateSummary

logical function IsSolutionPhaseDormant(iSolnPhaseIndex)

    USE ModuleThermoIO, ONLY: lPhaseDormancyActive, lDormantSolnPhases

    implicit none

    integer, intent(in) :: iSolnPhaseIndex

    IsSolutionPhaseDormant = .FALSE.

    if (.NOT. lPhaseDormancyActive) return
    if (.NOT. allocated(lDormantSolnPhases)) return
    if ((iSolnPhaseIndex < 1) .OR. (iSolnPhaseIndex > size(lDormantSolnPhases))) return

    IsSolutionPhaseDormant = lDormantSolnPhases(iSolnPhaseIndex)

    return
end function IsSolutionPhaseDormant

logical function IsPureConSpeciesDormant(iSpeciesIndex)

    USE ModuleThermoIO, ONLY: lPhaseDormancyActive, lDormantSpecies

    implicit none

    integer, intent(in) :: iSpeciesIndex

    IsPureConSpeciesDormant = .FALSE.

    if (.NOT. lPhaseDormancyActive) return
    if (.NOT. allocated(lDormantSpecies)) return
    if ((iSpeciesIndex < 1) .OR. (iSpeciesIndex > size(lDormantSpecies))) return

    IsPureConSpeciesDormant = lDormantSpecies(iSpeciesIndex)

    return
end function IsPureConSpeciesDormant

logical function IsSpeciesDormant(iSpeciesIndex)

    USE ModuleThermo, ONLY: iPhase

    implicit none

    integer, intent(in) :: iSpeciesIndex
    logical :: IsPureConSpeciesDormant, IsSolutionPhaseDormant

    IsSpeciesDormant = .FALSE.

    if (iSpeciesIndex <= 0) return

    if (iPhase(iSpeciesIndex) == 0) then
        IsSpeciesDormant = IsPureConSpeciesDormant(iSpeciesIndex)
    else if (iPhase(iSpeciesIndex) > 0) then
        IsSpeciesDormant = IsSolutionPhaseDormant(iPhase(iSpeciesIndex))
    end if

    return
end function IsSpeciesDormant

logical function IsAssemblagePhaseDormant(iAssemblageIndex)

    USE ModuleThermo

    implicit none

    integer, intent(in) :: iAssemblageIndex
    logical :: IsPureConSpeciesDormant, IsSolutionPhaseDormant

    IsAssemblagePhaseDormant = .FALSE.

    if (iAssemblageIndex > 0) then
        IsAssemblagePhaseDormant = IsPureConSpeciesDormant(iAssemblageIndex)
    else if (iAssemblageIndex < 0) then
        IsAssemblagePhaseDormant = IsSolutionPhaseDormant(-iAssemblageIndex)
    end if

    return
end function IsAssemblagePhaseDormant

subroutine RecordDormantPhaseCandidateSummary(iPureSpeciesIndex, dPureDrivingForce, iSolnPhaseIndex, dSolnDrivingForce)

    USE ModuleThermo, ONLY: nSolnPhasesSys
    USE ModuleThermoIO, ONLY: iStrongestDormantPureConPhaseIndex, dStrongestDormantPureConDrivingForce, &
        iStrongestDormantSolnPhaseIndex, dStrongestDormantSolnDrivingForce, iDormantPhaseAddOrder

    implicit none

    integer, intent(in) :: iPureSpeciesIndex, iSolnPhaseIndex
    real(8), intent(in) :: dPureDrivingForce, dSolnDrivingForce

    iStrongestDormantPureConPhaseIndex = 0
    dStrongestDormantPureConDrivingForce = 0D0
    if (iPureSpeciesIndex > 0) then
        iStrongestDormantPureConPhaseIndex = nSolnPhasesSys + iPureSpeciesIndex
        dStrongestDormantPureConDrivingForce = dPureDrivingForce
    end if

    iStrongestDormantSolnPhaseIndex = 0
    dStrongestDormantSolnDrivingForce = 0D0
    if (iSolnPhaseIndex > 0) then
        iStrongestDormantSolnPhaseIndex = iSolnPhaseIndex
        dStrongestDormantSolnDrivingForce = dSolnDrivingForce
    end if

    iDormantPhaseAddOrder = 0
    if ((iPureSpeciesIndex > 0) .AND. (iSolnPhaseIndex > 0)) then
        if (dPureDrivingForce < dSolnDrivingForce) then
            iDormantPhaseAddOrder = 1
        else
            iDormantPhaseAddOrder = 2
        end if
    else if (iPureSpeciesIndex > 0) then
        iDormantPhaseAddOrder = 1
    else if (iSolnPhaseIndex > 0) then
        iDormantPhaseAddOrder = 2
    end if

    return
end subroutine RecordDormantPhaseCandidateSummary

subroutine ResolveDormantPhaseMasks

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer :: i, j, species_start, species_end
    logical :: lMatch

    call ResetDormantPhaseCandidateSummary

    if (allocated(lDormantSolnPhases)) deallocate(lDormantSolnPhases)
    if (allocated(lDormantSpecies)) deallocate(lDormantSpecies)

    if (nSolnPhasesSys > 0) allocate(lDormantSolnPhases(nSolnPhasesSys))
    if (nSpecies > 0) allocate(lDormantSpecies(nSpecies))

    if (allocated(lDormantSolnPhases)) lDormantSolnPhases = .FALSE.
    if (allocated(lDormantSpecies)) lDormantSpecies = .FALSE.

    if (.NOT. lPhaseDormancyActive) return

    if ((nPhasesDormant <= 0) .AND. (nPhasesDormantExcept <= 0)) return

    if (allocated(lDormantSolnPhases) .AND. (nPhasesDormantExcept > 0)) lDormantSolnPhases = .TRUE.
    if (allocated(lDormantSpecies) .AND. (nPhasesDormantExcept > 0)) then
        do i = nSpeciesPhase(nSolnPhasesSys) + 1, nSpecies
            if (iPhase(i) == 0) lDormantSpecies(i) = .TRUE.
        end do
    end if

    do i = 1, nSolnPhasesSys
        lMatch = .FALSE.
        do j = 1, nPhasesDormant
            if (trim(adjustl(cSolnPhaseName(i))) == trim(adjustl(cPhasesDormant(j)))) then
                lMatch = .TRUE.
                exit
            end if
        end do
        if (lMatch .AND. allocated(lDormantSolnPhases)) lDormantSolnPhases(i) = .TRUE.

        lMatch = .FALSE.
        do j = 1, nPhasesDormantExcept
            if (trim(adjustl(cSolnPhaseName(i))) == trim(adjustl(cPhasesDormantExcept(j)))) then
                lMatch = .TRUE.
                exit
            end if
        end do
        if (lMatch .AND. allocated(lDormantSolnPhases)) lDormantSolnPhases(i) = .FALSE.
    end do

    species_start = nSpeciesPhase(nSolnPhasesSys) + 1
    do i = species_start, nSpecies
        if (iPhase(i) /= 0) cycle

        lMatch = .FALSE.
        do j = 1, nPhasesDormant
            if (trim(adjustl(cSpeciesName(i))) == trim(adjustl(cPhasesDormant(j)))) then
                lMatch = .TRUE.
                exit
            end if
        end do
        if (lMatch .AND. allocated(lDormantSpecies)) lDormantSpecies(i) = .TRUE.

        lMatch = .FALSE.
        do j = 1, nPhasesDormantExcept
            if (trim(adjustl(cSpeciesName(i))) == trim(adjustl(cPhasesDormantExcept(j)))) then
                lMatch = .TRUE.
                exit
            end if
        end do
        if (lMatch .AND. allocated(lDormantSpecies)) lDormantSpecies(i) = .FALSE.
    end do

    do i = 1, nSolnPhasesSys
        if (.NOT. allocated(lDormantSolnPhases)) exit
        if (.NOT. lDormantSolnPhases(i)) cycle
        species_start = nSpeciesPhase(i - 1) + 1
        species_end = nSpeciesPhase(i)
        if (allocated(lDormantSpecies)) lDormantSpecies(species_start:species_end) = .TRUE.
    end do

    return
end subroutine ResolveDormantPhaseMasks

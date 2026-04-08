
!-------------------------------------------------------------------------------
!
!> \file    GetPureConPhaseDrivingForce.f90
!> \brief   Get direct driving-force output for a pure condensed phase.
!
!-------------------------------------------------------------------------------


subroutine GetPureConPhaseDrivingForce(cPureConOut, dDrivingForceOut, INFO)

    USE ModuleThermo, ONLY: nElements, nSpecies, nDummySpecies, nSolnPhasesSys, nSpeciesPhase, &
        cSpeciesName, dStdGibbsEnergy, dStoichSpecies, dElementPotential, dSpeciesTotalAtoms
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer,       intent(out)   :: INFO
    integer                      :: i, j
    real(8),       intent(out)   :: dDrivingForceOut
    real(8)                      :: dElementPotentialPhase
    character(*),  intent(in)    :: cPureConOut
    character(30)                :: cPhaseTemp


    INFO             = 0
    dDrivingForceOut = 0D0

    if ((INFOThermo /= 0) .AND. (INFOThermo /= 12)) then
        INFO = -1
        return
    end if

    cPhaseTemp = cPureConOut(1:min(30, len(cPureConOut)))
    cPhaseTemp = trim(adjustl(cPhaseTemp))

    i = 0
    do j = nSpeciesPhase(nSolnPhasesSys) + 1, nSpecies - nDummySpecies
        if (cPhaseTemp == trim(adjustl(cSpeciesName(j)))) then
            i = j
            exit
        end if
    end do

    if (i == 0) then
        INFO = 1
        return
    end if

    dElementPotentialPhase = 0D0
    do j = 1, nElements
        dElementPotentialPhase = dElementPotentialPhase + dElementPotential(j) * dStoichSpecies(i, j)
    end do

    if (dSpeciesTotalAtoms(i) <= 0D0) then
        INFO = 1
        return
    end if

    dDrivingForceOut = dStdGibbsEnergy(i) - dElementPotentialPhase
    dDrivingForceOut = dDrivingForceOut / dSpeciesTotalAtoms(i)

    return

end subroutine GetPureConPhaseDrivingForce

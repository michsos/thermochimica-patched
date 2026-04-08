
!-------------------------------------------------------------------------------
!
!> \file    GetPureConPhaseDrivingForce.f90
!> \brief   Get direct driving-force output for a pure condensed phase.
!
!-------------------------------------------------------------------------------


subroutine GetPureConPhaseDrivingForce(cPureConOut, dDrivingForceOut, INFO)

    USE ModuleParseCS, ONLY: nSpeciesPhaseCS, nSolnPhasesSysCS, nSpeciesCS, cSpeciesNameCS, &
        nGibbsEqSpecies, dGibbsCoeffSpeciesTemp, dGibbsMagneticCS, dStoichSpeciesCS, iParticlesPerMoleCS
    USE ModuleThermo, ONLY: nElemOrComp, iElementSystem, dElementPotential, dIdealConstant
    USE ModuleThermoIO, ONLY: INFOThermo, dTemperature

    implicit none

    integer,       intent(out)   :: INFO
    integer                      :: i, j, k, l1, l2, iCounterGibbsEqn
    real(8),       intent(out)   :: dDrivingForceOut
    real(8)                      :: dStdGibbsEnergyPhase, dElementPotentialPhase
    real(8)                      :: dSpeciesTotalAtomsPhase, dTemp, dLogT
    real(8)                      :: B, D, p, invpmone, tau, Tcritical, g
    real(8)                      :: structure_factor
    real(8),       dimension(6)  :: dGibbsCoeff
    character(*),  intent(in)    :: cPureConOut
    character(30)                :: cPhaseTemp, cSpeciesTemp


    INFO             = 0
    dDrivingForceOut = 0D0

    if ((INFOThermo /= 0) .AND. (INFOThermo /= 12)) then
        INFO = -1
        return
    end if

    cPhaseTemp = cPureConOut(1:min(30, len(cPureConOut)))
    cPhaseTemp = trim(adjustl(cPhaseTemp))

    i = 0
    do j = nSpeciesPhaseCS(nSolnPhasesSysCS) + 1, nSpeciesCS
        cSpeciesTemp = trim(adjustl(cSpeciesNameCS(j)))
        if (cPhaseTemp == cSpeciesTemp) then
            i = j
            exit
        end if
    end do

    if (i == 0) then
        INFO = 1
        return
    end if

    dTemp          = 1D0 / (dIdealConstant * dTemperature)
    dGibbsCoeff(1) = 1D0
    dGibbsCoeff(2) = dTemperature
    dGibbsCoeff(3) = dTemperature * DLOG(dTemperature)
    dGibbsCoeff(4) = dTemperature**2
    dGibbsCoeff(5) = dTemperature**3
    dGibbsCoeff(6) = 1D0 / dTemperature
    dLogT          = DLOG(dTemperature)

    iCounterGibbsEqn = 0
    do j = 1, i - 1
        iCounterGibbsEqn = iCounterGibbsEqn + nGibbsEqSpecies(j)
    end do

    l1 = 0
    do k = 1, nGibbsEqSpecies(i)
        if ((dTemperature <= dGibbsCoeffSpeciesTemp(1, iCounterGibbsEqn + k)) .AND. (l1 == 0)) then
            l1 = k
        end if
    end do

    l2 = l1
    if (l2 == 0) l2 = nGibbsEqSpecies(i)
    l2 = l2 + iCounterGibbsEqn

    dStdGibbsEnergyPhase = 0D0
    do k = 2, 7
        dStdGibbsEnergyPhase = dStdGibbsEnergyPhase + dGibbsCoeffSpeciesTemp(k, l2) * dGibbsCoeff(k - 1)
    end do

    if (dGibbsCoeffSpeciesTemp(9, l2) == 99D0) then
        dStdGibbsEnergyPhase = dStdGibbsEnergyPhase + dGibbsCoeffSpeciesTemp(8, l2) * dLogT
    else
        dStdGibbsEnergyPhase = dStdGibbsEnergyPhase + dGibbsCoeffSpeciesTemp(8, l2) * dTemperature**dGibbsCoeffSpeciesTemp(9, l2)
    end if

    if (dGibbsCoeffSpeciesTemp(11, l2) == 99D0) then
        dStdGibbsEnergyPhase = dStdGibbsEnergyPhase + dGibbsCoeffSpeciesTemp(10, l2) * dLogT
    else
        dStdGibbsEnergyPhase = dStdGibbsEnergyPhase + dGibbsCoeffSpeciesTemp(10, l2) * dTemperature**dGibbsCoeffSpeciesTemp(11, l2)
    end if

    if (dGibbsCoeffSpeciesTemp(13, l2) == 99D0) then
        dStdGibbsEnergyPhase = dStdGibbsEnergyPhase + dGibbsCoeffSpeciesTemp(12, l2) * dLogT
    else
        dStdGibbsEnergyPhase = dStdGibbsEnergyPhase + dGibbsCoeffSpeciesTemp(12, l2) * dTemperature**dGibbsCoeffSpeciesTemp(13, l2)
    end if

    if (dGibbsMagneticCS(i, 1) /= 0D0) then
        Tcritical       = dGibbsMagneticCS(i, 1)
        B               = dGibbsMagneticCS(i, 2)
        structure_factor = dGibbsMagneticCS(i, 3)
        p               = dGibbsMagneticCS(i, 4)
        invpmone        = 1D0 / p - 1D0

        if (Tcritical < 0D0) then
            Tcritical = -Tcritical * structure_factor
            B         = -B * structure_factor
        end if

        tau = dTemperature / Tcritical
        D   = 518D0 / 1125D0 + (11692D0 / 15975D0) * invpmone

        if (tau > 1D0) then
            g = -((tau**(-5)) / 10D0 + (tau**(-15)) / 315D0 + (tau**(-25)) / 1500D0) / D
        else
            g = 1D0 - (79D0 / (140D0 * p * tau) + (474D0 / 497D0) * invpmone * ((tau**3) / 6D0 + (tau**9) / 135D0 + (tau**15) / 600D0)) / D
        end if

        dStdGibbsEnergyPhase = dStdGibbsEnergyPhase + dIdealConstant * dTemperature * DLOG(B + 1D0) * g
    end if

    dStdGibbsEnergyPhase = dStdGibbsEnergyPhase * dTemp

    dElementPotentialPhase = 0D0
    dSpeciesTotalAtomsPhase = 0D0
    j = 0
    do k = 1, nElemOrComp
        if (iElementSystem(k) == 0) then
            if (DABS(dStoichSpeciesCS(i, k)) > 0D0) then
                INFO = 1
                return
            end if
            cycle
        end if

        j = j + 1
        dElementPotentialPhase = dElementPotentialPhase + dElementPotential(j) * dStoichSpeciesCS(i, k)
        dSpeciesTotalAtomsPhase = dSpeciesTotalAtomsPhase + DABS(dStoichSpeciesCS(i, k))
    end do

    if (dSpeciesTotalAtomsPhase <= 0D0) then
        INFO = 1
        return
    end if

    dDrivingForceOut = (dStdGibbsEnergyPhase - dElementPotentialPhase) / dSpeciesTotalAtomsPhase
    dDrivingForceOut = dDrivingForceOut * dIdealConstant * dTemperature

    return

end subroutine GetPureConPhaseDrivingForce

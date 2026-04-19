
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckMiscibilityGap.f90
    !> \brief   Check for a miscibility gap.
    !> \author  M.H.A. Piro
    !> \date    Aug. 30, 2012
    !> \sa      Subminimization.f90
    !> \sa      CheckSolnPhaseAdd.f90
    !> \sa      CheckConvergence.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   08/30/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check whether a particular non-ideal solution phase
    !! containing a miscibility gap should be added to the system.  The subminimization routine determines
    !! whether the driving force associated with a local minima is positive/negative.  There may be multiple
    !! local minima for a particular solution phase and different local minima may be discovered depending on
    !! the initial estimates of the mole fractions of this phase.  The approach taken here performs
    !! subminimization multiple times where the initial estimates for each case begin at the extremums of the
    !! domain space.  Specifically, the mole fractions of all constituents are set to an arbitrarily small
    !! value (e.g., 1D-3) except for one constituent where the sum of the mole fractions equals unity.  Each
    !! constituent in this phase is systematically initialized as being dominant.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iSolnPhaseIndex  An integer scalar representing the absolute index of the solution phase that
    !!                                is being considered.
    !> \param[out]  lAddPhase        A logical scalar indicating whether the phase should be added (i.e., TRUE)
    !                                 to the system or not (i.e., FALSE).
    !
    ! iFirst            Absolute index of first species in this phase.
    ! iLast             Absolute index of last species in this phase.
    ! nConstituents     Total number of constituents in this phase.
    ! lSubMinConverged  A logical scalar indicating convergence (true).
    ! dMinMoleFraction  A double real scalar representing the minimum mole fraction to be used in this routine.
    ! dMaxMoleFraction  A double real scalar representing the maximum mole fraction to be used in this routine.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckMiscibilityGap(iSolnPhaseIndex,lAddPhase)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::  i, iFirst, iLast, iFirstOther, iLastOther, nConstituents, iSolnPhaseIndex, iOther, iMaxOther
    real(8)::  dMinMoleFraction, dMaxMoleFraction, dBestDrivingForce, dSum
    logical::  lAddPhase, lPhasePass
    real(8), allocatable :: dBestMolFraction(:)


    ! Initialize variables:
    iterLastMiscGapCheck = iterGlobal
    iFirst           = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    iLast            = nSpeciesPhase(iSolnPhaseIndex)
    nConstituents    = iLast - iFirst + 1
    dMinMoleFraction = 1D-3
    dMaxMoleFraction = 1D0 - dMinMoleFraction * DFLOAT(nConstituents-1)
    dMaxMoleFraction = DMAX1(dMaxMoleFraction, 0.9D0)
    lAddPhase        = .FALSE.
    dBestDrivingForce = 0D0

    if (allocated(dBestMolFraction)) deallocate(dBestMolFraction)
    allocate(dBestMolFraction(nConstituents))
    dBestMolFraction = 0D0

    ! First probe an interior composition.  Corner-only starts are a poor match
    ! for liquid minima that are stabilized by mixing and sit away from the
    ! simplex extremums.
    dMolFraction(iFirst:iLast) = 1D0 / DFLOAT(nConstituents)
    lPhasePass = .FALSE.
    call Subminimization(iSolnPhaseIndex, lPhasePass)

    if (lPhasePass) then
        if (dDrivingForceSoln(iSolnPhaseIndex) < dBestDrivingForce) then
            dBestDrivingForce = dDrivingForceSoln(iSolnPhaseIndex)
            dBestMolFraction = dMolFraction(iFirst:iLast)
            lAddPhase = .TRUE.
        end if
    end if

    ! For duplicate liquid phases, also probe a perturbed version of the
    ! companion liquid's current composition.  The second basin is often close to
    ! the first stable liquid rather than at a simplex corner.
    if (lMiscibility(iSolnPhaseIndex)) then
        iOther = 0
        if ((iSolnPhaseIndex > 1) .AND. (cSolnPhaseName(iSolnPhaseIndex) == cSolnPhaseName(iSolnPhaseIndex - 1))) then
            iOther = iSolnPhaseIndex - 1
        else
            do i = 1, nSolnPhasesSys
                if (i == iSolnPhaseIndex) cycle
                if (cSolnPhaseName(i) == cSolnPhaseName(iSolnPhaseIndex)) then
                    iOther = i
                    exit
                end if
            end do
        end if

        if (iOther > 0) then
            iFirstOther = nSpeciesPhase(iOther - 1) + 1
            iLastOther  = nSpeciesPhase(iOther)
            if (iLastOther - iFirstOther == iLast - iFirst) then
                dMolFraction(iFirst:iLast) = dMolFraction(iFirstOther:iLastOther)
                iMaxOther = MAXVAL(MAXLOC(dMolFraction(iFirst:iLast)))
                dMolFraction(iFirst + iMaxOther - 1) = DMAX1(dMinMoleFraction, 0.5D0 * dMolFraction(iFirst + iMaxOther - 1))
                dSum = SUM(dMolFraction(iFirst:iLast))
                if (dSum > 0D0) dMolFraction(iFirst:iLast) = dMolFraction(iFirst:iLast) / dSum

                lPhasePass = .FALSE.
                call Subminimization(iSolnPhaseIndex, lPhasePass)

                if (lPhasePass) then
                    if (dDrivingForceSoln(iSolnPhaseIndex) < dBestDrivingForce) then
                        dBestDrivingForce = dDrivingForceSoln(iSolnPhaseIndex)
                        dBestMolFraction = dMolFraction(iFirst:iLast)
                        lAddPhase = .TRUE.
                    end if
                end if
            end if
        end if
    end if

    ! Perform subminimization multiple times by initializing from all extremums of the domain space:
    LOOP_Constituents: do i = 1, nConstituents

        ! Initialize the mole fractions:
        dMolFraction(iFirst:iLast) = dMinMoleFraction
        dMolFraction(iFirst+i-1)   = dMaxMoleFraction

        ! Perform subminimization:
        lPhasePass = .FALSE.
        call Subminimization(iSolnPhaseIndex, lPhasePass)

        if (lPhasePass) then
            if (dDrivingForceSoln(iSolnPhaseIndex) < dBestDrivingForce) then
                dBestDrivingForce = dDrivingForceSoln(iSolnPhaseIndex)
                dBestMolFraction = dMolFraction(iFirst:iLast)
                lAddPhase = .TRUE.
            end if
        end if

    end do LOOP_Constituents

    if (lAddPhase) then
        dMolFraction(iFirst:iLast) = dBestMolFraction
        dDrivingForceSoln(iSolnPhaseIndex) = dBestDrivingForce
    end if

    if (allocated(dBestMolFraction)) deallocate(dBestMolFraction)

    return

end subroutine CheckMiscibilityGap

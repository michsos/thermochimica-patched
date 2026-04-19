

!-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataBlockQKTO.f90
    !> \brief   Parse the data block section corresponding to a QKTO phase of a ChemSage data-file.
    !> \author  M.H.A. Piro
    !> \date    Dec. 21, 2012
    !> \sa      ParseCSDataFile.f90
    !> \sa      ParseCSDataBlock.f90
    !> \sa      ParseCSDataBlockGibbs.f90
    !> \sa      ParseCSInterpolationOverrides.f90
    !
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified and is completely
    ! independent of ChemApp and related products, including Solgas, Solgasmix, Fact, FactSage
    ! and ChemSage.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !   10/06/2011      M.H.A. Piro     Original code
    !   12/21/2012      M.H.A. Piro     Relocate to independent f90 file.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to parse the "data block" section of a ChemSage data-file
    !! containing a "QKTO" phase (Quasi-chemical Kohler-TOop).  The molar excess Gibbs energy of mixing is
    !! generally given as \f$ \Delta g_i^{ex} = \sum_{p=1} x_i^m x_j^n (^pL) \f$, where i and j are the
    !! constituent indices, m and n are the exponents (provided by the model) and \f$ ^pL \f$ is the
    !! temperature dependent mixing parameter.
    !!
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO                      A scalar integer that indicates a successful exit or identifies an error.
    ! nSpeciesCS                Number of species in the system (combined solution species and pure
    !                            separate phases).
    ! nGibbsEqSpecies           Number of Gibbs energy equations for a particular species.
    ! iSpeciesAtomsCS           Integer matrix representing the number of atoms of a particular
    !                            elements in a species.
    ! iParticlesPerMoleCS       An integer vector containing the number of particles per mole of the
    !                            constituent species formula mass.  The default value is 1.
    ! cSolnPhaseNameCS          The name of a solution phase.
    ! cSolnPhaseTypeCS          The type of a solution phase.
    ! cSolnPhaseTypeSupport     A character array representing solution phase types that are supported.
    !
!-------------------------------------------------------------------------------------------------------------


subroutine ParseCSDataBlockQKTO( i )

    USE ModuleParseCS

    implicit none

    integer :: i, j, k

    ! QKTOM phases include magnetic mixing terms before non-ideal mixing terms,
    ! analogous to RKMPM in ParseCSDataBlockRKMP.
    ! The end of the list of magnetic mixing terms is indicated by a "0".
    if (cSolnPhaseTypeCS(i) == 'QKTOM') then
        LOOP_MagneticMixingQKTO: do

            ! Grow arrays if capacity is exceeded:
            if (nMagParamCS + 1 > size(iMagneticParamCS, 1)) call GrowMagneticParamArrays

            ! Read in number of constituents involved in parameter:
            read (1,*,IOSTAT = INFO) iMagneticParamCS(nMagParamCS+1,1)

            ! The end of the parameter listing is marked by "0":
            if (iMagneticParamCS(nMagParamCS+1,1) == 0) exit LOOP_MagneticMixingQKTO

            ! Check if the parameter is binary or ternary:
            if (iMagneticParamCS(nMagParamCS+1,1) == 2) then

                ! Binary magnetic mixing terms:
                read (1,*,IOSTAT = INFO) iMagneticParamCS(nMagParamCS+1,2:4)

                k = iMagneticParamCS(nMagParamCS+1,4)

                do j = 1, k
                    nMagParamCS = nMagParamCS + 1
                    iMagneticParamCS(nMagParamCS,1:3) = iMagneticParamCS(nMagParamCS-j+1,1:3)
                    iMagneticParamCS(nMagParamCS,4)   = j - 1
                    read (1,*,IOSTAT = INFO) dMagneticParamCS(nMagParamCS,1:2)
                end do

            elseif (iMagneticParamCS(nMagParamCS+1,1) == 3) then

                ! Ternary magnetic mixing terms:
                read (1,*,IOSTAT = INFO) iMagneticParamCS(nMagParamCS+1,2:5)

                k = iMagneticParamCS(nMagParamCS+1,5)

                do j = 1, k
                    nMagParamCS = nMagParamCS + 1
                    iMagneticParamCS(nMagParamCS,1:5) = iMagneticParamCS(nMagParamCS-j+1,1:5)
                    iMagneticParamCS(nMagParamCS,5)   = iMagneticParamCS(nMagParamCS,1+j)
                    read (1,*,IOSTAT = INFO) dMagneticParamCS(nMagParamCS,1:2)
                end do

            elseif (iMagneticParamCS(nMagParamCS+1,1) == 4) then

                ! Quaternary magnetic mixing terms:
                read (1,*,IOSTAT = INFO) iMagneticParamCS(nMagParamCS+1,2:6)

                k = iMagneticParamCS(nMagParamCS+1,6)

                do j = 1, k
                    nMagParamCS = nMagParamCS + 1
                    iMagneticParamCS(nMagParamCS,1:6) = iMagneticParamCS(nMagParamCS-j+1,1:6)
                    iMagneticParamCS(nMagParamCS,6)   = iMagneticParamCS(nMagParamCS,1+j)
                    read (1,*,IOSTAT = INFO) dMagneticParamCS(nMagParamCS,1:2)
                end do

            else
                ! This parameter is not recognized; record an error.
                INFO = 43
                return
            end if

        end do LOOP_MagneticMixingQKTO
    end if

    ! Loop through excess parameters:
    LOOP_ExcessMixingQKTO: do

        ! Grow arrays if capacity is exceeded:
        if (nParamCS + 1 > size(iRegularParamCS, 1)) call GrowRegularParamArrays

        ! Read in number of constituents involved in parameter:
        read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,1)

        ! The end of the parameter listing is marked by "0",
        ! or a negative number indicating the number of extra parameter lines.
        ! These lines indicate overwriting of default interpolation schemes.
        if (iRegularParamCS(nParamCS+1,1) <= 0) then
            call ParseCSInterpolationOverrides(i)
            exit LOOP_ExcessMixingQKTO
        end if

        nParamCS = nParamCS + 1

        if (iRegularParamCS(nParamCS,1) == 2) then
            ! Binary mixing terms:
            read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS,2:5), dRegularParamCS(nParamCS,1:6)
        elseif (iRegularParamCS(nParamCS,1) == 3) then
            ! Ternary mixing terms:
            read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS,2:7), dRegularParamCS(nParamCS,1:6)
        else
            ! This is not recognized.
            INFO = 1600 + i
            return
        end if

    end do LOOP_ExcessMixingQKTO

    ! Report an error if necessary:
    if (INFO /= 0) then
        INFO = 1600 + i
        return
    end if

    return

end subroutine ParseCSDataBlockQKTO

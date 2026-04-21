!> @file DumpDatabaseJSON.f90
!> @brief Serialize parsed ChemSage database to JSON for external solvers.
!>
!> After ParseCSDataFile has populated ModuleParseCS globals, this routine
!> writes all structural data (elements, species, phases, Gibbs coefficients,
!> sublattice structures, excess parameters, magnetic parameters, coordination
!> numbers, etc.) to a JSON file that Julia or other consumers can load.

subroutine DumpDatabaseJSON(cOutPath)

    USE ModuleParseCS

    implicit none

    character(*), intent(in) :: cOutPath
    integer :: iUnit, i, j, k, s, m, n, p
    integer :: iPhaseForSub(MAX(1,nCountSublatticeCS))
    integer :: nSp, iGStart, nInts, iSolnIdx, nMeaningful
    character(1) :: cTypeChar

    iUnit = 99
    open(unit=iUnit, file=TRIM(cOutPath), status='replace', action='write')

    ! ---------------------------------------------------------------
    ! Top-level object
    ! ---------------------------------------------------------------
    write(iUnit, '(A)') '{'

    ! Scalars
    write(iUnit, '(A,I0,A)') '  "n_elements": ', nElementsCS, ','
    write(iUnit, '(A,I0,A)') '  "n_species": ', nSpeciesCS, ','
    write(iUnit, '(A,I0,A)') '  "n_solution_phases": ', nSolnPhasesSysCS, ','
    write(iUnit, '(A,I0,A)') '  "n_params": ', nParamCS, ','
    write(iUnit, '(A,I0,A)') '  "n_mag_params": ', nMagParamCS, ','
    write(iUnit, '(A,I0,A)') '  "n_count_sublattice": ', nCountSublatticeCS, ','

    ! ---------------------------------------------------------------
    ! Element names
    ! ---------------------------------------------------------------
    call jsWriteStrArr3(iUnit, 'element_names', cElementNameCS, nElementsCS, .true.)

    ! ---------------------------------------------------------------
    ! Atomic masses
    ! ---------------------------------------------------------------
    call jsWriteRealArr(iUnit, 'atomic_masses', dAtomicMassCS, nElementsCS, .true.)

    ! ---------------------------------------------------------------
    ! Species data (all nSpeciesCS species in order)
    ! ---------------------------------------------------------------
    call jsWriteStrArr128(iUnit, 'species_names', cSpeciesNameCS, nSpeciesCS, .true.)
    call jsWriteIntArr(iUnit, 'species_phase', iPhaseCS, nSpeciesCS, .true.)
    call jsWriteIntArr(iUnit, 'particles_per_mole', iParticlesPerMoleCS, &
                       nSpeciesCS, .true.)

    ! Stoichiometry matrix [nSpecies][nElements]
    write(iUnit, '(A)') '  "stoichiometry": ['
    do j = 1, nSpeciesCS
        write(iUnit, '(A)', advance='no') '    ['
        do k = 1, nElementsCS
            call jsWriteReal(iUnit, dStoichSpeciesCS(j, k), k < nElementsCS)
        end do
        if (j < nSpeciesCS) then
            write(iUnit, '(A)') '],'
        else
            write(iUnit, '(A)') ']'
        end if
    end do
    write(iUnit, '(A)') '  ],'

    ! Number of Gibbs intervals per species
    call jsWriteIntArr(iUnit, 'n_gibbs_intervals', nGibbsEqSpecies, nSpeciesCS, .true.)

    ! Gibbs coefficients per species: nested array [nSpecies][nIntervals][21]
    write(iUnit, '(A)') '  "gibbs_coefficients": ['
    iGStart = 0
    do j = 1, nSpeciesCS
        nInts = nGibbsEqSpecies(j)
        write(iUnit, '(A)') '    ['
        do i = 1, nInts
            write(iUnit, '(A)', advance='no') '      ['
            do k = 1, nGibbsCoeff
                call jsWriteReal(iUnit, dGibbsCoeffSpeciesTemp(k, iGStart + i), &
                                 k < nGibbsCoeff)
            end do
            if (i < nInts) then
                write(iUnit, '(A)') '],'
            else
                write(iUnit, '(A)') ']'
            end if
        end do
        iGStart = iGStart + nInts
        if (j < nSpeciesCS) then
            write(iUnit, '(A)') '    ],'
        else
            write(iUnit, '(A)') '    ]'
        end if
    end do
    write(iUnit, '(A)') '  ],'

    ! Gibbs magnetic data [nSpecies][4]
    write(iUnit, '(A)') '  "gibbs_magnetic": ['
    do j = 1, nSpeciesCS
        write(iUnit, '(A)', advance='no') '    ['
        do k = 1, 4
            call jsWriteReal(iUnit, dGibbsMagneticCS(j, k), k < 4)
        end do
        if (j < nSpeciesCS) then
            write(iUnit, '(A)') '],'
        else
            write(iUnit, '(A)') ']'
        end if
    end do
    write(iUnit, '(A)') '  ],'

    ! ---------------------------------------------------------------
    ! Solution phase metadata
    ! ---------------------------------------------------------------
    call jsWriteStrArr128(iUnit, 'soln_phase_names', cSolnPhaseNameCS, &
                          nSolnPhasesSysCS, .true.)
    call jsWriteStrArr8(iUnit, 'soln_phase_types', cSolnPhaseTypeCS, &
                        nSolnPhasesSysCS, .true.)

    ! Cumulative species counts (include the 0-th element)
    write(iUnit, '(A)', advance='no') '  "n_species_phase": ['
    do i = 0, nSolnPhasesSysCS
        write(iUnit, '(I0)', advance='no') nSpeciesPhaseCS(i)
        if (i < nSolnPhasesSysCS) write(iUnit, '(A)', advance='no') ', '
    end do
    write(iUnit, '(A)') '],'

    ! Phase-to-sublattice index mapping
    call jsWriteIntArr(iUnit, 'phase_sublattice_idx', iPhaseSublatticeCS, &
                       nSolnPhasesSysCS, .true.)

    ! nPairsSRO per sublattice phase [nCountSublattice][2]
    write(iUnit, '(A)') '  "n_pairs_sro": ['
    do s = 1, nCountSublatticeCS
        write(iUnit, '(A,I0,A,I0,A)', advance='no') '    [', &
            nPairsSROCS(s, 1), ', ', nPairsSROCS(s, 2), ']'
        if (s < nCountSublatticeCS) then
            write(iUnit, '(A)') ','
        else
            write(iUnit, '(A)') ''
        end if
    end do
    write(iUnit, '(A)') '  ],'

    ! ---------------------------------------------------------------
    ! Cumulative parameter counts per phase
    ! ---------------------------------------------------------------
    write(iUnit, '(A)', advance='no') '  "n_param_phase": ['
    do i = 0, nSolnPhasesSysCS
        write(iUnit, '(I0)', advance='no') nParamPhaseCS(i)
        if (i < nSolnPhasesSysCS) write(iUnit, '(A)', advance='no') ', '
    end do
    write(iUnit, '(A)') '],'

    ! ---------------------------------------------------------------
    ! Regular (excess) mixing parameters
    ! ---------------------------------------------------------------
    write(iUnit, '(A)') '  "regular_params": ['
    do p = 1, nParamCS
        write(iUnit, '(A)') '    {'
        ! Write indices — only meaningful positions, zero the rest
        ! Position 1 = nConstituents (j), 2..j+1 = constituent indices,
        ! j+2 = mixing order. For SUBG/SUBQ: positions 2:9 + 10:11.
        write(iUnit, '(A)', advance='no') '      "indices": ['
        nMeaningful = iRegularParamCS(p, 1) + 2  ! j + 2 meaningful positions
        if (nMeaningful > 11) nMeaningful = 11
        ! For SUBQ/SUBG, all 11 positions are meaningful
        if (cRegularParamCS(p) == 'Q' .OR. cRegularParamCS(p) == 'G' .OR. &
            cRegularParamCS(p) == 'R' .OR. cRegularParamCS(p) == 'B') then
            nMeaningful = 11
        end if
        do k = 1, 11
            if (k <= nMeaningful) then
                write(iUnit, '(I0)', advance='no') iRegularParamCS(p, k)
            else
                write(iUnit, '(A)', advance='no') '0'
            end if
            if (k < 11) write(iUnit, '(A)', advance='no') ', '
        end do
        write(iUnit, '(A)') '],'
        write(iUnit, '(A)', advance='no') '      "coeffs": ['
        do k = 1, 6
            call jsWriteReal(iUnit, dRegularParamCS(p, k), k < 6)
        end do
        write(iUnit, '(A)') '],'
        ! Sanitize null bytes — SUBL params don't set this field
        cTypeChar = cRegularParamCS(p)
        if (IACHAR(cTypeChar) < 32 .OR. IACHAR(cTypeChar) > 126) cTypeChar = ' '
        write(iUnit, '(A,A,A)') '      "type": "', cTypeChar, '"'
        if (p < nParamCS) then
            write(iUnit, '(A)') '    },'
        else
            write(iUnit, '(A)') '    }'
        end if
    end do
    write(iUnit, '(A)') '  ],'

    ! ---------------------------------------------------------------
    ! Magnetic mixing parameters
    ! ---------------------------------------------------------------
    write(iUnit, '(A)', advance='no') '  "n_mag_param_phase": ['
    do i = 0, nSolnPhasesSysCS
        write(iUnit, '(I0)', advance='no') nMagParamPhaseCS(i)
        if (i < nSolnPhasesSysCS) write(iUnit, '(A)', advance='no') ', '
    end do
    write(iUnit, '(A)') '],'

    write(iUnit, '(A)') '  "magnetic_params": ['
    do p = 1, nMagParamCS
        write(iUnit, '(A)') '    {'
        write(iUnit, '(A)', advance='no') '      "indices": ['
        nMeaningful = iMagneticParamCS(p, 1) + 2
        if (nMeaningful > 11) nMeaningful = 11
        do k = 1, 11
            if (k <= nMeaningful) then
                write(iUnit, '(I0)', advance='no') iMagneticParamCS(p, k)
            else
                write(iUnit, '(A)', advance='no') '0'
            end if
            if (k < 11) write(iUnit, '(A)', advance='no') ', '
        end do
        write(iUnit, '(A)') '],'
        write(iUnit, '(A)', advance='no') '      "coeffs": ['
        call jsWriteReal(iUnit, dMagneticParamCS(p, 1), .true.)
        call jsWriteReal(iUnit, dMagneticParamCS(p, 2), .false.)
        write(iUnit, '(A)') ']'
        if (p < nMagParamCS) then
            write(iUnit, '(A)') '    },'
        else
            write(iUnit, '(A)') '    }'
        end if
    end do
    write(iUnit, '(A)') '  ],'

    ! ---------------------------------------------------------------
    ! Sublattice data (one entry per sublattice phase)
    ! ---------------------------------------------------------------
    ! Build reverse map: sublattice index -> solution phase index
    iPhaseForSub = 0
    do i = 1, nSolnPhasesSysCS
        if (iPhaseSublatticeCS(i) > 0) then
            iPhaseForSub(iPhaseSublatticeCS(i)) = i
        end if
    end do

    write(iUnit, '(A)') '  "sublattice_phases": ['
    do s = 1, nCountSublatticeCS
        iSolnIdx = iPhaseForSub(s)
        nSp = nSpeciesPhaseCS(iSolnIdx) - nSpeciesPhaseCS(iSolnIdx - 1)

        write(iUnit, '(A)') '    {'
        write(iUnit, '(A,I0,A)') '      "soln_phase_index": ', iSolnIdx, ','
        write(iUnit, '(A,I0,A)') '      "n_sublattices": ', nSublatticePhaseCS(s), ','

        ! Site ratios
        n = nSublatticePhaseCS(s)
        write(iUnit, '(A)', advance='no') '      "site_ratios": ['
        do k = 1, n
            call jsWriteReal(iUnit, dStoichSublatticeCS(s, k), k < n)
        end do
        write(iUnit, '(A)') '],'

        ! Constituent counts per sublattice
        write(iUnit, '(A)', advance='no') '      "constituent_counts": ['
        do k = 1, n
            write(iUnit, '(I0)', advance='no') nConstituentSublatticeCS(s, k)
            if (k < n) write(iUnit, '(A)', advance='no') ', '
        end do
        write(iUnit, '(A)') '],'

        ! Constituent names [nSublattices][nConstituents]
        write(iUnit, '(A)') '      "constituent_names": ['
        do m = 1, n
            write(iUnit, '(A)', advance='no') '        ['
            do k = 1, nConstituentSublatticeCS(s, m)
                write(iUnit, '(A)', advance='no') '"'
                write(iUnit, '(A)', advance='no') &
                    TRIM(ADJUSTL(cConstituentNameSUBCS(s, m, k)))
                write(iUnit, '(A)', advance='no') '"'
                if (k < nConstituentSublatticeCS(s, m)) &
                    write(iUnit, '(A)', advance='no') ', '
            end do
            if (m < n) then
                write(iUnit, '(A)') '],'
            else
                write(iUnit, '(A)') ']'
            end if
        end do
        write(iUnit, '(A)') '      ],'

        ! Charges [nSublattices][nConstituents]
        write(iUnit, '(A)') '      "charges": ['
        do m = 1, n
            write(iUnit, '(A)', advance='no') '        ['
            do k = 1, nConstituentSublatticeCS(s, m)
                call jsWriteReal(iUnit, dSublatticeChargeCS(s, m, k), &
                                 k < nConstituentSublatticeCS(s, m))
            end do
            if (m < n) then
                write(iUnit, '(A)') '],'
            else
                write(iUnit, '(A)') ']'
            end if
        end do
        write(iUnit, '(A)') '      ],'

        ! Chemical groups [nSublattices][nConstituents]
        write(iUnit, '(A)') '      "chemical_groups": ['
        do m = 1, n
            write(iUnit, '(A)', advance='no') '        ['
            do k = 1, nConstituentSublatticeCS(s, m)
                write(iUnit, '(I0)', advance='no') iChemicalGroupCS(s, m, k)
                if (k < nConstituentSublatticeCS(s, m)) &
                    write(iUnit, '(A)', advance='no') ', '
            end do
            if (m < n) then
                write(iUnit, '(A)') '],'
            else
                write(iUnit, '(A)') ']'
            end if
        end do
        write(iUnit, '(A)') '      ],'

        ! Constituent mapping [nSublattices][nSpeciesInPhase]
        write(iUnit, '(A)') '      "constituent_map": ['
        do m = 1, n
            write(iUnit, '(A)', advance='no') '        ['
            do k = 1, nSp
                write(iUnit, '(I0)', advance='no') iConstituentSublatticeCS(s, m, k)
                if (k < nSp) write(iUnit, '(A)', advance='no') ', '
            end do
            if (m < n) then
                write(iUnit, '(A)') '],'
            else
                write(iUnit, '(A)') ']'
            end if
        end do
        write(iUnit, '(A)') '      ],'

        ! Coordination numbers [nSpeciesInPhase][4]
        write(iUnit, '(A)') '      "coordination_numbers": ['
        do k = 1, nSp
            write(iUnit, '(A)', advance='no') '        ['
            do m = 1, 4
                call jsWriteReal(iUnit, dCoordinationNumberCS(iSolnIdx, k, m), m < 4)
            end do
            if (k < nSp) then
                write(iUnit, '(A)') '],'
            else
                write(iUnit, '(A)') ']'
            end if
        end do
        write(iUnit, '(A)') '      ],'

        ! Pair IDs [nSpeciesInPhase][4]
        write(iUnit, '(A)') '      "pair_ids": ['
        do k = 1, nSp
            write(iUnit, '(A)', advance='no') '        ['
            do m = 1, 4
                write(iUnit, '(I0)', advance='no') iPairIDCS(iSolnIdx, k, m)
                if (m < 4) write(iUnit, '(A)', advance='no') ', '
            end do
            if (k < nSp) then
                write(iUnit, '(A)') '],'
            else
                write(iUnit, '(A)') ']'
            end if
        end do
        write(iUnit, '(A)') '      ],'

        ! Constituent coefficients [nSpeciesInPhase][5]
        write(iUnit, '(A)') '      "constituent_coefficients": ['
        do k = 1, nSp
            write(iUnit, '(A)', advance='no') '        ['
            do m = 1, 5
                call jsWriteReal(iUnit, dConstituentCoefficientsCS(s, k, m), m < 5)
            end do
            if (k < nSp) then
                write(iUnit, '(A)') '],'
            else
                write(iUnit, '(A)') ']'
            end if
        end do
        write(iUnit, '(A)') '      ],'

        ! Zeta values [nSpeciesInPhase]
        write(iUnit, '(A)', advance='no') '      "zeta": ['
        do k = 1, nSp
            call jsWriteReal(iUnit, dZetaSpeciesCS(s, k), k < nSp)
        end do
        write(iUnit, '(A)') '],'

        ! Pair stoichiometry [nSpeciesInPhase][nElements]
        write(iUnit, '(A)') '      "stoich_pairs": ['
        do k = 1, nSp
            write(iUnit, '(A)', advance='no') '        ['
            do m = 1, nElementsCS
                call jsWriteReal(iUnit, dStoichPairsCS(s, k, m), m < nElementsCS)
            end do
            if (k < nSp) then
                write(iUnit, '(A)') '],'
            else
                write(iUnit, '(A)') ']'
            end if
        end do
        write(iUnit, '(A)') '      ]'

        if (s < nCountSublatticeCS) then
            write(iUnit, '(A)') '    },'
        else
            write(iUnit, '(A)') '    }'
        end if
    end do
    write(iUnit, '(A)') '  ]'

    ! Close top-level object
    write(iUnit, '(A)') '}'

    close(iUnit)

contains

    ! ===============================================================
    ! Helper: write a real value with optional trailing comma
    ! ===============================================================
    subroutine jsWriteReal(iu, dVal, lComma)
        integer, intent(in) :: iu
        real(8), intent(in) :: dVal
        logical, intent(in) :: lComma
        character(30) :: cBuf
        write(cBuf, '(ES23.15E3)') dVal
        write(iu, '(A)', advance='no') TRIM(ADJUSTL(cBuf))
        if (lComma) write(iu, '(A)', advance='no') ', '
    end subroutine jsWriteReal

    ! ===============================================================
    ! Helper: write 1D integer array
    ! ===============================================================
    subroutine jsWriteIntArr(iu, cKey, iArr, nn, lComma)
        integer, intent(in) :: iu, nn
        integer, intent(in) :: iArr(:)
        character(*), intent(in) :: cKey
        logical, intent(in) :: lComma
        integer :: kk
        write(iu, '(A,A,A)', advance='no') '  "', TRIM(cKey), '": ['
        do kk = 1, nn
            write(iu, '(I0)', advance='no') iArr(kk)
            if (kk < nn) write(iu, '(A)', advance='no') ', '
        end do
        if (lComma) then
            write(iu, '(A)') '],'
        else
            write(iu, '(A)') ']'
        end if
    end subroutine jsWriteIntArr

    ! ===============================================================
    ! Helper: write 1D real(8) array
    ! ===============================================================
    subroutine jsWriteRealArr(iu, cKey, dArr, nn, lComma)
        integer, intent(in) :: iu, nn
        real(8), intent(in) :: dArr(:)
        character(*), intent(in) :: cKey
        logical, intent(in) :: lComma
        integer :: kk
        write(iu, '(A,A,A)', advance='no') '  "', TRIM(cKey), '": ['
        do kk = 1, nn
            call jsWriteReal(iu, dArr(kk), kk < nn)
        end do
        if (lComma) then
            write(iu, '(A)') '],'
        else
            write(iu, '(A)') ']'
        end if
    end subroutine jsWriteRealArr

    ! ===============================================================
    ! Helper: write 1D character(3) array
    ! ===============================================================
    subroutine jsWriteStrArr3(iu, cKey, cArr, nn, lComma)
        integer, intent(in) :: iu, nn
        character(3), intent(in) :: cArr(:)
        character(*), intent(in) :: cKey
        logical, intent(in) :: lComma
        integer :: kk
        write(iu, '(A,A,A)', advance='no') '  "', TRIM(cKey), '": ['
        do kk = 1, nn
            write(iu, '(A)', advance='no') '"'
            write(iu, '(A)', advance='no') TRIM(ADJUSTL(cArr(kk)))
            write(iu, '(A)', advance='no') '"'
            if (kk < nn) write(iu, '(A)', advance='no') ', '
        end do
        if (lComma) then
            write(iu, '(A)') '],'
        else
            write(iu, '(A)') ']'
        end if
    end subroutine jsWriteStrArr3

    ! ===============================================================
    ! Helper: write 1D character(128) array
    ! ===============================================================
    subroutine jsWriteStrArr128(iu, cKey, cArr, nn, lComma)
        integer, intent(in) :: iu, nn
        character(128), intent(in) :: cArr(:)
        character(*), intent(in) :: cKey
        logical, intent(in) :: lComma
        integer :: kk
        write(iu, '(A,A,A)') '  "', TRIM(cKey), '": ['
        do kk = 1, nn
            write(iu, '(A)', advance='no') '    "'
            write(iu, '(A)', advance='no') TRIM(ADJUSTL(cArr(kk)))
            write(iu, '(A)', advance='no') '"'
            if (kk < nn) then
                write(iu, '(A)') ','
            else
                write(iu, '(A)') ''
            end if
        end do
        if (lComma) then
            write(iu, '(A)') '  ],'
        else
            write(iu, '(A)') '  ]'
        end if
    end subroutine jsWriteStrArr128

    ! ===============================================================
    ! Helper: write 1D character(8) array
    ! ===============================================================
    subroutine jsWriteStrArr8(iu, cKey, cArr, nn, lComma)
        integer, intent(in) :: iu, nn
        character(8), intent(in) :: cArr(:)
        character(*), intent(in) :: cKey
        logical, intent(in) :: lComma
        integer :: kk
        write(iu, '(A,A,A)', advance='no') '  "', TRIM(cKey), '": ['
        do kk = 1, nn
            write(iu, '(A)', advance='no') '"'
            write(iu, '(A)', advance='no') TRIM(ADJUSTL(cArr(kk)))
            write(iu, '(A)', advance='no') '"'
            if (kk < nn) write(iu, '(A)', advance='no') ', '
        end do
        if (lComma) then
            write(iu, '(A)') '],'
        else
            write(iu, '(A)') ']'
        end if
    end subroutine jsWriteStrArr8

end subroutine DumpDatabaseJSON

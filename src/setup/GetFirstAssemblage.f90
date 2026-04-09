
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GetFirstAssemblage.f90
    !> \brief   Determine the first phase assemblage for testing.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      LevelingSolver.f90
    !
    !
    ! References:
    ! ===========
    !
    ! For further information regarding this software, refer to the following material:
    !
    !        M.H.A. Piro, "Computation of Thermodynamic Equilibria Pertinent to Nuclear Materials
    !        in Multi-Physics Codes," PhD Dissertation, Royal Military College of Canada, 2011.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code
    !   07/31/2011      M.H.A. Piro         Clean up code: remove unnecessary variables, update variable names
    !   10/21/2011      M.H.A. Piro         Clean up code: modules and simplify
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to estimate the very first phase assemblage for the Leveling
    !! subroutine.  Only the pure species will be considered here, since any possible combination of
    !! the pure species will yield a positive number of moles.  In more mathematical terms, the stoichiometry
    !! matrix is a diagonal matrix when the assemblage is comprised of only pure species.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! iAssemblage               Integer vector containing the indices of species estimated to be part of
    !                            the equilibrium phase assemblage.
    ! dAtomFractionSpecies      Double matrix representing the atomic fraction of species
    ! dChemicalPotential        The estimated chemical potential of each species.
    ! dLevel                    The adjustment to the element potentials as a result of leveling.
    ! iAtomFractionSpecies      A temporary integer vector representing the corresponding integer from
    !                            dAtomFractionSpecies.
    ! dTempArr                  A temporary double real array representing the species as rows and
    !                            elements as columns.  This variable will equal zero everywhere except
    !                            for pure species.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine GetFirstAssemblage

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer::                               i, j, k, kDormant
    integer,dimension(nSpecies,nElements):: iAtomFractionSpecies
    real(8),dimension(nSpecies,nElements):: dTempArr
    logical                                :: IsSpeciesDormant


    ! The following is a temporary integer matrix with equivalent dimensions to dAtomFractionSpecies
    ! (a real(8) matrix).  Since all coefficients of dAtomFractionSpecies are between 0 and 1, the
    ! conversion of coefficients to an integer results in converting all coefficients to either 0 or 1.
    ! Therefore, the only non-zero coefficients of the following variable correspond to pure species.

    iAtomFractionSpecies = INT(dAtomFractionSpecies)

    ! The temporary real(8) matrix dTempArr has the same dimensions as iAtomFractionSpecies.  All
    ! coefficients of this matrix are zero except for the rows representing pure species and the
    ! column corresponding to its respective element.  Therefore, each column (representing an element)
    ! will only contain the standard Gibbs energy of a pure species.

    dTempArr = 0D0
    do i = 1,nSpecies
        dTempArr(i,1:nElements) = dChemicalPotential(i) * DFLOAT(iAtomFractionSpecies(i,1:nElements))
    end do

    ! The following integer vector contains the indices of species with the
    ! lowest standard Gibbs energy for its respective element.  It is possible
    ! for the lowest standard Gibbs energy of a pure species from the database
    ! to be positive.  The MINLOC function would therefore return the index of
    ! the first zero coefficient for a particular column instead of the pure
    ! species with a positive standard Gibbs energy.  The MASK feature is utilized
    ! to ensure that pure species with a positive standard Gibbs energy are considered.

    iAssemblage = 0
    do j = 1, nElements
        k = 0
        kDormant = 0
        do i = 1, nSpecies
            if (iAtomFractionSpecies(i,j) == 0) cycle
            if (IsSpeciesDormant(i)) then
                if (kDormant == 0) then
                    kDormant = i
                else if (dTempArr(i,j) < dTempArr(kDormant,j)) then
                    kDormant = i
                end if
                cycle
            end if
            if (k == 0) then
                k = i
            else if (dTempArr(i,j) < dTempArr(k,j)) then
                k = i
            end if
        end do
        ! Dormant pure species can still provide a temporary Leveling carrier when
        ! no active pure carrier exists for an element. They are removed again before
        ! the active GEM assemblage is constructed.
        if (k == 0) k = kDormant
        if (k == 0) then
            INFOThermo = 10
            return
        end if
        iAssemblage(j) = k
    end do

    ! Establish the first adjustment to the element potentials:
    do i = 1, nElements
        if (iAssemblage(i) == 0) return
        dLevel(i) = dChemicalPotential(iAssemblage(i))
    end do

    return

end subroutine GetFirstAssemblage

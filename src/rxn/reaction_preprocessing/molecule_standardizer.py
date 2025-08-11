from typing import List, Optional, Tuple

from rxn.chemutils.conversion import canonicalize_smiles
from rxn.chemutils.exceptions import InvalidSmiles
from rxn.chemutils.reaction_equation import ReactionEquation

from rxn.reaction_preprocessing.annotations.missing_annotation_detector import (
    MissingAnnotationDetector,
)
from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    MoleculeAnnotation,
)
from rxn.reaction_preprocessing.annotations.molecule_replacer import MoleculeReplacer
from rxn.reaction_preprocessing.annotations.rejected_molecules_filter import (
    RejectedMoleculesFilter,
)
from rxn.reaction_preprocessing.cleaner import remove_isotope_information


class MoleculeStandardizationError(ValueError):
    """Base class for standardization exceptions."""


class RejectedMolecule(MoleculeStandardizationError):
    """Exception raised when standardizing a molecule annotated as "Rejected"."""

    def __init__(self, smiles: str):
        """
        Args:
            smiles: rejected SMILES string.
        """
        super().__init__(f'Cannot standardize: rejected molecule "{smiles}"')


class MissingAnnotation(MoleculeStandardizationError):
    """Exception raised when standardizing a molecule that should be annotated."""

    def __init__(self, smiles: str):
        """
        Args:
            smiles: rejected SMILES string.
        """
        super().__init__(f'Cannot standardize: molecule "{smiles}" must be annotated.')
        self.smiles = smiles


class MoleculeStandardizer:
    """
    Class to standardize standalone molecules (reactions are standardized with
    the Standardizer class).

    Note that the standardization of one molecule may lead to a combination
    of molecules, hence the functions return lists of strings.
    """

    def __init__(
        self,
        annotations: Optional[List[MoleculeAnnotation]] = None,
        discard_missing_annotations: bool = False,
        canonicalize: bool = True,
    ):
        """
        Args:
            annotations: A list of MoleculeAnnotation objects used to perform
                the substitutions /rejections. Defaults to an empty list.
            discard_missing_annotations: whether reactions containing unannotated
                molecules that should be must be rejected.
            canonicalize: whether to canonicalize the compounds.
        """
        if annotations is None:
            annotations = []

        self.discard_unannotated_metals = discard_missing_annotations
        self.canonicalize = canonicalize

        self.rejection_filter = RejectedMoleculesFilter.from_molecule_annotations(
            annotations
        )
        self.missing_annotation_detector = (
            MissingAnnotationDetector.from_molecule_annotations(annotations)
        )
        self.molecule_replacer = MoleculeReplacer.from_molecule_annotations(annotations)

    def __call__(self, smiles: str) -> List[str]:
        """See doc for standardize()."""
        return self.standardize(smiles)

    def standardize(self, smiles: str) -> List[str]:
        """
        Standardize a molecule.

        The returned value is a list, because in some cases standardization
        returns two independent molecules.

        Args:
            smiles: SMILES string to standardize. Use dots for fragment bonds!

        Raises:
            SanitizationError of one of its subclasses: error in sanitization.
            InvalidSmiles: Invalid SMILES.
            ValueError: "~" being used for fragment bonds.

        Returns:
            Standardized SMILES string.
        """
        if "~" in smiles:
            raise ValueError(f'MoleculeStandardizer must be used without "~": {smiles}')

        # Discard isotope information
        smiles = remove_isotope_information(smiles)

        # Check validity of SMILES (may raise InvalidSmiles), and
        # overwrite if canonicalization required
        canonical_smiles = canonicalize_smiles(smiles)
        if self.canonicalize:
            smiles = canonical_smiles

        # Check for rejected molecules
        if not self.rejection_filter.is_valid_molecule_smiles(smiles):
            raise RejectedMolecule(smiles)

        # Check for non-annotated molecules
        if self.discard_unannotated_metals:
            if self.missing_annotation_detector.molecule_needs_annotation(smiles):
                raise MissingAnnotation(smiles)

        # Replace annotated molecules
        return self.molecule_replacer.replace_molecule_smiles(smiles)

    def standardize_in_equation(self, reaction: ReactionEquation) -> ReactionEquation:
        """
        Do the molecule-wise standardization for a reaction equation.

        Relies on standardize_in_equation_with_errors(), for modularity purposes.
        Will propagate the exceptions raised in that function.
        """

        # Ignoring the lists of SMILES returned in the tuple (which, by construction,
        # will always be empty: if not, an exception will have been raised earlier).
        reaction, *_ = self.standardize_in_equation_with_errors(
            reaction, propagate_exceptions=True
        )
        return reaction

    def standardize_in_equation_with_errors(
        self, reaction: ReactionEquation, propagate_exceptions: bool = False
    ) -> Tuple[ReactionEquation, List[str], List[str], List[str]]:
        """
        Do the molecule-wise standardization for a reaction equation, and get the reasons for
        potential failures.

        This function was originally implemented in Standardizer, and then moved here for more
        modularity.

        Args:
            reaction: reaction to standardize.
            propagate_exceptions: if True, will stop execution and raise directly
                instead of collecting the SMILES leading to the failure. Not ideal,
                but probably the only way (?) to not have duplicated code in the
                function standardize_in_equation().

        Returns:
            Tuple:
                - the standardized reaction equation (or an empty one if there was a failure).
                - list of invalid SMILES in the reaction.
                - list of rejected SMILES in the reaction.
                - list of missing annotations in the reaction.
        """

        missing_annotations = []
        invalid_smiles = []
        rejected_smiles = []

        # Iterate over the reactants, agents, products and update the
        # standardized reaction at the same time
        standardized_reaction = ReactionEquation([], [], [])
        for original_role_group, new_role_group in zip(reaction, standardized_reaction):
            for smiles in original_role_group:
                try:
                    new_role_group.extend(self.standardize(smiles))
                except InvalidSmiles:
                    if propagate_exceptions:
                        raise
                    invalid_smiles.append(smiles)
                except RejectedMolecule:
                    if propagate_exceptions:
                        raise
                    rejected_smiles.append(smiles)
                except MissingAnnotation:
                    if propagate_exceptions:
                        raise
                    missing_annotations.append(smiles)
                except Exception:
                    # in case of unexpected exceptions we log the SMILES
                    # as invalid
                    if propagate_exceptions:
                        raise
                    invalid_smiles.append(smiles)

        # If there was any error: replace by empty reaction equation (">>")
        if invalid_smiles or rejected_smiles or missing_annotations:
            standardized_reaction = ReactionEquation([], [], [])

        return (
            standardized_reaction,
            invalid_smiles,
            rejected_smiles,
            missing_annotations,
        )

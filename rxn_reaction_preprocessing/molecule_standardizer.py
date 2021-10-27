from typing import List
from typing import Optional

from rxn_chemutils.conversion import canonicalize_smiles

from rxn_reaction_preprocessing.annotations.missing_annotation_detector import \
    MissingAnnotationDetector
from rxn_reaction_preprocessing.annotations.molecule_annotation import MoleculeAnnotation
from rxn_reaction_preprocessing.annotations.molecule_replacer import MoleculeReplacer
from rxn_reaction_preprocessing.annotations.rejected_molecules_filter import RejectedMoleculesFilter
from rxn_reaction_preprocessing.cleaner import remove_isotope_information


class MoleculeStandardizationError(ValueError):
    """Base class for standardization exceptions."""


class RejectedMolecule(MoleculeStandardizationError):
    """Exception raised when standardizing a molecule annotated as "Rejected"."""

    def __init__(self, smiles):
        """
        Args:
            smiles: rejected SMILES string.
        """
        super().__init__(f'Cannot standardize: rejected molecule "{smiles}"')


class MissingAnnotation(MoleculeStandardizationError):
    """Exception raised when standardizing a molecule that should be annotated."""

    def __init__(self, smiles):
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
        canonicalize: bool = True
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

        self.rejection_filter = RejectedMoleculesFilter.from_molecule_annotations(annotations)
        self.missing_annotation_detector = MissingAnnotationDetector.from_molecule_annotations(
            annotations
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
        if '~' in smiles:
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

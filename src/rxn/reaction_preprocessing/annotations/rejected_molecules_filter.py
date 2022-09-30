from typing import Iterable, Optional

from rxn.chemutils.reaction_equation import ReactionEquation

from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    AnnotationDecision,
    MoleculeAnnotation,
)


class RejectedMoleculesFilter:
    """
    Class to filter rejected molecules (potentially based on annotations).
    """

    def __init__(self, rejected_molecules: Iterable[str]):
        """
        Args:
            rejected_molecules: Molecules to reject. Use dots for fragment bonds!
        """
        self.rejected_molecules = set(rejected_molecules)

    def is_valid_molecule_smiles(self, smiles: str) -> bool:
        """
        Whether a molecule SMILES is considered to be valid.

        Args:
            smiles: molecule SMILES. Fragment bonds must be given with a dot!
        """
        return smiles not in self.rejected_molecules

    def is_valid_reaction_smiles(
        self, smiles: str, fragment_bond: Optional[str] = None
    ) -> bool:
        """
        Whether a reaction SMILES is considered to be valid.

        Args:
            smiles: reaction SMILES.
            fragment_bond: fragment bond used in the reaction SMILES.
        """
        return self.is_valid_reaction_equation(
            ReactionEquation.from_string(smiles, fragment_bond)
        )

    def is_valid_reaction_equation(self, reaction_equation: ReactionEquation) -> bool:
        """
        Whether a reaction equation is considered to be valid.

        Args:
            reaction_equation: reaction equation instance.
        """
        return all(
            self.is_valid_molecule_smiles(smiles)
            for smiles in reaction_equation.iter_all_smiles()
        )

    @classmethod
    def from_molecule_annotations(
        cls, molecule_annotations: Iterable[MoleculeAnnotation]
    ) -> "RejectedMoleculesFilter":
        """
        Instantiate from existing molecule annotations.

        Args:
            molecule_annotations: existing molecule annotations.
        """
        rejected_molecules = (
            annotation.original_without_fragment_bond
            for annotation in molecule_annotations
            if annotation.decision is AnnotationDecision.REJECT
        )
        return cls(rejected_molecules)

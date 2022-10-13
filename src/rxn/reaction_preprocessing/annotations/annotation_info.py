from typing import Dict, Iterable

from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    AnnotationDecision,
    MoleculeAnnotation,
)


class MoleculeNotAnnotated(ValueError):
    def __init__(self, smiles: str):
        super().__init__(f'The molecule "{smiles}" was not annotated.')


class AnnotationInfo:
    """Get rapid inforamtion about the annotation status of given SMILES."""

    def __init__(self, annotations: Iterable[MoleculeAnnotation]):
        self.annotations: Dict[str, MoleculeAnnotation] = {
            annotation.original_without_fragment_bond: annotation
            for annotation in annotations
        }

    def is_annotated(self, smiles: str) -> bool:
        """
        Whether a molecule SMILES is present in the annotations.

        Args:
            smiles: molecule SMILES. Fragment bonds must be given with a dot!
        """
        return smiles in self.annotations

    def is_accepted(self, smiles: str) -> bool:
        """
        Whether a molecule SMILES has been annotated as accepted.

        Raises:
            MoleculeNotAnnotated: when the given smiles is not in the annotations.

        Args:
            smiles: molecule SMILES. Fragment bonds must be given with a dot!
        """
        try:
            return self.annotations[smiles].decision is AnnotationDecision.ACCEPT
        except KeyError as e:
            raise MoleculeNotAnnotated(smiles) from e

    def is_rejected(self, smiles: str) -> bool:
        """
        Whether a molecule SMILES has been annotated as rejected.

        Raises:
            MoleculeNotAnnotated: when the given smiles is not in the annotations.

        Args:
            smiles: molecule SMILES. Fragment bonds must be given with a dot!
        """
        try:
            return self.annotations[smiles].decision is AnnotationDecision.REJECT
        except KeyError as e:
            raise MoleculeNotAnnotated(smiles) from e

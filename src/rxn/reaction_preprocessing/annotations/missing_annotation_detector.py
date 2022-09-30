from typing import Callable, Generator, Iterable, Optional, Set, Union

from rxn.chemutils.reaction_equation import ReactionEquation

from rxn.reaction_preprocessing.annotations.annotation_criterion import (
    AnnotationCriterion,
)
from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    MoleculeAnnotation,
)


class MissingAnnotationDetector:
    """
    Find reactions with molecules that should be annotated, taking into account
    a set of already-annotated molecules.
    """

    def __init__(
        self,
        annotated_molecules: Set[str],
        requires_annotation_fn: Optional[Callable[[str], bool]] = None,
    ):
        """
        Args:
            annotated_molecules: set of already-annotated molecules.
            requires_annotation_fn: function with which to decide whether a molecule
                needs an annotation. Defaults to AnnotationCriterion().
        """
        self.annotated_molecules = annotated_molecules

        if requires_annotation_fn is None:
            requires_annotation_fn = AnnotationCriterion()
        self.requires_annotation_fn = requires_annotation_fn

    def molecule_needs_annotation(self, smiles: str) -> bool:
        """
        Whether a molecule needs annotation.

        Checks the overlap between the elements in the molecule and the extended
        transition metals, and then looks in the annotated molecules if necessary.
        """
        if not self.requires_annotation_fn(smiles):
            return False
        else:
            return smiles not in self.annotated_molecules

    def missing_in_reaction_equation(
        self, reaction_equation: ReactionEquation
    ) -> Generator[str, None, None]:
        """In a reaction equation, find the molecules requiring annotation."""
        for smiles in reaction_equation.iter_all_smiles():
            if self.molecule_needs_annotation(smiles):
                yield smiles

    def missing_in_reaction_equations(
        self, reaction_equations: Iterable[ReactionEquation]
    ) -> Generator[str, None, None]:
        """In multiple reaction equations, find the molecules requiring annotation."""
        for reaction_equation in reaction_equations:
            yield from self.missing_in_reaction_equation(reaction_equation)

    def missing_in_reaction_smiles(
        self,
        reaction_smiles: Union[Iterable[str], str],
        fragment_bond: Optional[str] = None,
    ) -> Generator[str, None, None]:
        """
        In one or multiple reaction SMILES, find the molecules requiring annotation.

        Args:
            reaction_smiles: One reaction SMILES (str), or multiple reaction SMILES.
            fragment_bond: fragment bond used in the reaction SMILES.
        """
        if isinstance(reaction_smiles, str):
            reaction_smiles = [reaction_smiles]

        reaction_equations = (
            ReactionEquation.from_string(reaction_smile, fragment_bond)
            for reaction_smile in reaction_smiles
        )
        return self.missing_in_reaction_equations(reaction_equations)

    @classmethod
    def from_molecule_annotations(
        cls,
        molecule_annotations: Iterable[MoleculeAnnotation],
        requires_annotation_fn: Optional[Callable[[str], bool]] = None,
    ) -> "MissingAnnotationDetector":
        """
        Create a MissingAnnotationDetector instance from existing molecule annotations.

        Args:
            molecule_annotations: existing molecule annotations.
            requires_annotation_fn: function with which to decide whether a molecule
                needs an annotation. Defaults to AnnotationCriterion().
        """
        original_smiles = {
            annotation.original_without_fragment_bond
            for annotation in molecule_annotations
        }
        # Also consider the updated SMILES, but only if they consist in exactly one molecule.
        updated_smiles = {
            annotation.updated_without_fragment_bond[0]
            for annotation in molecule_annotations
            if annotation.updated_smiles is not None
            and len(annotation.updated_without_fragment_bond) == 1
        }
        return cls(
            annotated_molecules=original_smiles | updated_smiles,
            requires_annotation_fn=requires_annotation_fn,
        )

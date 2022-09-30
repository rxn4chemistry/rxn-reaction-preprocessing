from typing import Iterable, List, Mapping, Optional, Union

from rxn.chemutils.reaction_equation import ReactionEquation

from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    AnnotationDecision,
    MoleculeAnnotation,
)


class MoleculeReplacer:
    """
    Class to replace SMILES strings by their improved alternatives (for instance
    from the catalyst annotations).
    """

    def __init__(self, replacements: Mapping[str, Union[str, List[str]]]):
        """
        Args:
            replacements: mapping between SMILES strings to replace, and what
                to replace them with (given as a string or list of strings).
                Use dots for fragment bonds!
        """

        # The mapped values given as strings are converted to lists
        self.replacements = {
            key: ([value] if isinstance(value, str) else value)
            for key, value in replacements.items()
        }

    def replace_molecule_smiles(self, smiles: str) -> List[str]:
        """
        Do the molecule replacements in a molecule SMILES.

        Returns:
            List of replaced molecules (will have length of one in most cases)
        """
        try:
            return self.replacements[smiles]
        except KeyError:
            # If not found: return the original SMILES (as a list!)
            return [smiles]

    def replace_in_reaction_smiles(
        self, smiles: str, fragment_bond: Optional[str] = None
    ) -> str:
        """
        Do the molecule replacements in a reaction SMILES.

        Args:
            smiles: reaction SMILES.
            fragment_bond: fragment bond used in the reaction SMILES.
        """
        reaction_equation = ReactionEquation.from_string(smiles, fragment_bond)
        return self.replace_in_reaction_equation(reaction_equation).to_string(
            fragment_bond
        )

    def replace_in_reaction_equation(
        self, reaction_equation: ReactionEquation
    ) -> ReactionEquation:
        """
        Do the molecule replacements in a ReactionEquation instance.
        """
        groups = (self._replace_in_molecule_list(group) for group in reaction_equation)
        return ReactionEquation(*groups)

    def _replace_in_molecule_list(self, molecules: List[str]) -> List[str]:
        """
        Do the replacements in a list of SMILES, potentially leading to a larger list.
        """
        return [
            replaced_molecule
            for molecule in molecules
            for replaced_molecule in self.replace_molecule_smiles(molecule)
        ]

    @classmethod
    def from_molecule_annotations(
        cls, molecule_annotations: Iterable[MoleculeAnnotation]
    ) -> "MoleculeReplacer":
        """Instantiate from a list of molecule annotations."""

        def criterion(annotation: MoleculeAnnotation) -> bool:
            """Whether to add a replacement rule for a given annotation."""
            is_accepted = annotation.decision == AnnotationDecision.ACCEPT
            has_updated_smiles = annotation.updated_smiles is not None
            return is_accepted and has_updated_smiles

        replacements = {
            annotation.original_without_fragment_bond: annotation.updated_without_fragment_bond
            for annotation in molecule_annotations
            if criterion(annotation)
        }
        return cls(replacements)

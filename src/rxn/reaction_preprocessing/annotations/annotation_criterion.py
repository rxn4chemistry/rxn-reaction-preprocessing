import itertools
from typing import Generator, Iterable, Optional, Set

from rdkit.Chem import GetPeriodicTable
from rxn.chemutils.miscellaneous import atom_type_counter


class AnnotationCriterion:
    """
    Determine what molecules need an annotation, independently of the molecules
    that have been annotated so far.

    To do so, this class relies on whether the molecules contain (extended)
    transition metals or not.
    """

    def __init__(
        self,
        additional_elements_to_consider: Optional[Iterable[str]] = None,
        elements_not_to_consider: Optional[Iterable[str]] = None,
    ):
        """
        Args:
            additional_elements_to_consider: elements for which to require an
                annotation, in addition to the extended transition metals.
            elements_not_to_consider: elements for which not to require an
                annotation, even if they are extended transition metals.
        """
        self.elements_requiring_annotation = set(
            AnnotationCriterion.extended_transition_metals()
        )

        if additional_elements_to_consider is not None:
            self.elements_requiring_annotation.update(additional_elements_to_consider)
        if elements_not_to_consider is not None:
            self.elements_requiring_annotation.difference_update(
                elements_not_to_consider
            )

    def __call__(self, smiles: str) -> bool:
        """
        Function making the object callable, falling back to "requires_annotation".

        Args:
            smiles: molecule SMILES. Use dots for fragment bonds!
        """
        return self.requires_annotation(smiles)

    def requires_annotation(self, smiles: str) -> bool:
        """
        Whether a given SMILES string requires an annotation.

        Args:
            smiles: molecule SMILES. Use dots for fragment bonds!
        """
        return bool(
            AnnotationCriterion.elements_in_smiles(smiles)
            & self.elements_requiring_annotation
        )

    @staticmethod
    def elements_in_smiles(smiles: str) -> Set[str]:
        return set(atom_type_counter(smiles).keys())

    @staticmethod
    def extended_transition_metals() -> Generator[str, None, None]:
        """
        Atomic symbols for the extended transition metals.
        """
        periodic_table = GetPeriodicTable()
        return (
            periodic_table.GetElementSymbol(atomic_number)
            for atomic_number in AnnotationCriterion.extended_transition_metal_numbers()
        )

    @staticmethod
    def extended_transition_metal_numbers() -> Generator[int, None, None]:
        """
        Atomic numbers for the extended transition metals.
        """
        yield from itertools.chain(
            # Al
            (13,),
            # first-row transition metals + Ga
            range(21, 32),
            # second-row transition metals + In
            range(39, 50),
            # lanthanides + third-row transition metals + Tl + Pb + Bi + Po
            range(57, 85),
            # actinides + fourth-row transition metals
            range(89, 113),
        )

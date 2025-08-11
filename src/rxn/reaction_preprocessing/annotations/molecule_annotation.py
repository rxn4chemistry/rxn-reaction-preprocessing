import json
from enum import auto
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Union

import attr
from rxn.chemutils.multicomponent_smiles import multicomponent_smiles_to_list
from rxn.utilities.types import RxnEnum


class AnnotationDecision(RxnEnum):
    ACCEPT = auto()
    REJECT = auto()


@attr.s(auto_attribs=True, init=False)
class MoleculeAnnotation:
    """
    Specifies a molecule annotation, i.e. a SMILES string that may have an
    updated SMILES, whether to keep it or not, etc.
    """

    original_smiles: str
    updated_smiles: Optional[str]
    decision: AnnotationDecision
    categories: List[str]
    extra_info: Dict[str, Any]

    def __init__(
        self,
        original_smiles: str,
        updated_smiles: Optional[str],
        decision: str,
        categories: List[str],
        **extra_info: Any
    ):
        """
        Args:
            original_smiles: original SMILES that is potentially present in a
                data set. Fragment bonds are indicated by a tilde '~'.
            updated_smiles: if specified, SMILES with which to replace
                original_smiles. Also uses '~' for fragment bonds, and dots
                '.' may be used to separate compounds from the solvent in which
                they are solved.
            decision: "accept" or "reject".
            categories: categories to which the annotation belongs to.
            **extra_info: additional information not covered by the other
                variables.
        """

        decision_enum = AnnotationDecision.from_string(decision)
        self.__attrs_init__(
            original_smiles=original_smiles,
            updated_smiles=updated_smiles,
            decision=decision_enum,
            categories=categories,
            extra_info=extra_info,
        )

    @property
    def original_without_fragment_bond(self) -> str:
        """Get the original SMILES with dots instead of tildes to delimit fragments."""
        return self.original_smiles.replace("~", ".")

    @property
    def updated_without_fragment_bond(self) -> List[str]:
        """
        Get the updated SMILES with dots instead of tildes to delimit fragments.

        Since dots may be used to delimit solvents from compounds, a list must be returned.
        """
        if self.updated_smiles is None:
            raise RuntimeError("No updated SMILES!")
        return multicomponent_smiles_to_list(self.updated_smiles, "~")


def load_annotations(json_file: Union[Path, str]) -> List[MoleculeAnnotation]:
    """
    Load the molecule annotations from a JSON file.

    Args:
        json_file: path to the JSON file containing the annotations.

    Returns:
        List of annotations.
    """

    with open(json_file, "rt") as f:
        json_content = json.load(f)

    return [MoleculeAnnotation(**block) for block in json_content]


def load_annotations_multiple(
    json_files: Iterable[Union[Path, str]],
) -> List[MoleculeAnnotation]:
    """
    Load the molecule annotations from multiple JSON files.

    Args:
        json_files: paths to the JSON file containing the annotations.

    Returns:
        List of annotations.
    """
    annotations = []
    for json_file in json_files:
        annotations.extend(load_annotations(json_file))
    return annotations

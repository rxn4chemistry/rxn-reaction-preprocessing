import copy
import logging
from typing import Any, Dict

from rxn.chemutils.reaction_smiles import parse_any_reaction_smiles
from rxn.chemutils.utils import remove_atom_mapping

from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    load_annotations_multiple,
)
from rxn.reaction_preprocessing.config import PreprocessConfig, StandardizeConfig
from rxn.reaction_preprocessing.mixed_reaction_filter import MixedReactionFilter
from rxn.reaction_preprocessing.molecule_standardizer import MoleculeStandardizer
from rxn.reaction_preprocessing.reaction_standardizer import ReactionStandardizer

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class PistachioRecordStandardizer:
    """
    Class to standardize the reactions in Pistachio records.

    It combines some aspects of the STANDARDIZE and PREPROCESS steps.

    The main goal is to standardize the reaction SMILES and filter invalid
    reaction records. The individual SMILES strings for the components are
    also standardized, but the error is not propagated if this fails.

    Duplicates are not removed here, as the criterion to do so is unclear, and
    also because this class does not store previously seen reactions.

    The initialization from StandardizeConfig and PreprocessConfig is not optimal,
    as it (currently) contains more information than is needed, such as the location
    of the CSV files to load from and write to.
    """

    def __init__(
        self, cfg_standardize: StandardizeConfig, cfg_preprocess: PreprocessConfig
    ):
        self.fragment_bond = cfg_standardize.fragment_bond.value
        self.molecule_standardizer = MoleculeStandardizer(
            annotations=load_annotations_multiple(
                cfg_standardize.annotation_file_paths
            ),
            discard_missing_annotations=cfg_standardize.discard_unannotated_metals,
            canonicalize=True,
        )
        self.reaction_standardizer = ReactionStandardizer()

        self.reaction_filter = MixedReactionFilter(
            max_reactants=cfg_preprocess.max_reactants,
            max_agents=cfg_preprocess.max_agents,
            max_products=cfg_preprocess.max_products,
            min_reactants=cfg_preprocess.min_reactants,
            min_agents=cfg_preprocess.min_agents,
            min_products=cfg_preprocess.min_products,
            max_reactants_tokens=cfg_preprocess.max_reactants_tokens,
            max_agents_tokens=cfg_preprocess.max_agents_tokens,
            max_products_tokens=cfg_preprocess.max_products_tokens,
            max_absolute_formal_charge=cfg_preprocess.max_absolute_formal_charge,
        )

        # in earlier versions, this was "reactionSmiles"
        self.reaction_smiles_field = "smiles"

    def standardize(self, reaction_record: Dict[str, Any]) -> Dict[str, Any]:
        """
        Standardize a Pistachio reaction record.

        Args:
            reaction_record:

        Returns:

        """
        # Do not do changes in-place
        reaction_record = copy.deepcopy(reaction_record)

        reaction_smiles = reaction_record["data"][self.reaction_smiles_field]
        reaction_smiles = remove_atom_mapping(reaction_smiles)
        reaction = parse_any_reaction_smiles(reaction_smiles)
        reaction = self.molecule_standardizer.standardize_in_equation(reaction)
        reaction = self.reaction_standardizer.standardize(reaction)

        # This will raise if the filters do not pass
        self.reaction_filter.validate(reaction)

        for component_dict in reaction_record.get("components", []):
            self._try_standardize_component(component_dict)

        for action_dict in reaction_record.get("actions", []):
            for component_dict in action_dict.get("components", []):
                self._try_standardize_component(component_dict)

        reaction_record["data"][self.reaction_smiles_field] = reaction.to_string(
            self.fragment_bond
        )

        return reaction_record

    def _try_standardize_component(self, component_dict: Dict[str, Any]) -> None:
        """Standardize one component in the Pistachio record. Does not raise
        any exception if this fails."""
        try:
            smiles = component_dict["smiles"]
            standardized_smiles = self.molecule_standardizer.standardize(smiles)

            # IMPORTANT: as we need a string, not a list, we must convert
            # the list to a string again, i.e. use a fragment bond even if
            # the molecule standardizer provides multiple separate molecules.
            component_dict["smiles"] = ".".join(standardized_smiles)
        except Exception:
            # if there was any error (KeyError, standardization error, etc.), we
            # do nothing
            return

import copy
from typing import List

from rxn_chemutils.extended_reaction_smiles import parse_extended_reaction_smiles
from rxn_chemutils.reaction_equation import (
    merge_reactants_and_agents,
    sort_compounds,
    canonicalize_compounds,
)


class ReactionSmilesExtractor:
    """
    !!! DEPRECATED! Use ReactionSmilesConverter and ReactionEquation instead. !!!
    Enables extraction of various information from the reaction SMILES.
    Determines three variables:
    self.reactants: RDKit-canonicalized SMILES strings for the reactants (left of the reaction arrow)
    self.agents: RDKit-canonicalized SMILES strings for the products (above the reaction arrow)
    self.products RDKit-canonicalized SMILES strings for the products (right of the reaction arrow)
    """

    def __init__(self, reaction_smiles: str):
        reaction = parse_extended_reaction_smiles(reaction_smiles)
        self.reaction = canonicalize_compounds(reaction)

    @property
    def reactants(self) -> List[str]:
        return self.reaction.reactants

    @property
    def agents(self) -> List[str]:
        return self.reaction.agents

    @property
    def products(self) -> List[str]:
        return self.reaction.products

    def get_fragmented_reaction_smiles(
        self,
        mixed: bool = False,
        order_molecules: bool = False,
        fragment_bond: str = None,
    ) -> str:
        """Get a reaction smiles string, where the fragments belonging
        together are grouped.
        Keyword Arguments:
            mixed {bool} -- [make no disctinction between reactant and agents and mix them] (default: {False})
            order_molecules {bool} -- [order the molecules alphabetically (RDKit default) to have a canonical string] (default: {False})
            fragment_bond {str} -- [keep the fragment information by replacing the '.' within groups with] (default: {None})
        Returns:
            str -- [reaction smiles with fragment information]
        """

        # NB: the compounds were already canonicalized in the constructor, not necessary to do it here
        reaction = copy.deepcopy(self.reaction)

        if mixed:
            reaction = merge_reactants_and_agents(reaction)

        if order_molecules:
            reaction = sort_compounds(reaction)

        return reaction.to_string(fragment_bond=fragment_bond)
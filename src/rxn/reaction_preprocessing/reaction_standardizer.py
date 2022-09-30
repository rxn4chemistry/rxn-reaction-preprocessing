from rxn.chemutils.reaction_equation import (
    ReactionEquation,
    merge_reactants_and_agents,
    remove_duplicate_compounds,
    remove_precursors_from_products,
    sort_compounds,
)


class ReactionStandardizer:
    """
    Standardization on the reaction level.

    I.e. this class does not care about canonical SMILES strings, but rather
    about removing duplicate reagents, merging reactants and reagents, etc.
    To be effective, this class relies on the molecules in the reaction
    SMILES to be canonical already.

    NB: at the moment (Oct 2021), this class is called not in the standardization
    step, but in the preprocess step. This may change in the future.
    """

    def __call__(self, reaction: ReactionEquation) -> ReactionEquation:
        """See doc for standardize()."""
        return self.standardize(reaction)

    def standardize(self, reaction: ReactionEquation) -> ReactionEquation:
        """
        Standardize a reaction.

        This function will make a copy and not modify the initial instance.

        Args:
            reaction: the reaction to standardize.

        Returns:
            A standardized reaction.
        """

        # NB: this produces a copy, so no explicit deep copy necessary
        reaction = merge_reactants_and_agents(reaction)
        reaction = remove_precursors_from_products(reaction)
        reaction = remove_duplicate_compounds(reaction)
        reaction = sort_compounds(reaction)

        return reaction

# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
"""A class encapsulating filtering functionality for chemical reactions"""
import itertools
from functools import partial
from typing import Callable, Generator, Iterable, List, Tuple, Union

from rxn.chemutils.exceptions import InvalidSmiles
from rxn.chemutils.reaction_equation import ReactionEquation
from rxn.chemutils.tokenization import to_tokens

from .utils import MolEquation, get_atoms_for_mols, get_formal_charge_for_mols

_POLYMER_HEAD_AND_TAIL_PLACEHOLDER_ATOMS = {"Kr", "Rn", "Xe"}
_ATOM_TYPES_ALLOWED_IN_PRODUCT = _POLYMER_HEAD_AND_TAIL_PLACEHOLDER_ATOMS | {"H"}

SmilesBasedCheck = Callable[[ReactionEquation], bool]
MolBasedCheck = Callable[[MolEquation], bool]


class ReactionFilterError(ValueError):
    """Exception raised when calling validate() on reactions not passing one
    or several filters."""

    def __init__(self, reaction: ReactionEquation, reasons: Iterable[str]):
        # Store information just in case
        self.reaction = reaction
        self.reasons = list(reasons)

        super().__init__(
            f'Reaction "{self.reaction.to_string("~")}" did not pass the '
            f'filters: {"; ".join(self.reasons)}'
        )


class MixedReactionFilter:
    def __init__(
        self,
        max_reactants: int = 10,
        max_agents: int = 0,
        max_products: int = 1,
        min_reactants: int = 2,
        min_agents: int = 0,
        min_products: int = 1,
        max_reactants_tokens: int = 300,
        max_agents_tokens: int = 0,
        max_products_tokens: int = 200,
        max_absolute_formal_charge: int = 2,
    ):
        """Creates a new instance of the type MixedReactionFilter.

        Args:
            max_reactants: The maximum number of reactant molecules.
            max_agents: The maximum number of agent molcules.
            max_products: The maximum number of product molecules.
            min_reactants: The minimum number of reactant molecules.
            min_agents: The minium number of agent molecules.
            min_products: The minimum number of product molecules.
            max_reactants_tokens: The maximum number of precursor tokens.
            max_agents_tokens: The maximum number of agent tokens.
            max_products_tokens: The maximum number of product tokens.
            max_absolute_formal_charge: The maximum formal charge (for
                reactants, agents, or products).
        """
        self.max_reactants = max_reactants
        self.max_agents = max_agents
        self.max_products = max_products
        self.min_reactants = min_reactants
        self.min_agents = min_agents
        self.min_products = min_products
        self.max_reactants_tokens = max_reactants_tokens
        self.max_agents_tokens = max_agents_tokens
        self.max_products_tokens = max_products_tokens
        self.max_absolute_formal_charge = max_absolute_formal_charge

        self.smiles_based_checks: List[Tuple[SmilesBasedCheck, str]] = [
            (self.max_reactants_exceeded, "max_reactants_exceeded"),
            (self.max_agents_exceeded, "max_agents_exceeded"),
            (self.max_products_exceeded, "max_products_exceeded"),
            (self.min_reactants_subceeded, "min_reactants_subceeded"),
            (self.min_agents_subceeded, "min_agents_subceeded"),
            (self.min_products_subceeded, "min_products_subceeded"),
            (self.products_subset_of_reactants, "products_subset_of_reactants"),
            (self.max_reactant_tokens_exceeded, "max_reactant_tokens_exceeded"),
            (self.max_agent_tokens_exceeded, "max_agent_tokens_exceeded"),
            (self.max_product_tokens_exceeded, "max_product_tokens_exceeded"),
        ]
        self.mol_based_checks: List[Tuple[MolBasedCheck, str]] = [
            (self.products_single_atoms, "products_single_atoms"),
            (self.formal_charge_exceeded, "formal_charge_exceeded"),
            (self.invalid_atom_type, "invalid_atom_type"),
            (self.different_atom_types, "different_atom_types"),
        ]

    def validate(self, reaction: ReactionEquation) -> None:
        """
        Make sure that the given reaction is valid; if not, an exception will
        be raised.

        Raises:
            ReactionFilterError: if the reaction does not pass the filters.

        Args:
            reaction: reaction to validate.
        """
        valid, reasons = self.validate_reasons(reaction)

        if not valid:
            raise ReactionFilterError(reaction, reasons)

    def is_valid(self, reaction: ReactionEquation) -> bool:
        """
        Whether a reaction is valid based on the rules set on the instance of this
        MixedReactionFilter class.

        Args:
            reaction: The reaction to validate.

        Returns:
            bool: Whether or not the reaction is valid according to the rules
            set on the instance of this MixedReactionFilter class.
        """

        def callbacks() -> Generator[Callable[[], bool], None, None]:
            """Generator function for providing the checks to make as
            callable objects.

            Formulating it as a generator makes it efficient; for instance,
            the mol_equation object will not be generated if any of the
            SMILES-based checks fails.
            """
            for smiles_based_fn, _ in self.smiles_based_checks:
                yield partial(smiles_based_fn, reaction)

            try:
                mol_equation = MolEquation.from_reaction_equation(reaction)
            except InvalidSmiles:
                # If there is an invalid SMILES, we yield a final callback that
                # will then return `True` (meaning: invalid)
                yield lambda: True
                return

            for mol_based_fn, _ in self.mol_based_checks:
                yield partial(mol_based_fn, mol_equation)

        return not any(callback() for callback in callbacks())

    def validate_reasons(self, reaction: ReactionEquation) -> Tuple[bool, List[str]]:
        reasons = []

        for smiles_based_fn, error_message in self.smiles_based_checks:
            if smiles_based_fn(reaction):
                reasons.append(error_message)

        try:
            mol_equation = MolEquation.from_reaction_equation(reaction)
        except InvalidSmiles:
            reasons.append("rdkit_molfromsmiles_failed")
        else:
            for mol_based_fn, error_message in self.mol_based_checks:
                if mol_based_fn(mol_equation):
                    reasons.append(error_message)

        valid = len(reasons) == 0

        return valid, reasons

    def max_reactants_exceeded(self, reaction: ReactionEquation) -> bool:
        """Checks whether the number of reactants exceeds the maximum.

        Args:
            reaction: The reaction to test.

        Returns:
            bool: Whether the number of reactants exceeds the maximum.
        """

        return len(reaction.reactants) > self.max_reactants

    def max_agents_exceeded(self, reaction: ReactionEquation) -> bool:
        """Checks whether the number of agents exceeds the maximum.

        Args:
            reaction: The reaction to test.

        Returns:
            bool: Whether the number of agents exceeds the maximum.
        """

        return len(reaction.agents) > self.max_agents

    def max_products_exceeded(self, reaction: ReactionEquation) -> bool:
        """Checks whether the number of products exceeds the maximum.

        Args:
            reaction: The reaction to test.

        Returns:
            bool: Whether the number of products exceeds the maximum.
        """

        return len(reaction.products) > self.max_products

    def min_reactants_subceeded(self, reaction: ReactionEquation) -> bool:
        """Checks whether the number of reactants exceeds the maximum.

        Args:
            reaction: The reaction to test.

        Returns:
            bool: Whether the number of reactants exceeds the maximum.
        """

        return len(reaction.reactants) < self.min_reactants

    def min_agents_subceeded(self, reaction: ReactionEquation) -> bool:
        """Checks whether the number of agents exceeds the maximum.

        Args:
            reaction: The reaction to test.

        Returns:
            bool: Whether the number of agents exceeds the maximum.
        """

        return len(reaction.agents) < self.min_agents

    def min_products_subceeded(self, reaction: ReactionEquation) -> bool:
        """Checks whether the number of products exceeds the maximum.

        Args:
            reaction: The reaction to test.

        Returns:
            bool: Whether the number of products exceeds the maximum.
        """

        return len(reaction.products) < self.min_products

    def products_subset_of_reactants(self, reaction: ReactionEquation) -> bool:
        """Checks whether the set of products is a subset of the set of reactants.

        Args:
            reaction: The reaction to test.

        Returns:
            bool: Whether the set of products is a subset of the set of reactants.
        """
        products = set(reaction.products)
        reactants = set(reaction.reactants)

        return len(products) > 0 and products.issubset(reactants)

    def products_single_atoms(
        self, reaction: Union[MolEquation, ReactionEquation]
    ) -> bool:
        """Checks whether the products solely consist of single atoms.

        Args:
            reaction: The reaction to test.

        Returns:
            bool: Whether the products solely consist of single atoms.
        """
        if isinstance(reaction, ReactionEquation):
            reaction = MolEquation.from_reaction_equation(reaction)

        return len(reaction.products) > 0 and all(
            [product.GetNumAtoms() == 1 for product in reaction.products]
        )

    def max_reactant_tokens_exceeded(self, reaction: ReactionEquation) -> bool:
        """Check whether the number of reactant tokens exceeds the maximum.

        Args:
            reaction: The reaction to test.

        Returns:
            [type]: Whether the number of reactant tokens exceeds the maximum.
        """
        return self._group_tokens_exceeded(
            reaction.reactants, self.max_reactants_tokens
        )

    def max_agent_tokens_exceeded(self, reaction: ReactionEquation) -> bool:
        """Check whether the number of agent tokens exceeds the maximum.

        Args:
            reaction: The reaction to test.

        Returns:
            [type]: Whether the number of agent tokens exceeds the maximum.
        """
        return self._group_tokens_exceeded(reaction.agents, self.max_agents_tokens)

    def max_product_tokens_exceeded(self, reaction: ReactionEquation) -> bool:
        """Check whether the number of product tokens exceeds the maximum.

        Args:
            reaction: The reaction to test.

        Returns:
            [type]: Whether the number of product tokens exceeds the maximum.
        """
        return self._group_tokens_exceeded(reaction.products, self.max_products_tokens)

    def _group_tokens_exceeded(self, smiles_list: List[str], threshold: int) -> bool:
        """Check whether the number of SMILES tokens in a group exceeds the maximum.

        Returns:
            [type]: Whether the number of product tokens exceeds the maximum.
        """

        smiles = ".".join(
            smiles_list
        )  # NB: we use '.' here, but '~' would be the same.
        return len(to_tokens(smiles)) > threshold

    def formal_charge_exceeded(
        self, reaction: Union[MolEquation, ReactionEquation]
    ) -> bool:
        """Check whether the absolute formal charge of the reactants, agents,
        or products exceeds a maximum.

        Args:
            reaction: The reaction to test.

        Returns:
            bool: Whether the absolute formal charge of the reactants, agents,
            or products exceeds a maximum.
        """
        if isinstance(reaction, ReactionEquation):
            reaction = MolEquation.from_reaction_equation(reaction)

        return (
            abs(get_formal_charge_for_mols(reaction.reactants))
            > self.max_absolute_formal_charge
            or abs(get_formal_charge_for_mols(reaction.agents))
            > self.max_absolute_formal_charge
            or abs(get_formal_charge_for_mols(reaction.products))
            > self.max_absolute_formal_charge
        )

    def invalid_atom_type(self, reaction: Union[MolEquation, ReactionEquation]) -> bool:
        """
        Check whether the reaction contains atoms with invalid atom types such as the asterisk "*".

        Args:
            reaction: The reaction to test.

        Returns:
            bool: Whether the reaction contains invalid atom types.
        """
        if isinstance(reaction, ReactionEquation):
            reaction = MolEquation.from_reaction_equation(reaction)

        # So far, the only invalid atom type is "*"; this function can be
        # reformulated to account for additional ones if some appear later on.
        mols = itertools.chain(reaction.reactants, reaction.agents, reaction.products)
        return "*" in get_atoms_for_mols(mols)

    def different_atom_types(
        self, reaction: Union[MolEquation, ReactionEquation]
    ) -> bool:
        """Check whether the products contain atom types not found in the agents or reactants.

        It handles the presence of placeholders for polymer head and tail representations.

        Args:
            reaction: The reaction to test.

        Returns:
            bool: Whether the products contain atom types not found in the agents or reactants.
        """
        if isinstance(reaction, ReactionEquation):
            reaction = MolEquation.from_reaction_equation(reaction)

        products_atoms = get_atoms_for_mols(reaction.products)
        # ignore H atom (because usually implicit) and atoms used in polymer representations
        products_atoms -= _ATOM_TYPES_ALLOWED_IN_PRODUCT
        agents_atoms = get_atoms_for_mols(reaction.agents)
        reactants_atoms = get_atoms_for_mols(reaction.reactants)

        return len(products_atoms - (reactants_atoms | agents_atoms)) != 0

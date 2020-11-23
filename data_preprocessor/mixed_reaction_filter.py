""" A class encapsulating filtering functionality for chemical reactions """

from typing import List, Tuple
from rdkit.Chem import AllChem as rdk
from .reaction import Reaction
from .smiles_tokenizer import SmilesTokenizer


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
        # filter_single_atom_products: bool = True,
        # filter_purifications: bool = True,
        # filter_alchemy: bool = True,
    ):
        """Creates a new instance of the type MixedReactionFilter.

        Args:
            max_reactants (int, optional): The maximum number of reactant molecules. Defaults to 10.
            max_agents (int, optional): The maximum number of agent molcules. Defaults to 0.
            max_products (int, optional): The maximum number of product molecules. Defaults to 1.
            min_reactants (int, optional): The minimum number of reactant molecules. Defaults to 2.
            min_agents (int, optional): The minium number of agent molecules. Defaults to 0.
            min_products (int, optional): The minimum number of product molecules. Defaults to 1.
            max_reactants_tokens (int, optional): The maximum number of precursor tokens. Defaults to 300.
            max_agents_tokens (int, optional): The maximum number of agent tokens. Defaults to 0.
            max_products_tokens (int, optional): The maximum number of product tokens. Defaults to 200.
            max_absolute_formal_charge (int, optional): The maximum formal charge (for reactants, agents, or products). Defaults to 2.
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

        self.tokenizer = SmilesTokenizer()

    def validate(self, reaction: Reaction) -> bool:
        """Validate a reaction using the rules set on the instance of this MixedReactionFilter class.

        Args:
            reaction (Reaction): The reaction to validate.

        Returns:
            bool: Whether or not the reaction is valid according to the rules set on the instance of this MixedReactionFilter class.
        """

        return not any(
            [
                self.max_reactants_exceeded(reaction),
                self.max_agents_exceeded(reaction),
                self.max_products_exceeded(reaction),
                self.min_reactants_subceeded(reaction),
                self.min_agents_subceeded(reaction),
                self.min_products_subceeded(reaction),
                self.products_subset_of_reactants(reaction),
                self.products_single_atoms(reaction),
                self.max_reactant_tokens_exceeded(reaction),
                self.max_agent_tokens_exceeded(reaction),
                self.max_product_tokens_exceeded(reaction),
                self.formal_charge_exceeded(reaction),
                self.different_atom_types(reaction),
            ]
        )

    def validate_reasons(self, reaction: Reaction) -> Tuple[bool, List[str]]:
        valid = True
        reasons = []

        if self.max_reactants_exceeded(reaction):
            valid = False
            reasons.append("max_reactants_exceeded")
        if self.max_agents_exceeded(reaction):
            valid = False
            reasons.append("max_agents_exceeded")
        if self.max_products_exceeded(reaction):
            valid = False
            reasons.append("max_products_exceeded")
        if self.min_reactants_subceeded(reaction):
            valid = False
            reasons.append("min_reactants_subceeded")
        if self.min_agents_subceeded(reaction):
            valid = False
            reasons.append("min_agents_subceeded")
        if self.min_products_subceeded(reaction):
            valid = False
            reasons.append("min_products_subceeded")
        if self.products_subset_of_reactants(reaction):
            valid = False
            reasons.append("products_subset_of_reactants")
        if self.products_single_atoms(reaction):
            valid = False
            reasons.append("products_single_atoms")
        if self.max_reactant_tokens_exceeded(reaction):
            valid = False
            reasons.append("max_reactant_tokens_exceeded")
        if self.max_agent_tokens_exceeded(reaction):
            valid = False
            reasons.append("max_agent_tokens_exceeded")
        if self.max_product_tokens_exceeded(reaction):
            valid = False
            reasons.append("max_product_tokens_exceeded")
        if self.formal_charge_exceeded(reaction):
            valid = False
            reasons.append("formal_charge_exceeded")
        if self.different_atom_types(reaction):
            valid = False
            reasons.append("different_atom_types")

        return (valid, reasons)

    def max_reactants_exceeded(self, reaction: Reaction) -> bool:
        """Checks whether the number of reactants exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the number of reactants exceeds the maximum.
        """

        return len(reaction.reactants) > self.max_reactants

    def max_agents_exceeded(self, reaction: Reaction) -> bool:
        """Checks whether the number of agents exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the number of agents exceeds the maximum.
        """

        return len(reaction.agents) > self.max_agents

    def max_products_exceeded(self, reaction: Reaction) -> bool:
        """Checks whether the number of products exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the number of products exceeds the maximum.
        """

        return len(reaction.products) > self.max_products

    def min_reactants_subceeded(self, reaction: Reaction) -> bool:
        """Checks whether the number of reactants exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the number of reactants exceeds the maximum.
        """

        return len(reaction.reactants) < self.min_reactants

    def min_agents_subceeded(self, reaction: Reaction) -> bool:
        """Checks whether the number of agents exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the number of agents exceeds the maximum.
        """

        return len(reaction.agents) < self.min_agents

    def min_products_subceeded(self, reaction: Reaction) -> bool:
        """Checks whether the number of products exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the number of products exceeds the maximum.
        """

        return len(reaction.products) < self.min_products

    def products_subset_of_reactants(self, reaction: Reaction) -> bool:
        """Checks whether the set of products is a subset of the set of reactants.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the set of products is a subset of the set of reactants.
        """

        return set(reaction.get_products_as_smiles()).issubset(
            set(reaction.get_reactants_as_smiles())
        )

    def products_single_atoms(self, reaction: Reaction) -> bool:
        """Checks whether the products solely consist of single atoms.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the products solely consist of single atoms.
        """
        if not all([product.GetNumAtoms() == 1 for product in reaction.products]):
            return False
        return True

    def max_reactant_tokens_exceeded(self, reaction: Reaction) -> bool:
        """Check whether the number of reactant tokens exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            [type]: Whether the number of reactant tokens exceeds the maximum.
        """

        smiles = reaction.get_reactants_as_smiles()

        if len(smiles) > self.max_reactants_tokens:
            return True

        return (
            len(self.tokenizer.tokenize(".".join(smiles)).split(" "))
            > self.max_reactants_tokens
        )

    def max_agent_tokens_exceeded(self, reaction: Reaction) -> bool:
        """Check whether the number of agent tokens exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            [type]: Whether the number of agent tokens exceeds the maximum.
        """

        smiles = reaction.get_agents_as_smiles()

        # In this and the other token counting methods, this is done to capture
        # the empty set, which will not be empty when joined
        if len(smiles) <= self.max_agents_tokens:
            return False

        return (
            len(self.tokenizer.tokenize(".".join(smiles)).split(" "))
            > self.max_agents_tokens
        )

    def max_product_tokens_exceeded(self, reaction: Reaction) -> bool:
        """Check whether the number of product tokens exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            [type]: Whether the number of product tokens exceeds the maximum.
        """

        smiles = reaction.get_products_as_smiles()

        if len(smiles) > self.max_products_tokens:
            return True

        return (
            len(self.tokenizer.tokenize(".".join(smiles)).split(" "))
            > self.max_products_tokens
        )

    def formal_charge_exceeded(self, reaction: Reaction) -> bool:
        """Check whether the absolute formal charge of the reactants, agents, or products exceeds a maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the absolute formal charge of the reactants, agents, or products exceeds a maximum.
        """

        # Fragment bonds should, if the user choses to, be dealt with
        # before filtering reactions
        return (
            abs(
                rdk.GetFormalCharge(
                    rdk.MolFromSmiles(".".join(reaction.get_reactants_as_smiles()))
                )
            )
            > self.max_absolute_formal_charge
            or abs(
                rdk.GetFormalCharge(
                    rdk.MolFromSmiles(".".join(reaction.get_agents_as_smiles()))
                )
            )
            > self.max_absolute_formal_charge
            or abs(
                rdk.GetFormalCharge(
                    rdk.MolFromSmiles(".".join(reaction.get_products_as_smiles()))
                )
            )
            > self.max_absolute_formal_charge
        )

    def different_atom_types(self, reaction: Reaction) -> bool:
        """Check whether the products contain atom types not found in the agents or reactants.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the products contain atom types not found in the agents or reactants.
        """

        reactants = rdk.MolFromSmiles(".".join(reaction.get_reactants_as_smiles()))
        agents = rdk.MolFromSmiles(".".join(reaction.get_agents_as_smiles()))
        products = rdk.MolFromSmiles(".".join(reaction.get_products_as_smiles()))

        products_atoms = set([atom.GetSymbol() for atom in products.GetAtoms()])
        agents_atoms = set([atom.GetSymbol() for atom in agents.GetAtoms()])
        reactants_atoms = set([atom.GetSymbol() for atom in reactants.GetAtoms()])

        return len(products_atoms - (reactants_atoms | agents_atoms)) != 0

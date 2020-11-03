""" A class encapsulating filtering functionality for chemical reactions """

from typing import List
from rdkit.Chem import AllChem as rdk
from .reaction import Reaction
from .smiles_tokenizer import SmilesTokenizer


class MixedReactionFilter:
    def __init__(
        self,
        max_number_of_reactants: int = 10,
        max_number_of_products: int = 1,
        max_number_of_agents: int = 0,
        max_number_of_precursor_tokens: int = 300,
        max_number_of_agent_tokens: int = 0,
        max_number_of_product_tokens=200,
        max_absolute_formal_charge: int = 2,
        # filter_single_atom_products: bool = True,
        # filter_purifications: bool = True,
        # filter_alchemy: bool = True,
    ):
        self.max_number_of_products = max_number_of_products
        self.max_number_of_reactants = max_number_of_reactants
        self.max_number_of_agents = max_number_of_agents
        self.max_number_of_product_tokens = max_number_of_product_tokens
        self.max_number_of_agent_tokens = max_number_of_agent_tokens
        self.max_number_of_precursor_tokens = max_number_of_precursor_tokens
        self.max_absolute_formal_charge = max_absolute_formal_charge

        self.tokenizer = SmilesTokenizer()

    def validate(self, reaction: Reaction) -> bool:
        return any(
            [
                self.max_reactants_exceeded(reaction),
                self.max_agents_exceeded(reaction),
                self.max_products_exceeded(reaction),
                self.products_subset_of_reactants(reaction),
                self.products_single_atoms(reaction),
                self.max_n_reactant_tokens_exceeded(reaction),
                self.max_n_agent_tokens_exceeded(reaction),
                self.max_n_product_tokens_exceeded(reaction),
                self.formal_charge_exceeded(reaction),
                self.different_atom_types(reaction),
            ]
        )

    def max_reactants_exceeded(self, reaction: Reaction) -> bool:
        """Checks whether the number of reactants exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the number of reactants exceeds the maximum.
        """

        return len(reaction.reactants) > self.max_number_of_reactants

    def max_agents_exceeded(self, reaction: Reaction) -> bool:
        """Checks whether the number of agents exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the number of agents exceeds the maximum.
        """

        return len(reaction.agents) > self.max_number_of_agents

    def max_products_exceeded(self, reaction: Reaction) -> bool:
        """Checks whether the number of products exceeds the maximum.

        Args:
            reaction (Reaction): The reaction to test.

        Returns:
            bool: Whether the number of products exceeds the maximum.
        """

        return len(reaction.products) > self.max_number_of_products

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

        for product in reaction.products:
            if product.GetNumAtoms() > 1:
                return False
        return True

    def max_n_reactant_tokens_exceeded(self, reaction):
        return (
            len(
                self.tokenizer.tokenize(
                    ".".join(reaction.get_reactants_as_smiles())
                ).split(" ")
            )
            > self.max_number_of_precursor_tokens
        )

    def max_n_agent_tokens_exceeded(self, reaction):
        return (
            len(
                self.tokenizer.tokenize(
                    ".".join(reaction.get_agents_as_smiles())
                ).split(" ")
            )
            > self.max_number_of_agent_tokens
        )

    def max_n_product_tokens_exceeded(self, reaction):
        return (
            len(
                self.tokenizer.tokenize(
                    ".".join(reaction.get_products_as_smiles())
                ).split(" ")
            )
            > self.max_number_of_product_tokens
        )

    def formal_charge_exceeded(self, reaction: Reaction) -> bool:
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
        reactants = rdk.MolFromSmiles(".".join(reaction.get_reactants_as_smiles()))
        agents = rdk.MolFromSmiles(".".join(reaction.get_agents_as_smiles()))
        products = rdk.MolFromSmiles(".".join(reaction.get_products_as_smiles()))

        products_atoms = set([atom.GetSymbol() for atom in reactants.GetAtoms()])
        agents_atoms = set([atom.GetSymbol() for atom in agents.GetAtoms()])
        reactants_atoms = set([atom.GetSymbol() for atom in products.GetAtoms()])

        return len(products_atoms - (reactants_atoms | agents_atoms)) != 0

# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
""" Contains the class Reaction representing unidirectional reactions. """
from enum import auto
from enum import Enum
from typing import Any
from typing import Iterable
from typing import List
from typing import Set
from typing import Tuple

from rdkit.Chem import AllChem as rdk
from rdkit.Chem.rdchem import Mol


class ReactionPart(Enum):
    reactants = auto()
    agents = auto()
    products = auto()


class Reaction:

    def __init__(
        self,
        reaction_smarts: str,
        remove_duplicates: bool = False,
        fragment_bond='.',
        **kwargs: Any,
    ):
        """Creates a new instance of type Reaction based on a reaction SMARTs.

        Args:
            reaction_smarts: A reaction smarts
            remove_duplicates: Whether to remove duplicates from within reactants, agents and products. Defaults to False.
            kwargs: Keyword arguments supplied to rdkit MolToSmiles.
            fragment_bond: Token for the modeling of fragment bonds. Defaults to None.
        """
        kwargs.setdefault('canonical', True)

        self.__reaction_smarts = reaction_smarts
        self.__remove_duplicates = remove_duplicates
        self.__smiles_to_mol_kwargs = kwargs
        self.__fragment_bond = fragment_bond
        self.reactants, self.agents, self.products = self.__reaction_to_mols(
            self.__reaction_smarts
        )

    #
    # Overwrites / Virtuals
    #

    def __len__(self):
        """Returns the number of molecules participating in this reaction (reactants, agents, and products).

        Returns:
            The number of molecules participating in this reaction.
        """
        return len(self.reactants) + len(self.agents) + len(self.products)

    def __str__(self):
        """Returns the reaction SMARTS of this instance (reactants>agents>products).

        Returns:
            The reaction SMARTS representing this instance.
        """
        return (
            '.'.join(
                [
                    rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs
                                    ).replace('.', self.__fragment_bond)
                    for m in self.reactants if m
                ]
            ) + '>' + '.'.join(
                [
                    rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs
                                    ).replace('.', self.__fragment_bond) for m in self.agents if m
                ]
            ) + '>' + '.'.join(
                [
                    rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs
                                    ).replace('.', self.__fragment_bond)
                    for m in self.products if m
                ]
            )
        )

    def __eq__(self, other) -> bool:
        """Compares the count, order, and SMILES string of each molecule in this reaction.

        Args:
            other: A reaction to be compared with this instance.

        Returns:
            Whether this instance is equal to another.
        """
        if len(self) != len(other):
            return False

        if len(self.reactants) != len(other.reactants):
            return False

        if len(self.agents) != len(other.agents):
            return False

        if len(self.products) != len(other.products):
            return False

        # We care what the output of MolToSmiles is, so check for equality
        # on this. Inefficient, but allows us to mutate the molecule at
        # any time
        molecules_self = (self.reactants, self.agents, self.products)
        molecules_other = (other.reactants, other.agents, other.products)
        for i in range(len(molecules_self)):
            for a, b in zip(molecules_self[i], molecules_other[i]):
                if not self.__mols_equal(a, b):
                    return False

        return True

    #
    # Private Methods
    #

    def __reaction_to_mols(
        self,
        reaction_smarts: str,
    ) -> Tuple[List[Mol], List[Mol], List[Mol]]:
        """Creates a tuple of lists of reactants, agents, and products as rdkit Mol instances from a reaction SMARTS.

        Args:
            reaction_smarts: A reaction SMARTS.

        Raises:
            ValueError: This error is raised if a non-valid reaction SMARTS is provided.

        Returns:
            A tuple of lists of reactants, agents, and products representing the reaction.
        """
        if reaction_smarts.count('>') != 2:
            raise ValueError("A valid SMARTS reaction must contain two '>'.")

        raw_reactants, raw_agents, raw_products = tuple(reaction_smarts.split('>'))

        reactants = raw_reactants.split('.')
        agents = raw_agents.split('.')
        products = raw_products.split('.')

        if self.__remove_duplicates:
            reactants = list(dict.fromkeys(reactants))
            agents = list(dict.fromkeys(agents))
            products = list(dict.fromkeys(products))

        return (
            [
                rdk.MolFromSmiles(reactant.replace(self.__fragment_bond, '.'))
                for reactant in reactants if reactant != ''
            ],
            [
                rdk.MolFromSmiles(agent.replace(self.__fragment_bond, '.')) for agent in agents
                if agent != ''
            ],
            [
                rdk.MolFromSmiles(product.replace(self.__fragment_bond, '.'))
                for product in products if product != ''
            ],
        )

    def __mol_to_smiles(self, mol: Mol) -> str:
        """Encodes a molecule as a SMILES string by applying the rdkit MolToSmiles arguments supplied to this instantce.

        Args:
            mol: An rdkit Mol instance.

        Returns:
            The SMILES encoding of the input Mol.
        """
        return rdk.MolToSmiles(mol, **self.__smiles_to_mol_kwargs)

    def __mols_equal(self, first_mol: Mol, second_mol: Mol) -> bool:
        return self.__mol_to_smiles(first_mol) == self.__mol_to_smiles(second_mol)

    def __get_mols_as_smiles(self, mols: Iterable[Mol]) -> List[str]:
        return [rdk.MolToSmiles(mol, **self.__smiles_to_mol_kwargs) for mol in mols if mol]

    def __get_formal_charge_for_mols(self, mols: Iterable[Mol]) -> int:
        return sum(rdk.GetFormalCharge(mol) for mol in mols)

    def __get_atoms_for_mols(self, mols: Iterable[Mol]) -> Set[str]:
        """Returns the list of atoms, non repetitive, for a list of RDKit Mols.

        Returns:
            A list of tokens for the atoms.
        """
        return {atom.GetSymbol() for mol in mols for atom in mol.GetAtoms()}

    #
    # Pubilc Methods
    #

    def get_reactants_as_smiles(self) -> List[str]:
        """Returns the reactants of this reactions as a list of SMILES.

        Returns:
            A list of SMILES of the reactants.
        """
        return self.__get_mols_as_smiles(self.reactants)

    def get_agents_as_smiles(self) -> List[str]:
        """Returns the agents of this reactions as a list of SMILES.

        Returns:
            A list of SMILES of the agents.
        """
        return self.__get_mols_as_smiles(self.agents)

    def get_products_as_smiles(self) -> List[str]:
        """Returns the products of this reactions as a list of SMILES.

        Returns:
            A list of SMILES of the products.
        """
        return self.__get_mols_as_smiles(self.products)

    def get_reactants_formal_charge(self) -> int:
        """Returns the the formal charge for the reactants of this reaction.

        Returns:
            Formal charge of the reactants.
        """
        return self.__get_formal_charge_for_mols(self.reactants)

    def get_agents_formal_charge(self) -> int:
        """Returns the the formal charge for the agents of this reaction.

        Returns:
            Formal charge of the agents.
        """
        return self.__get_formal_charge_for_mols(self.agents)

    def get_products_formal_charge(self) -> int:
        """Returns the the formal charge for the products of this reaction.

        Returns:
            Formal charge of the products.
        """
        return self.__get_formal_charge_for_mols(self.products)

    def get_reactants_atoms(self) -> Set[str]:
        """Returns the list of atoms, non repetitive, for the reactants of this reaction as a Set of strings.

        Returns:
            A list of tokens for the atoms.
        """
        return self.__get_atoms_for_mols(self.reactants)

    def get_agents_atoms(self) -> Set[str]:
        """Returns the list of atoms, non repetitive, for the agents of this reaction as a Set of strings.

        Returns:
            A list of tokens for the atoms.
        """
        return self.__get_atoms_for_mols(self.agents)

    def get_products_atoms(self) -> Set[str]:
        """Returns the list of atoms, non repetitive, for the products of this reaction as a Set of strings.

        Returns:
            A list of tokens for the atoms.
        """
        return self.__get_atoms_for_mols(self.products)

    def find(self, pattern: str) -> Tuple[List[int], List[int], List[int]]:
        """Find the occurences of a SMARTS pattern within the reaction and returns a tuple of lists of indices in the reactants, agents, and products.

        Args:
            pattern: A SMARTS pattern.

        Returns:
            A tuple of lists of indices from the lists of reactants, agents, and products.
        """
        p = rdk.MolFromSmarts(pattern)

        # Avoid three method calls and do it directly
        return (
            [
                i for i, m in enumerate(self.reactants)
                if m and len(list(m.GetSubstructMatch(p))) > 0
            ],
            [i for i, m in enumerate(self.agents) if m and len(list(m.GetSubstructMatch(p))) > 0],
            [
                i for i, m in enumerate(self.products)
                if m and len(list(m.GetSubstructMatch(p))) > 0
            ],
        )

    def find_in(self, pattern: str, reaction_part: ReactionPart) -> List[int]:
        """Finds a SMARTS pattern in a part (reactants, agents, products) of the reaction.

        Args:
            pattern: A SMARTS pattern.
            reaction_part: The reaction part to search.

        Returns:
            A list of indices from the list of molecules representing the chosen reaction part.
        """
        p = rdk.MolFromSmarts(pattern)

        if reaction_part == ReactionPart.reactants:
            return [
                i for i, m in enumerate(self.reactants)
                if m and len(list(m.GetSubstructMatch(p))) > 0
            ]

        if reaction_part == ReactionPart.agents:
            return [
                i for i, m in enumerate(self.agents) if m and len(list(m.GetSubstructMatch(p))) > 0
            ]

        if reaction_part == ReactionPart.products:
            return [
                i for i, m in enumerate(self.products)
                if m and len(list(m.GetSubstructMatch(p))) > 0
            ]

        return []

    # Annotate the return type as Reaction... Python...
    def remove(self, indices: Tuple[List[int], List[int], List[int]]):
        """Remove reactants, agents and products based on their index within the respective lists.

        Args:
            indices: The indices of the molecules to be removed from the reaction.

        Returns:
            Itself with changes applied.
        """
        if len(indices) > 0:
            for idx in sorted(indices[0], reverse=True):
                del self.reactants[idx]

        if len(indices) > 1:
            for idx in sorted(indices[1], reverse=True):
                del self.agents[idx]

        if len(indices) > 2:
            for idx in sorted(indices[2], reverse=True):
                del self.products[idx]

        return self

    def filter(self, indices: Tuple[List[int], List[int], List[int]]):
        """Filter for reactants, agents and products based on their index within the respective lists. This is the complement to remove.

        Args:
            indices: The indices of the molecules to not be removed from the reaction.

        Returns:
            Itself with changes applied.
        """
        if len(indices) > 0 and len(indices[0]) > 0:
            for idx in range(len(self.reactants) - 1, -1, -1):
                if idx not in indices[0]:
                    del self.reactants[idx]

        if len(indices) > 1 and len(indices[1]) > 0:
            for idx in range(len(self.agents) - 1, -1, -1):
                if idx not in indices[1]:
                    del self.agents[idx]

        if len(indices) > 2 and len(indices[2]) > 0:
            for idx in range(len(self.products) - 1, -1, -1):
                if idx not in indices[2]:
                    del self.products[idx]

        return self

    def sort(self, sort_reactants=True, sort_agents=True, sort_products=True):
        """Order the molecules participating in this reaction based on their SMILES strings.
           The rdkit MolToSmiles argument supplied to this instance will be applied.

        Args:
            sort_reactants: Whether to sort the reactants. Defaults to True.
            sort_agents: Whether to sort the agents. Defaults to True.
            sort_products: Whether to sort the products. Defaults to True.

        Returns:
            Itself with changes applied.
        """
        # Mol to SMILES here is rather inefficient, but this allows for
        # changes to the Mol objects at any time in the lifecycle of the instance
        if sort_reactants:
            self.reactants = sorted(
                self.reactants,
                key=lambda m: rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs),
            )

        if sort_agents:
            self.agents = sorted(
                self.agents,
                key=lambda m: rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs),
            )

        if sort_products:
            self.products = sorted(
                self.products,
                key=lambda m: rdk.MolToSmiles(m, **self.__smiles_to_mol_kwargs),
            )

        return self

    def remove_precursors_from_products(self):
        """Removes prodcuts that are also found in reactants or agents.

        Returns:
            Itself with prodcuts occuring as reactants or agents removed.
        """

        reactants_smiles = self.get_reactants_as_smiles()
        agents_smiles = self.get_agents_as_smiles()
        products_smiles = self.get_products_as_smiles()

        for i, product in reversed(list(enumerate(products_smiles))):
            if product in reactants_smiles or product in agents_smiles:
                del self.products[i]

        return self

    def has_none(self) -> bool:
        """Checks whether the reactants, agents, or products contain None (usually due to failed rdkit MolFromSmiles).

        Returns:
            Whether the reactants, agents, or products contain None.
        """

        return None in self.reactants or None in self.agents or None in self.products

    def remove_none(self):
        """Removes all None values from the reactants, agents, and products.

        Returns:
            Itself with None values removed.
        """
        self.reactants = [m for m in self.reactants if m is not None]
        self.agents = [m for m in self.agents if m is not None]
        self.products = [m for m in self.products if m is not None]

        return self

    #
    # Static Methods
    #

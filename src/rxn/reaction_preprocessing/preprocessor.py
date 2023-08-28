# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
""" The preprocessor class abstracts the workflow for preprocessing reaction data sets. """
import collections
import logging
from pathlib import Path
from typing import Counter, Iterator, List, Tuple

from rxn.chemutils.reaction_equation import ReactionEquation
from rxn.chemutils.reaction_smiles import parse_any_reaction_smiles
from rxn.utilities.containers import iterate_unique_values
from rxn.utilities.csv import CsvIterator, StreamingCsvEditor
from rxn.utilities.files import PathLike
from tabulate import tabulate

from .config import PreprocessConfig
from .mixed_reaction_filter import MixedReactionFilter
from .reaction_standardizer import ReactionStandardizer

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class Preprocessor:
    def __init__(
        self,
        mixed_reaction_filter: MixedReactionFilter,
        reaction_column_name: str,
        fragment_bond: str = ".",
    ):
        """
        Args:
            mixed_reaction_filter: mixed reaction filter.
            reaction_column_name: The name of the DataFrame column containing the reaction SMARTS.
            fragment_bond: The token that represents fragment bonds in the reaction SMILES.
        """
        self.reaction_standardizer = ReactionStandardizer()
        self.mixed_reaction_filter = mixed_reaction_filter
        self.rxn_column = reaction_column_name
        self.fragment_bond = fragment_bond

    def process(self, input_file_path: PathLike, output_file_path: PathLike) -> None:
        with open(input_file_path, "rt") as f_in, open(output_file_path, "wt") as f_out:
            csv_iterator = CsvIterator.from_stream(f_in)

            # first deduplication
            csv_iterator = self.remove_duplicate_reactions(csv_iterator)

            # reaction standardization
            rxn_standardization_editor = StreamingCsvEditor(
                [self.rxn_column], [self.rxn_column], self.standardize_rxn_smiles
            )
            csv_iterator = rxn_standardization_editor.process(csv_iterator)

            # second deduplication
            csv_iterator = self.remove_duplicate_reactions(csv_iterator)

            # filtering
            csv_iterator, error_counter = self.validate(
                csv_iterator=csv_iterator,
            )

            csv_iterator.to_stream(f_out)

            self.print_stats(
                error_counter=error_counter,
                original_count=-1,
                valid_count=-1,
                invalid_count=-1,
            )

    def remove_duplicate_reactions(self, csv_iterator: CsvIterator) -> CsvIterator:
        rxn_idx = csv_iterator.column_index(self.rxn_column)

        # The key for determining what is a duplicate is the value from the rxn column
        def key(row: List[str]) -> str:
            return row[rxn_idx]

        return CsvIterator(
            csv_iterator.columns, iterate_unique_values(csv_iterator.rows, key=key)
        )

    def standardize_rxn_smiles(self, rxn_smiles: str) -> str:
        """Function standardizing the reaction SMILES directly,
        to pass to pandas.apply()."""
        reaction = parse_any_reaction_smiles(rxn_smiles)
        reaction = self.reaction_standardizer(reaction)
        return reaction.to_string(self.fragment_bond)

    def validate(
        self,
        csv_iterator: CsvIterator,
    ) -> Tuple[CsvIterator, Counter[str]]:
        rxn_idx = csv_iterator.column_index(self.rxn_column)
        error_counter: Counter[str] = collections.Counter()

        def internal() -> Iterator[List[str]]:
            for row in csv_iterator.rows:
                reaction = ReactionEquation.from_string(
                    row[rxn_idx], fragment_bond=self.fragment_bond
                )
                valid, reasons = self.mixed_reaction_filter.validate_reasons(reaction)
                if valid:
                    yield row
                else:
                    for reason in reasons:
                        error_counter[reason] += 1

        return CsvIterator(columns=csv_iterator.columns, rows=internal()), error_counter

    def print_stats(
        self,
        error_counter: Counter[str],
        original_count: int,
        valid_count: int,
        invalid_count: int,
    ) -> None:
        """Prints statistics of the filtration to the logger."""
        logger.info(f"- {original_count} total reactions.")
        logger.info(f"- {valid_count} valid reactions.")
        logger.info(f"- {invalid_count} invalid reactions removed.")

        if sum(error_counter.values()) == 0:
            return

        headers = ["Reason", "Number of Reactions"]
        logger.info(
            f"- The {invalid_count} reactions were removed for the following reasons:\n"
            f'{tabulate(list(error_counter.items()), headers, tablefmt="fancy_grid")}'
        )


def preprocess(cfg: PreprocessConfig) -> None:
    output_file_path = Path(cfg.output_file_path)
    input_file_path = Path(cfg.input_file_path)
    if not input_file_path.exists():
        raise ValueError(
            f"Input file for preprocessing does not exist: {input_file_path}"
        )

    # Create a instance of the mixed reaction filter with default values.
    # Make arguments for all properties in script
    mrf = MixedReactionFilter(
        max_reactants=cfg.max_reactants,
        max_agents=cfg.max_agents,
        max_products=cfg.max_products,
        min_reactants=cfg.min_reactants,
        min_agents=cfg.min_agents,
        min_products=cfg.min_products,
        max_reactants_tokens=cfg.max_reactants_tokens,
        max_agents_tokens=cfg.max_agents_tokens,
        max_products_tokens=cfg.max_products_tokens,
        max_absolute_formal_charge=cfg.max_absolute_formal_charge,
    )

    rxn_column = cfg.reaction_column_name
    preprocessor = Preprocessor(
        mixed_reaction_filter=mrf,
        reaction_column_name=rxn_column,
        fragment_bond=cfg.fragment_bond.value,
    )

    preprocessor.process(input_file_path, output_file_path)

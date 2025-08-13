# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2022
# ALL RIGHTS RESERVED
"""The preprocessor class abstracts the workflow for preprocessing reaction data sets."""
import collections
import logging
from pathlib import Path
from typing import Callable, Counter, Iterable, Iterator, List

import attr
from attr import define
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


@define
class _Stats:
    initial_count: int = 0
    first_dedup_count: int = 0
    second_dedup_count: int = 0
    final_count: int = 0
    error_counter: Counter[str] = attr.Factory(collections.Counter)


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
        self.stats = _Stats()

    def process_file(
        self, input_file_path: PathLike, output_file_path: PathLike
    ) -> None:
        """Process the reactions in a CSV file.

        Args:
            input_file_path: CSV with the reactions to standardize.
            output_file_path: CSV where to save the standardized reactions.
        """

        with open(input_file_path, "rt") as f_in, open(output_file_path, "wt") as f_out:
            csv_iterator = CsvIterator.from_stream(f_in)
            csv_iterator = self.process_iterator(csv_iterator)
            csv_iterator.to_stream(f_out)

    def process_iterator(self, csv_iterator: CsvIterator) -> CsvIterator:
        """Process the reactions in a CSV iterator.

        Same as ``process_file``, except that it acts directly on the iterator.

        Args:
            csv_iterator: input CSV iterator for the reactions to process.

        Returns:
            CsvIterator with reactions after the processor step.
        """
        # Reset the stats
        self.stats = _Stats()

        csv_iterator = self._increment_initial_count(csv_iterator)

        # first deduplication
        csv_iterator = self._remove_duplicate_reactions(csv_iterator)
        csv_iterator = self._increment_first_dedup_count(csv_iterator)

        # reaction standardization
        rxn_standardization_editor = StreamingCsvEditor(
            [self.rxn_column], [self.rxn_column], self._standardize_rxn_smiles
        )
        csv_iterator = rxn_standardization_editor.process(csv_iterator)

        # second deduplication
        csv_iterator = self._remove_duplicate_reactions(csv_iterator)
        csv_iterator = self._increment_second_dedup_count(csv_iterator)

        # filtering
        csv_iterator = self._validate(csv_iterator=csv_iterator)

        csv_iterator = self._increment_final_count(csv_iterator)

        return csv_iterator

    def _increment_initial_count(self, csv_iterator: CsvIterator) -> CsvIterator:
        def fn() -> None:
            self.stats.initial_count += 1

        return self._call_for_each(csv_iterator, fn)

    def _increment_first_dedup_count(self, csv_iterator: CsvIterator) -> CsvIterator:
        def fn() -> None:
            self.stats.first_dedup_count += 1

        return self._call_for_each(csv_iterator, fn)

    def _increment_second_dedup_count(self, csv_iterator: CsvIterator) -> CsvIterator:
        def fn() -> None:
            self.stats.second_dedup_count += 1

        return self._call_for_each(csv_iterator, fn)

    def _increment_final_count(self, csv_iterator: CsvIterator) -> CsvIterator:
        def fn() -> None:
            self.stats.final_count += 1

        return self._call_for_each(csv_iterator, fn)

    def _call_for_each(
        self, csv_iterator: CsvIterator, fn: Callable[[], None]
    ) -> CsvIterator:
        """Allows to call a callback before accessing the CSV row.

        Useful for counting, for instance.
        """

        def iterate(values: Iterable[List[str]]) -> Iterator[List[str]]:
            for value in values:
                fn()
                yield value

        return CsvIterator(
            columns=csv_iterator.columns,
            rows=iterate(csv_iterator.rows),
        )

    def _remove_duplicate_reactions(self, csv_iterator: CsvIterator) -> CsvIterator:
        rxn_idx = csv_iterator.column_index(self.rxn_column)

        # The key for determining what is a duplicate is the value from the rxn column
        def key(row: List[str]) -> str:
            return row[rxn_idx]

        return CsvIterator(
            csv_iterator.columns, iterate_unique_values(csv_iterator.rows, key=key)
        )

    def _standardize_rxn_smiles(self, rxn_smiles: str) -> str:
        """Standardizing the reaction SMILES.

        If there is an error, returns ">>" instead.
        """
        try:
            reaction = parse_any_reaction_smiles(rxn_smiles)
            reaction = self.reaction_standardizer(reaction)
            return reaction.to_string(self.fragment_bond)
        except Exception as e:
            logger.error(f'Cannot standardize reaction SMILES "{rxn_smiles}": {e}')
            return ">>"

    def _validate(self, csv_iterator: CsvIterator) -> CsvIterator:
        rxn_idx = csv_iterator.column_index(self.rxn_column)

        def filter_invalid(rows: Iterable[List[str]]) -> Iterator[List[str]]:
            for row in rows:
                reaction = ReactionEquation.from_string(
                    row[rxn_idx], fragment_bond=self.fragment_bond
                )
                valid, reasons = self.mixed_reaction_filter.validate_reasons(reaction)
                if valid:
                    yield row
                else:
                    for reason in reasons:
                        self.stats.error_counter[reason] += 1

        return CsvIterator(
            columns=csv_iterator.columns, rows=filter_invalid(csv_iterator.rows)
        )

    def print_stats(self) -> None:
        """Prints statistics of the filtration to the logger."""
        # define "s" to make expressions shorter
        s = self.stats

        logger.info(f"- {s.initial_count} total reactions.")
        logger.info(f"- {s.first_dedup_count} reactions after first deduplication.")
        logger.info(f"- {s.second_dedup_count} reactions after second deduplication.")
        logger.info(f"- {s.final_count} valid reactions.")
        invalid_count = s.second_dedup_count - s.final_count
        logger.info(f"- {invalid_count} invalid reactions removed.")

        if sum(s.error_counter.values()) == 0:
            return

        headers = ["Reason", "Number of Reactions"]
        logger.info(
            f"- The {invalid_count} reactions were removed for the following reasons:\n"
            f'{tabulate(s.error_counter.most_common(), headers, tablefmt="fancy_grid")}'
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

    preprocessor.process_file(input_file_path, output_file_path)
    preprocessor.print_stats()

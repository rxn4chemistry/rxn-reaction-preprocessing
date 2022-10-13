import click
from rxn.utilities.files import count_lines, iterate_lines_from_file


@click.command()
@click.option(
    "--tokenized_products",
    "-p",
    required=True,
    help="File containing the tokenized products",
)
@click.option("--class_ids", "-c", required=True, help="File containing the class ids.")
@click.option("--output_file", "-o", required=True, help="Resulting file name.")
@click.option(
    "--no_brackets",
    is_flag=True,
    help='If given, "[" and "]" will not be added around the class id.',
)
def main(
    tokenized_products: str, class_ids: str, output_file: str, no_brackets: bool
) -> None:
    """Combine the class tokens with the file with the tokenized products."""

    n_products = count_lines(tokenized_products)
    n_class_ids = count_lines(class_ids)

    if n_products != n_class_ids:
        raise ValueError(
            f"The number of products ({n_products}) and the number of class "
            f"ids ({n_class_ids}) must be identical."
        )

    print(f"Adding class id prefixes for {n_products} products.")

    product_iterator = iterate_lines_from_file(tokenized_products)
    class_id_iterator = iterate_lines_from_file(class_ids)

    with open(output_file, "wt") as f:
        for product, class_id in zip(product_iterator, class_id_iterator):
            if no_brackets:
                prefix = f"{class_id} "
            else:
                prefix = f"[{class_id}] "
            f.write(prefix + product + "\n")


if __name__ == "__main__":
    main()

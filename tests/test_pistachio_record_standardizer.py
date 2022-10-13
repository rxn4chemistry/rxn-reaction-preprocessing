from typing import Any, Dict

import pytest
from rxn.chemutils.exceptions import InvalidSmiles

from rxn.reaction_preprocessing.annotations.molecule_annotation import (
    MoleculeAnnotation,
)
from rxn.reaction_preprocessing.config import (
    FragmentBond,
    PreprocessConfig,
    StandardizeConfig,
)
from rxn.reaction_preprocessing.mixed_reaction_filter import ReactionFilterError
from rxn.reaction_preprocessing.molecule_standardizer import (
    MoleculeStandardizer,
    RejectedMolecule,
)
from rxn.reaction_preprocessing.pistachio_record_standardizer import (
    PistachioRecordStandardizer,
)


def create_dummy_record(
    reaction_smiles: str, component_0: str, component_2: str, action_component: str
) -> Dict[str, Any]:
    """
    Create a fake reaction record with the given reaction SMILES, two component
    SMILES strings, as well as an action component SMILES string.

    Args:
        reaction_smiles: reaction SMILES.
        component_0: component SMILES that will end up at index 0.
        component_2: component SMILES that will end up at index 2.
        action_component: SMILES for the component of an action.
    """
    # Reaction record, slightly stripped down - still includes the reaction SMILES,
    # some components, as well as some components linked to actions, and other unnecessary data.
    return {
        "_id": "5dd43d06a1af77d2f23351ad",
        "title": "US08865712B2_0629",
        "data": {
            "headingText": "b) 1-(6-Bromo-3-fluoro-4-triethylsilanyl-pyridin-2-yl)-ethanone",
            "smiles": reaction_smiles,
            "ipcCodes": ["A61K 31/5377", "C07D 413/14", "C07D 471/04", "C07D 487/04"],
        },
        "components": [
            {
                "role": "Product",
                "name": "1-(6-Bromo-3-fluoro-4-triethylsilanyl-pyridin-2-yl)-ethanone",
                "smiles": component_0,
                "quantities": [{"type": "Mass", "value": 58.5, "text": "58.5 g"}],
            },
            {
                "role": "Reactant",
                "name": "petroleum ether",  # NB: no SMILES
            },
            {
                "role": "Solvent",
                "name": "water",
                "smiles": component_2,
                "inchi": "InChI=1S/H2O/h1H2",
                "quantities": [{"type": "Volume", "value": 0.25, "text": "250 ml"}],
            },
        ],
        "actions": [
            {
                "type": "Add",
                "text": "BuLi (100 ml, 250 mmol, 2.5 M in hexanes) was added dropwise",
                "components": [
                    {
                        "role": "Agent",
                        "name": "BuLi",
                        "smiles": action_component,
                        "inchi": "InChI=1S/C4H9.Li/c1-3-4-2;/h1,3-4H2,2H3;",
                        "quantities": [
                            {"type": "Amount", "value": 0.25, "text": "250 mmol"}
                        ],
                    }
                ],
            }
        ],
        "source": "grants/2014/I20141021_txt_reactions_wbib.json",
    }


def test_pistachio_record_standardizer() -> None:
    # Set up the standardizer, and overwrite its annotations manually
    config_standardize = StandardizeConfig(
        fragment_bond=FragmentBond.TILDE, annotation_file_paths=[]
    )
    config_preprocess = PreprocessConfig(fragment_bond=FragmentBond.TILDE)
    standardizer = PistachioRecordStandardizer(config_standardize, config_preprocess)
    # overwrite the molecule standardizer to allow for custom annotations
    standardizer.molecule_standardizer = MoleculeStandardizer(
        annotations=[
            MoleculeAnnotation("[Cl-].[Na+]", "[Na]Cl", "accept", []),
            MoleculeAnnotation("O.[BH4-].[Na+]", "O.[BH4-]~[Na+]", "accept", []),
            MoleculeAnnotation("[Pd]", None, "reject", []),
        ],
        discard_missing_annotations=config_standardize.discard_unannotated_metals,
    )

    # Test 1:
    # Nothing to update (already in standard form) -> the record is identical before
    # and after standardization, while being a different Python object.
    record = create_dummy_record("CC.O>>CCO", "CC", "O", "CC")
    updated_record = standardizer.standardize(record)
    assert record is not updated_record
    assert record == updated_record

    # Test 2:
    # Standardizes the salt, and also updates it in the components
    record = create_dummy_record("[Na+]~[Cl-].CC.O>>CCO", "[Na+].[Cl-]", "O", "CC")
    updated_record = standardizer.standardize(record)
    assert updated_record == create_dummy_record(
        "CC.O.[Na]Cl>>CCO", "[Na]Cl", "O", "CC"
    )

    # Test 3:
    # same as test 2, with the input in extended SMILES format, and atom mapping
    record = create_dummy_record(
        "[Na+].[Cl-].[CH3:1][CH3:2].[OH2:3]>>[CH3:1][CH2:2][OH:3] |f:0.1|",
        "[Na+].[Cl-]",
        "O",
        "CC",
    )
    updated_record = standardizer.standardize(record)
    assert updated_record == create_dummy_record(
        "CC.O.[Na]Cl>>CCO", "[Na]Cl", "O", "CC"
    )

    # Test 4:
    # Similar to test 2, except that here one component is split into two.
    # Note how in the components it looks like they still belong together (not ideal).
    record = create_dummy_record("[Na+]~O~[BH4-].CC>>CCO", "[Na+].O.[BH4-]", "O", "OCC")
    updated_record = standardizer.standardize(record)
    assert updated_record == create_dummy_record(
        "CC.O.[BH4-]~[Na+]>>CCO", "O.[BH4-].[Na+]", "O", "CCO"
    )

    # Test 5:
    # One of the compounds in the reaction SMILES is rejected
    with pytest.raises(RejectedMolecule):
        _ = standardizer.standardize(create_dummy_record("[Pd].O.C>>CO", "C", "O", "C"))

    # Test 6:
    # One molecule in the reaction SMILES cannot be canonicalized -> fail
    with pytest.raises(InvalidSmiles):
        _ = standardizer.standardize(create_dummy_record("Q.O.C>>CO", "C", "O", "C"))

    # Test 7:
    # One molecule in the components cannot be canonicalized -> the molecule is not updated
    record = standardizer.standardize(create_dummy_record("O.C>>CO", "Q", "O", "C"))
    assert record["components"][0]["smiles"] == "Q"

    # Test 8:
    # No product -> fail
    with pytest.raises(ReactionFilterError):
        _ = standardizer.standardize(create_dummy_record("O.C>>", "C", "O", "C"))

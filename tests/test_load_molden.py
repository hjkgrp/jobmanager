import pytest
import json
import numpy as np
from jobmanager.io import load_molden


@pytest.mark.parametrize(
    "name",
    [
        "terachem_fe2_co_6_restricted",
        "terachem_fe2_co_6_unrestricted"
    ]
)
def test_load_molden(resource_path_root, name, atol=1e-5):
    molden = load_molden(
        resource_path_root / "io" / "molden" / f"{name}.molden"
        )

    with open(resource_path_root / "io" / "molden" / f"{name}.json",
              "r") as fin:
        ref = json.load(fin)

    assert molden.keys() == ref.keys()

    # Compare title
    assert molden["title"] == ref["title"]
    # Compare the numpy arrays
    for key, values in molden.items():
        if key == "obasis":
            assert values.keys() == ref[key].keys()
            for subkey, subvalues in values.items():
                np.testing.assert_allclose(
                    subvalues, ref[key][subkey], atol=atol
                    )
        else:
            if key not in ["title", "permutation"]:
                np.testing.assert_allclose(values, ref[key], atol=atol)

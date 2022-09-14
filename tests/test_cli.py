import subprocess
import pytest


def test_cli_call():
    # Expecting an "unknown machine" value error if the command line
    # script is working correctly
    with pytest.raises(ValueError):
        subprocess.run(['jobmanager'])
import subprocess
import pytest


def test_cli_call():
    output = subprocess.run(['jobmanager'])
    print(output)
    assert False

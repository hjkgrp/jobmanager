import subprocess


def test_cli_call():
    output = subprocess.run(['jobmanager'])
    assert output.returncode == 0
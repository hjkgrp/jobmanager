import subprocess


def test_cli_call():
    # TODO: there might be a more elegant way of checking that the entry
    # point is registered correctly.
    output = subprocess.run(["jobmanager"], capture_output=True)
    # This is expected to fail because the jobmanager does not
    # recognize the machine it is running on:
    assert b"ValueError: Machine Unknown to Job Manager" in output.stderr

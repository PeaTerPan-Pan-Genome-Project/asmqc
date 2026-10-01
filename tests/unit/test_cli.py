import pytest

from asmqc import cli


def test_version_exits_zero(capsys):
    assert cli.main(["version"]) == 0
    assert capsys.readouterr().out.startswith("asmqc ")


def test_command_required():
    with pytest.raises(SystemExit):
        cli.main([])

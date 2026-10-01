import json

import pytest

from asmqc import cli


def test_version_exits_zero(capsys):
    assert cli.main(["version"]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["tools"]["samtools"] == "1.24"


def test_command_required():
    with pytest.raises(SystemExit):
        cli.main([])

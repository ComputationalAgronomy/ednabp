import argparse
import os
from unittest.mock import patch

import pytest

from ednabp.cli import pipeline as cli
from ednabp.common.default_settings import SETTINGS


def setting_names():
    names = {"input_path", "output_path", "enabled_stages", "raw_suffix"}
    for group, values in SETTINGS.items():
        if group in ("enabled_stages", "dir_name", "suffix") or not isinstance(
            values, dict
        ):
            continue
        if group == "prog":
            names |= {f"{k}_prog" for k in values}
        else:
            names |= set(values)
    return names


def parse(argv):
    with patch("ednabp.bp.BioPipeline") as pipeline:
        cli.main(argv)
    return pipeline.call_args.kwargs


def defaults():
    out = {}
    for group, values in SETTINGS.items():
        if group in ("enabled_stages", "dir_name", "suffix") or not isinstance(
            values, dict
        ):
            continue
        for k, v in values.items():
            out[f"{k}_prog" if group == "prog" else k] = v
    return out


def test_options_reach_the_pipeline(tmp_dirs):
    in_dir, out_dir = tmp_dirs
    options = parse(
        ["-i", in_dir, "-o", out_dir, "-e", "qfilter", "--min-qscore", "25"]
    )
    assert options["input_path"] == in_dir
    assert options["output_path"] == out_dir
    assert options["enabled_stages"] == ["qfilter"]
    assert options["min_qscore"] == 25.0
    expected = defaults()
    for name, value in options.items():
        if name in expected and name != "min_qscore":
            assert value == expected[name], name


def test_every_option_is_a_setting(tmp_dirs):
    """Each CLI option must name a BioPipeline setting, so a renamed or
    removed setting can't leave a dead option behind."""
    in_dir, out_dir = tmp_dirs
    dests = set()
    real_add = argparse._ActionsContainer.add_argument

    def record(self, *args, **kwargs):
        action = real_add(self, *args, **kwargs)
        dests.add(action.dest)
        return action

    with (
        patch("argparse._ActionsContainer.add_argument", record),
        patch("ednabp.bp.BioPipeline"),
    ):
        cli.main(["-i", in_dir, "-o", out_dir])

    assert dests - {"help"} - setting_names() == set()


def test_help_exits_cleanly(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["--help"])
    assert exit_info.value.code == 0
    assert "Quality Filter Settings" in capsys.readouterr().out


def test_missing_input_is_a_usage_error(tmp_dirs, capsys):
    _, out_dir = tmp_dirs
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["-i", os.path.join(out_dir, "nope"), "-o", out_dir])
    assert exit_info.value.code == 2
    assert "input path does not exist" in capsys.readouterr().err


def test_taxids_and_dry_options(tmp_dirs):
    in_dir, out_dir = tmp_dirs
    options = parse(
        ["-i", in_dir, "-o", out_dir, "--taxids", "2759,10239", "--dry"]
    )
    assert options["taxids"] == "2759,10239"
    assert options["dry"] is True
    defaults = parse(["-i", in_dir, "-o", out_dir])
    assert defaults["taxids"] is None and defaults["dry"] is False

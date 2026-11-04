"""Tests for the standalone tools, their dispatch and the shell completion."""

import shutil
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits

from specproc.main import (_profile_from_inputs, _profile_from_metadata,
                           copy_default_config_to_cwd, main)
from specproc.utils import completion
from specproc.utils.tool_registry import TOOLS, format_tool_table, get_tool, tool_names


def test_registry_contents():
    assert list(tool_names()) == ['air2vac', 'vac2air', 'identify']
    assert get_tool('match-anchors') is None
    assert get_tool('identify').module == 'specproc.utils.identify'
    assert get_tool('nope') is None


def test_tool_table():
    table = format_tool_table()
    for tool in TOOLS:
        assert tool.name in table
        assert tool.summary in table
        assert tool.usage in table
    assert 'specproc.air2vac INPUT OUTPUT' in table
    assert 'tools' in table


def test_registry_import_is_lazy():
    """Listing tools must not drag in PyQt/matplotlib."""
    code = (
        "import sys, specproc.utils.tool_registry as r; "
        "assert 'PyQt5' not in sys.modules, 'PyQt5 imported eagerly'; "
        "print(','.join(r.tool_names()))"
    )
    proc = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == 'air2vac,vac2air,identify'


def test_help_lists_tools(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(argv=['--help'])
    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    assert 'Command-line tools' in out
    assert 'specproc --completion bash' in out
    for name in tool_names():
        assert name in out


def test_tools_subcommand(capsys):
    assert main(argv=['tools']) == 0
    out = capsys.readouterr().out
    assert 'identify' in out
    assert 'specproc.<tool> --help' in out


def test_dispatch_air2vac(tmp_path):
    source = tmp_path / 'lines.dat'
    source.write_text('5000.0 Fe I\ncomment line\n6000.0 Fe I 2\n')
    target = tmp_path / 'converted.dat'

    assert main(argv=['air2vac', str(source), str(target)]) == 0

    lines = [ln for ln in target.read_text().splitlines() if not ln.startswith('#')]
    assert lines[0].endswith('Fe I')
    assert float(lines[0].split()[0]) > 5000.0  # vacuum > air
    assert lines[1] == 'comment line'
    assert float(lines[2].split()[0]) > 6000.0


def test_complete_top_level():
    assert completion.complete(1, ['specproc', '']) == [
        'tools', 'air2vac', 'vac2air', 'identify']
    assert '--config' in completion.complete(1, ['specproc', '--'])


def test_complete_option_values():
    assert completion.complete(2, ['specproc', '--mode', '']) == ['gui', 'cli']
    assert completion.complete(2, ['specproc', '--completion', '']) == ['bash', 'zsh']
    profiles = completion.complete(2, ['specproc', '--config', ''])
    assert 'MSTG193_slit' in profiles and 'XL216_hrs' in profiles
    # a path-like word is left to filename completion
    assert completion.complete(2, ['specproc', '--config', './raw/']) == []


def test_complete_tool_arguments():
    options = completion.complete(2, ['specproc', 'identify', '--'])
    assert '--lamp' in options and '--ref' in options
    choices = completion.complete(3, ['specproc', 'identify', '--poly-order', ''])
    assert choices == [str(k) for k in range(1, 10)]
    types = completion.complete(3, ['specproc', 'identify', '--poly-type', ''])
    assert types == ['legendre', 'chebyshev', 'polynomial']


def test_complete_after_tools_namespace():
    assert 'air2vac' in completion.complete(2, ['specproc', 'tools', ''])
    assert '--lamp' in completion.complete(3, ['specproc', 'tools', 'identify', '--'])


def test_complete_does_not_shadow_data_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'air2vac').mkdir()
    # an existing path is pipeline input, not the tool: only the top-level
    # options are offered, never the tool's own arguments.
    candidates = completion.complete(2, ['specproc', 'air2vac', '--'])
    assert '--config' in candidates
    assert '--lamp' not in candidates


def test_complete_protocol(capsys):
    assert main(argv=['--_complete', '1', 'specproc', '']) == 0
    out = capsys.readouterr().out.split()
    assert 'identify' in out and 'tools' in out


def test_completion_scripts_are_valid(tmp_path):
    cases = [('bash', completion.bash_script())]
    if shutil.which('zsh'):
        cases.append(('zsh', completion.zsh_script()))
    for shell, script in cases:
        path = tmp_path / f'completion.{shell}'
        path.write_text(script)
        proc = subprocess.run([shell, '-n', str(path)], capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr

def test_profile_from_fixed_telescope_metadata(tmp_path):
    assert _profile_from_metadata('MOST    ', 'Camera Spectrum-R') == 'mstg193slit'
    assert _profile_from_metadata('MOST    ', 'Camera Spectrum-B') == 'mstg193slit'
    assert _profile_from_metadata('XINGLONG 2.16m', 'HRS') == 'xinglong216hrs'
    assert _profile_from_metadata('unknown', 'camera') is None

    path = tmp_path / 'lamp.fits'
    hdu = fits.PrimaryHDU(np.ones(8))
    hdu.header['TELESCOP'] = 'MOST    '
    hdu.header['INSTRUME'] = 'Camera Spectrum-R'
    hdu.writeto(path)
    assert _profile_from_inputs([str(tmp_path)]) == 'mstg193slit'


def test_nanshan12_profile_from_swowner(tmp_path):
    from specproc.config.config_manager import resolve_telescope_profile

    path = tmp_path / 'obj.fit'
    hdu = fits.PrimaryHDU(np.ones(8))
    hdu.header['SWOWNER'] = 'nanshan '
    hdu.writeto(path)
    assert _profile_from_inputs([str(tmp_path)]) == 'nanshan12'
    assert resolve_telescope_profile('nanshan12').name == 'default_config_NS12_hrs.cfg'


def test_mstg_anchor_file_loads(tmp_path):
    from specproc.config.config_manager import ConfigManager
    from specproc.core.slit_pipeline import SlitPipeline

    config = ConfigManager('MSTG193_slit')
    pipeline = SlitPipeline(config)
    pairs = pipeline._load_id_file()
    assert pairs is not None and len(pairs) >= 10
    assert np.all(np.isfinite(pairs[:, :2]))


def test_copy_mstg_profile_to_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    path = copy_default_config_to_cwd(profile='MSTG193_slit')
    assert path == tmp_path / 'default_config_MSTG193_slit.cfg'
    text = path.read_text()
    assert 'name = mstg193slit' in text
    assert 'mstg193slit_anchors_R.csv' in text
    # Never overwrite an existing working config.
    assert copy_default_config_to_cwd(profile='MSTG193_slit') == path


def test_copy_generic_config_keeps_working_name(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert copy_default_config_to_cwd() == tmp_path / 'specproc.cfg'


def test_complete_dotted_command():
    assert '--auto' in completion.complete(1, ['specproc.identify', '--'])
    assert completion.complete(1, ['specproc.nope', '--']) == []
    assert 'specproc.identify' in completion.bash_script()

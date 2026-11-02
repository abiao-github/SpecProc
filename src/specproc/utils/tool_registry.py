"""Registry of the standalone SpecProc command-line tools.

The data-reduction pipeline lives in :mod:`specproc.core` and is driven by the
``specproc`` GUI or ``specproc --mode cli``.  The tools listed here are *not*
part of that pipeline: they are small, self-contained commands that prepare
or inspect the data (wavelength conversion, lamp-line identification, anchor
matching).

Each tool can be run as ``specproc <name> ...``; for backward compatibility it
is also installed as its own ``specproc.<name>`` console script.
"""

from dataclasses import dataclass, field
from typing import Dict, Optional, Sequence, Tuple


@dataclass(frozen=True)
class Tool:
    """Description of one command-line tool.

    Attributes:
        name: Subcommand used as ``specproc <name>``.
        module: Import path of the module implementing the tool.
        entry: Callable in ``module`` with signature ``entry(argv=None, prog=None)``.
        summary: One-line description shown in ``specproc --help``.
        parser: Name of the factory that builds the tool's
            :class:`argparse.ArgumentParser` (used by shell completion).
        parser_kwargs: Extra keyword arguments for that factory.
    """

    name: str
    module: str
    entry: str
    summary: str
    usage: str = ''
    parser: str = 'build_parser'
    parser_kwargs: Dict[str, object] = field(default_factory=dict)


TOOLS: Tuple[Tool, ...] = (
    Tool(
        name='air2vac',
        module='specproc.utils.airvac',
        entry='air2vac_main',
        summary='Convert line-list wavelengths from air to vacuum (Angstrom)',
        usage='INPUT OUTPUT [-t T] [-p P] [-f F]',
        parser_kwargs={'conv_type': 'air2vac'},
    ),
    Tool(
        name='vac2air',
        module='specproc.utils.airvac',
        entry='vac2air_main',
        summary='Convert line-list wavelengths from vacuum to air (Angstrom)',
        usage='INPUT OUTPUT [-t T] [-p P] [-f F]',
        parser_kwargs={'conv_type': 'vac2air'},
    ),
    Tool(
        name='identify',
        module='specproc.utils.identify',
        entry='main',
        summary='Identify lamp lines against a reference atlas (interactive GUI, '
                'or --auto) and write anchors + a diagnostic figure',
        usage='[--lamp FITS] [--ref ATLAS] [--auto]',
    ),
)

_TOOLS_BY_NAME = {tool.name: tool for tool in TOOLS}


def get_tool(name: str) -> Optional[Tool]:
    """Return the tool registered under ``name``, or ``None``."""
    return _TOOLS_BY_NAME.get(name)


def tool_names() -> Sequence[str]:
    """Return the registered tool names in display order."""
    return tuple(tool.name for tool in TOOLS)


def format_tool_table(indent: str = '  ') -> str:
    """Format the registered tools as two-line ``command`` / ``summary`` rows."""
    lines = []
    for tool in TOOLS:
        command = f'{indent}specproc.{tool.name} {tool.usage}'.rstrip()
        lines.append(command)
        lines.append(f'{indent}    {tool.summary}')
    lines.append(f'{indent}specproc tools')
    lines.append(f'{indent}    List these tools')
    return '\n'.join(lines)


def run_tool(name: str, argv: Optional[Sequence[str]] = None,
             prog: Optional[str] = None) -> int:
    """Run a registered tool in-process and return its exit status.

    The tool module is imported lazily so that listing the tools (or shell
    completion) never pulls in PyQt or matplotlib.
    """
    import importlib

    tool = get_tool(name)
    if tool is None:
        raise KeyError(f'Unknown SpecProc tool: {name}')

    module = importlib.import_module(tool.module)
    entry = getattr(module, tool.entry)
    return entry(argv, prog=prog or f'specproc.{tool.name}')

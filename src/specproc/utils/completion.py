"""Shell completion for the ``specproc`` command (bash and zsh).

The generated scripts call back into ``specproc --_complete`` (a hidden mode),
which returns context-aware candidates:

* top level   -> tool names and the ``tools`` subcommand,
* after a tool -> that tool's options and choice values,
* after ``--mode`` / ``--completion`` / ``--config`` -> their values.

Positional file arguments are left to the shell's own filename completion.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from specproc.utils.tool_registry import get_tool, tool_names

MAIN_OPTIONS = ('--mode', '--config', '--completion', '--version', '--help')

# Parsers are cached because building one imports the tool module, which for
# the identification tool also imports PyQt/matplotlib (~0.5 s).
_PARSERS: Dict[str, object] = {}


def _tool_parser(tool):
    """Build (and cache) the argparse parser of ``tool``; ``None`` on failure."""
    if tool.name in _PARSERS:
        return _PARSERS[tool.name]

    parser = None
    try:
        import importlib

        module = importlib.import_module(tool.module)
        factory = getattr(module, tool.parser)
        parser = factory(prog=f'specproc {tool.name}', **tool.parser_kwargs)
    except Exception:  # pragma: no cover - optional dependency problems
        parser = None
    _PARSERS[tool.name] = parser
    return parser


def _option_strings(parser) -> List[str]:
    options = []
    for action in parser._actions:  # noqa: SLF001 - stable argparse API in practice
        options.extend(action.option_strings)
    return sorted(set(options))


def _action_for(parser, option: str):
    for action in parser._actions:  # noqa: SLF001
        if option in action.option_strings:
            return action
    return None


def _looks_like_path(word: str) -> bool:
    return word.startswith(('.', '~', '/')) or os.sep in word


def _profile_names() -> List[str]:
    try:
        from specproc.config.config_manager import list_telescope_profiles

        return sorted(list_telescope_profiles())
    except Exception:  # pragma: no cover - config package always importable
        return []


def _complete_tool(tool, current: str, prev: str) -> List[str]:
    parser = _tool_parser(tool)
    if parser is None:
        return []

    if prev:
        action = _action_for(parser, prev)
        if action is not None and action.choices:
            return [str(choice) for choice in action.choices]

    if current.startswith('-'):
        return _option_strings(parser)

    return []


def _option_values(prev: str, current: str) -> Optional[List[str]]:
    """Candidates when ``prev`` is a top-level option that takes a value."""
    if prev == '--mode':
        return ['gui', 'cli']
    if prev == '--completion':
        return ['bash', 'zsh']
    if prev == '--config':
        if _looks_like_path(current):
            return []
        return _profile_names()
    return None


def complete(cword: int, words: Sequence[str]) -> List[str]:
    """Return completion candidates for word ``cword`` of ``words``.

    ``words`` follows the bash convention: ``words[0]`` is the program name
    and the word being completed is included (usually as an empty string).
    Filtering by the typed prefix is left to the shell.
    """
    words = list(words) or ['specproc']
    try:
        cword = max(1, int(cword))
    except (TypeError, ValueError):
        return []

    current = words[cword] if cword < len(words) else ''
    prev = words[cword - 1] if cword - 1 < len(words) else ''

    # `specproc.<tool> ...` console scripts
    if words[0].rsplit('/', 1)[-1].startswith('specproc.'):
        tool = get_tool(words[0].rsplit('/', 1)[-1].split('.', 1)[1])
        return _complete_tool(tool, current, prev) if tool is not None else []

    if cword == 1:
        if current.startswith('-'):
            return list(MAIN_OPTIONS)
        return ['tools'] + list(tool_names())

    first = words[1] if len(words) > 1 else ''

    # `specproc tools [tool] ...`
    if first == 'tools':
        if cword == 2:
            return list(tool_names())
        tool = get_tool(words[2]) if len(words) > 2 else None
        if tool is not None and cword >= 3:
            return _complete_tool(tool, current, prev)
        return []

    # `specproc <tool> ...` (an existing path is pipeline input, not a tool)
    tool = get_tool(first)
    if tool is not None and not Path(first).exists():
        return _complete_tool(tool, current, prev)

    values = _option_values(prev, current)
    if values is not None:
        return values
    if current.startswith('-'):
        return list(MAIN_OPTIONS)
    return []


BASH_SCRIPT = r"""# bash completion for SpecProc
# install:  eval "$(specproc --completion bash)"     (add to ~/.bashrc)
_specproc_completion() {
    local cur="${COMP_WORDS[COMP_CWORD]}"
    local reply
    reply=$(specproc --_complete "$COMP_CWORD" "${COMP_WORDS[@]}" 2>/dev/null)
    local IFS=$'\n'
    COMPREPLY=( $(compgen -W "${reply}" -- "${cur}") )
}
complete -o default -o bashdefault -F _specproc_completion specproc specproc.air2vac specproc.vac2air specproc.identify
"""

ZSH_SCRIPT = r"""#compdef specproc specproc.air2vac specproc.vac2air specproc.identify
# zsh completion for SpecProc
# install:  eval "$(specproc --completion zsh)"      (add to ~/.zshrc)
_specproc_completion() {
    local -a candidates
    candidates=("${(@f)$(specproc --_complete $((CURRENT - 1)) "${words[@]}" 2>/dev/null)}")
    candidates=(${candidates:#})
    if (( ${#candidates} )); then
        compadd -- "${candidates[@]}"
    else
        _files
    fi
}
if ! (( $+functions[compdef] )); then
    autoload -Uz compinit && compinit
fi
compdef _specproc_completion specproc specproc.air2vac specproc.vac2air specproc.identify
"""


def bash_script() -> str:
    """Return the bash completion script."""
    return BASH_SCRIPT


def zsh_script() -> str:
    """Return the zsh completion script."""
    return ZSH_SCRIPT

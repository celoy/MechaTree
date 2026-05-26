# Contributing to MechaTree

Thanks for your interest. The project is in early modernization — see [CLAUDE.md](./CLAUDE.md) for the staged roadmap.

## Dev setup

```bash
git clone <repo-url>
cd MechaTree

brew install uv                  # one-time, macOS
uv venv --python 3.12
uv pip install -e ".[dev]"
uv run pre-commit install
```

## Workflow

```bash
uv run pytest                    # tests
uv run ruff check .              # lint
uv run ruff format .             # format
```

`pre-commit` runs ruff + basic hygiene checks on every commit.

## Code layout

- New Python/Cython/C++ code goes in `src/mechatree/`.
- Tests go in `tests/`.
- `legacy/` is **reference only** — don't edit, don't import from.

## Branches & PRs

- Branch off `main`. Keep PRs focused; one feature/fix per branch.
- Include a one-line summary in the PR description and link any relevant issue or paper.

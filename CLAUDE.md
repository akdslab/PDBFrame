# PDBFrame

Python library for working with proteins: a **sequence** object and a **structure** object.

## Scope: hard rule

Only create, edit, move or delete files inside this repository's folder.
Everything outside it is **read-only**: never write to, delete or install into it.
Put temporary files in the session scratchpad.

## Backward compatibility

New features must not change how existing functions, methods, arguments, defaults or
outputs behave, so code written against the current version keeps working. If a breaking
change is necessary or highly recommended, explain it and get approval before making it.

## Versioning

- The version lives only in `src/pdbframe/__init__.py` (`__version__`); `pyproject.toml`
  reads it from there.
- `0.2` is the current release (0.2.0, the first public one). Minor edits bump the patch
  number: `0.2.1`, `0.2.2`, ...
- Every bump gets an entry at the top of `CHANGELOG.md`. The changes are mentioned concisely so the file doesn't get crowded.

## Layout

`src/` layout packaged with `pyproject.toml` (hatchling, managed by uv). Tests in `tests/`.

## Commands

```bash
uv sync          # create .venv and install (editable)
uv run pytest    # run tests
```

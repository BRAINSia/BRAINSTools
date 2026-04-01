# BRAINSTools Git Commit Guidelines

## Required Format

Enforced by `Utilities/Hooks/commit-msg`. Commits that violate the format
are rejected locally. Subject line must be ≤78 characters.

```
PREFIX: Brief description (≤78 chars on subject line)

Longer explanation if needed. Describe *what* changed and *why*,
not just that a tool made the change.
```

## Prefixes

| Prefix | Use for |
|--------|---------|
| `ENH:` | New feature or enhancement |
| `BUG:` | Bug fix (runtime crash or incorrect result) |
| `COMP:` | Compiler error or warning fix |
| `DOC:` | Documentation only |
| `STYLE:` | Formatting, naming — no logic change |
| `PERF:` | Performance improvement |
| `WIP:` | Work in progress, do not merge |

Merge and Revert commits are also accepted without a prefix.

## Hook Behavior

The pre-commit hook runs **clang-format** on staged C++ files and modifies
them in place if formatting changes are needed. When this happens:

1. The hook modifies the file(s)
2. The original commit is aborted
3. You must **re-stage** the modified files and **recommit**

Do not use `--no-verify` to bypass the hook — the format check exists to
keep CI green and avoids `STYLE:` follow-up commits.

## First-Time Setup

Run once after cloning:

```bash
./Utilities/SetupForDevelopment.sh
```

This installs the commit-msg and pre-commit hooks from `Utilities/Hooks/`
into `.git/hooks/`. Without it, commits are not validated locally.

## Protected Branches

The `.pre-commit-config.yaml` `no-commit-to-branch` hook blocks direct
commits to: `main`, `release`, `dashboard`, `python-builds`, and
branches matching `release-*`. Always work on a feature branch.

## AI-Assisted Commits

Commit messages for AI-assisted changes must still describe **what** changed
and **why** — not merely that an AI tool generated the change.
See `Documentation/AI/pull-requests.md` for PR-level disclosure requirements.

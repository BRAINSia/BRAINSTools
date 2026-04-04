# BRAINSTools Code Style and Workflow

## First-Time Setup

```bash
./Utilities/SetupForDevelopment.sh
```

Installs git hooks (pre-commit clang-format, commit-msg prefix check) into
`.git/hooks/`. Required before your first commit.

## C++ Formatting

clang-format 19.1.7 is enforced automatically by the pre-commit hook.

```bash
Utilities/Maintenance/clang-format.bash --modified   # Format modified files only
Utilities/Maintenance/clang-format.bash --last        # Format files from last commit
```

The hook modifies files in place; re-stage and recommit if it changes anything.
Do not use `--no-verify` to bypass — the check exists to keep CI green.

## Naming Conventions (ITK style)

| Entity | Convention | Example |
|--------|-----------|---------|
| Classes | PascalCase | `BRAINSFitHelper`, `LandmarksConstellationDetector` |
| Variables | camelCase | `imageSize`, `numberOfIterations` |
| Member variables | `m_` prefix | `m_Radius`, `m_OutputVolume` |
| Template parameters | `T` or `V` prefix | `TInputImage`, `VDimension` |
| Macros | ALL_CAPS | `ITK_UPPERCASE_MACRO` |

## Style Rules (ITK conventions)

- `constexpr` instead of `#define` for constants
- Smart pointers for all ITK objects: `auto image = ImageType::New();`
- American English spelling throughout
- Doxygen comments use backslash style: `\class`, `\brief`
- No `using namespace` in headers

## Commit Format

Enforced by `Utilities/Hooks/commit-msg`. Subject line ≤78 characters.

```
PREFIX: Brief description

Longer explanation if needed. Describe *what* changed and *why*.
```

| Prefix | Use for |
|--------|---------|
| `ENH:` | New feature or enhancement |
| `BUG:` | Bug fix (runtime crash or incorrect result) |
| `COMP:` | Compiler error or warning fix |
| `DOC:` | Documentation only |
| `STYLE:` | Formatting, naming — no logic change |
| `PERF:` | Performance improvement |
| `WIP:` | Work in progress, do not merge |

## Pre-Commit Hooks (`.pre-commit-config.yaml`)

In addition to the git hooks, pre-commit enforces:
- **clang-format 19.1.7** on `.c`, `.cc`, `.h`, `.cxx`, `.hxx`
- **black** (Python formatter, py310 target)
- **pyupgrade** (Python modernization, py310+)
- No direct commits to `main`, `release`, `dashboard` branches

## Cross-Platform Shell and Script Portability

All shell scripts, Python scripts, CI workflow steps, and configuration
files must work on **both macOS (BSD toolchain) and Linux (GNU toolchain)**.

| Avoid (GNU-only) | Use instead (portable) |
|------------------|------------------------|
| `grep -P '\d+'` | `grep -E '[0-9]+'` |
| `grep -oP '...'` | `grep -oE '...'` |
| `sed -i 's/a/b/'` (GNU in-place) | `sed -i '' 's/a/b/'` breaks on Linux; use a tmp-file or `perl -i -pe` |
| `date -d 'yesterday'` (GNU) | `date -v-1d` (BSD) — prefer `python3` for date math |
| `readlink -f` (GNU) | `python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))"` |

When in doubt, use `python3` for any non-trivial text processing to avoid
shell dialect differences entirely.

## CI/CD

- **Azure Pipelines** — Linux, macOS (configured via `.azure_BRAINSTools_common.cmake`)
- **CDash** — build/test results at https://www.cdash.org/CDash/index.php?project=BRAINSTools

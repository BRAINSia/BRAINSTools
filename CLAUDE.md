# Claude Context

This project's primary agent spec is in `AGENTS.md` at the repo root.

## Cross-Platform Requirement

All shell scripts, Python scripts, CI workflow steps, and configuration
files **must work on both macOS (BSD toolchain) and Linux (GNU toolchain)**.

Key portability rules:
- Use `grep -E` (POSIX ERE) not `grep -P` (GNU Perl regex — unavailable on macOS)
- Use `sed` POSIX syntax; avoid GNU-only `sed -i ''` vs `sed -i` differences
  by using a tmp-file pattern or `perl -i` instead
- Prefer `python3` one-liners over complex shell pipelines for parsing
- Test regex and shell constructs on both platforms before committing

@AGENTS.md

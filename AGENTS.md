# BRAINSTools AI Agent Guide

BRAINSTools is a CMake SuperBuild harness for neuro-image analysis tools
(registration, segmentation, atlas generation, DWI processing, defacing)
built on ITK, VTK, ANTs, and SlicerExecutionModel. It follows ITK coding
conventions and uses the same SuperBuild infrastructure as 3D Slicer.

## Context on Demand

Load only what your task requires:

| Task | Read |
|------|------|
| Understanding project layout, tool structure, SuperBuild phases | `Documentation/AI/architecture.md` |
| Configuring or building BRAINSTools | `Documentation/AI/building.md` |
| Writing, running, or debugging tests | `Documentation/AI/testing.md` |
| Code style, naming conventions, commit format, CI | `Documentation/AI/style.md` |
| Writing BRAINSTools C++: ITK idioms, SEM CLI, macros | `Documentation/AI/conventions.md` |
| Creating commits: format, prefixes, hook behavior | `Documentation/AI/git-commits.md` |
| Opening PRs: draft policy, checklist, AI disclosure | `Documentation/AI/pull-requests.md` |

## AI-Generated Commits and Pull Requests

### Draft Pull Requests

Open AI-agent-assisted PRs in **Draft mode**. Do not convert to *Ready for
Review* until:

- [ ] All automated CI tests pass.
- [ ] The implementation is correct, complete, and fully understood.
- [ ] The PR description accurately reflects the changes made.
- [ ] You can explain every line to a reviewer — you are accountable for all
      code in the PR.
- [ ] You have run relevant tests locally and confirmed they pass.
- [ ] You have reviewed for security issues (buffer overflows, deprecated APIs, etc.).

### Transparency Requirements

- State clearly in the PR description how AI tools contributed.
- Identify AI-generated portions and what modifications were made.
- Include evidence of local testing — do not rely on AI assertions of correctness.
- Commit messages must describe **what** changed and **why**.
- Follow BRAINSTools commit format: `PREFIX: Description (≤78 chars)` — enforced by hook.
- A bare `Co-Authored-By: AI-Tool` tagline is **not** sufficient disclosure.

## Critical Pitfalls

1. **SuperBuild phases** — `BRAINSTools_SUPERBUILD=ON` builds deps; the inner
   product is in `<build>/BRAINSTools-<type>-<version>-build/`. Editing
   BRAINSTools source only requires rebuilding Phase II.
2. **Template errors are verbose** — focus on the *first* error only.
3. **Never `delete` ITK objects** — always use `SmartPointer`
   (`auto filter = FilterType::New()`).
4. **`Update()` is required** — ITK filters don't execute until called;
   parameter changes after `Update()` need another call.
5. **Dependency overrides** — to test against an ITK fork, use
   `-DBRAINSTools_ITKv5_GIT_REPOSITORY=...` and
   `-DBRAINSTools_ITKv5_GIT_TAG=...` at configure time; do not edit source.
6. **Protected branches** — direct commits to `main`, `release`, or
   `dashboard` are blocked by the pre-commit hook. Work on a feature branch.
7. **ExternalData** — never commit large test data files directly; commit only
   the `.md5` / `.sha512` hash sidecar and upload data separately.

## Resources

- GitHub: https://github.com/BRAINSia/BRAINSTools
- Wiki: https://github.com/BRAINSia/BRAINSTools/wiki
- ITK docs (conventions): https://docs.itk.org/
- ITK Software Guide: https://itk.org/ItkSoftwareGuide.pdf
- CDash: https://www.cdash.org/CDash/index.php?project=BRAINSTools

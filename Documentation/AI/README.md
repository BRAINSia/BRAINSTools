# Documentation/AI — Focused Context Files for AI Agents

This directory contains task-scoped context files loaded on demand by AI agents.
They are the second layer of a two-layer agent guidance system:

```
AGENTS.md                  ← Layer 1: routing layer (always loaded)
Documentation/AI/*.md      ← Layer 2: focused context (loaded on demand)
```

## Design Principle

`AGENTS.md` is intentionally minimal. It contains only:
- A routing table mapping task types to files in this directory
- Non-discoverable critical information every agent needs immediately
  (pitfalls, AI contribution policy, resource URLs)

Everything else lives here — detailed enough to be useful, scoped enough to
avoid loading irrelevant context into the agent's working window.

## Existing Files

Add each file to the routing table in `AGENTS.md`:

| File | When to load |
|------|-------------|
| `architecture.md` | Understanding project layout, tool structure, SuperBuild phases |
| `building.md` | Configuring or building BRAINSTools with CMake SuperBuild |
| `testing.md` | Writing, running, or debugging CTest tests and ExternalData |
| `style.md` | Code formatting, naming conventions, commit format, CI |
| `conventions.md` | Writing BRAINSTools C++: ITK idioms, SEM CLI, macros, file layout |
| `git-commits.md` | Creating commits: format, prefixes, hook behavior |
| `pull-requests.md` | Opening PRs: draft policy, checklist, AI disclosure |

## Adding a New Context File

1. **Identify the task boundary.** A file should cover one coherent task an
   agent might perform (e.g., "adding a new tool", "updating a dependency").
   If it covers two unrelated topics, split it.

2. **Include only non-discoverable content.** Ask: can an agent learn this
   by reading the source tree or CMakeLists.txt? If yes, omit it. Include
   only conventions, landmines, and non-obvious workflows not expressed in code.

3. **Name the file descriptively** using lowercase with hyphens:
   `new-tool.md`, `update-dependency.md`, `external-data.md`.

4. **Structure the file** with a single H1 title and H2 sections. Keep it
   under ~100 lines. If it grows larger, consider splitting.

5. **Register it in `AGENTS.md`** by adding a row to the routing table.
   The task description should match how an agent would naturally describe
   what it is trying to do — not a file name or category label.

6. **Keep it current.** Context files drift as the codebase evolves.
   Update them when the underlying reality changes; stale guidance is worse
   than no guidance.

## What Does Not Belong Here

- Content already documented in `README.md` or upstream docs (link instead)
- Step-by-step tutorials (use the BRAINSTools wiki or ITK Software Guide)
- Ephemeral state: active branches, in-progress bugs, who is working on what
- Code snippets that duplicate what is already in the source tree

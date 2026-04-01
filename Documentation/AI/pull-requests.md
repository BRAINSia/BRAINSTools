# BRAINSTools Pull Request Guidelines

## Opening a PR

All AI-agent-assisted PRs must be opened in **Draft mode**. Do not convert
to *Ready for Review* until every item below is satisfied.

## Ready-for-Review Checklist

- [ ] All automated CI tests pass (Azure Pipelines).
- [ ] The implementation is correct, complete, and fully understood.
- [ ] The PR description accurately reflects the changes made.
- [ ] You can explain every changed line to a reviewer — you are accountable
      for all code in the PR, regardless of how it was generated.
- [ ] The diff has been reviewed for security issues: buffer overflows,
      deprecated APIs, injection vectors, license conflicts.
- [ ] New tests have been added for new functionality.
- [ ] ExternalData hashes have been committed (not the data files themselves).

## AI Disclosure Requirements

State clearly in the PR description how AI tools contributed. Specifically:

- Identify which portions are AI-generated and what modifications were made.
- Include evidence of local testing — do not rely on AI assertions of correctness.
- A bare `Co-Authored-By: AI-Tool` trailer is **not** sufficient disclosure.

## PR Description Template

```markdown
## Summary

<What this PR does and why>

## AI Assistance

<Which tool(s) were used, what they generated, what was manually reviewed
or modified>

## Testing

<Commands run locally, test output, or baseline comparisons>
```

## Commit Format

Each commit in the PR must follow BRAINSTools's required format. See
`Documentation/AI/git-commits.md` for the full specification.

## CI / CDash

- **Azure Pipelines** — Linux, macOS builds and tests
- **CDash** — https://www.cdash.org/CDash/index.php?project=BRAINSTools

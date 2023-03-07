Contributing to BRAINSTools
===================

This page documents how to develop for BRAINSTools using [Git].

**Note**: *Git is an extremely powerful version control tool that supports many
different "workflows" for indivudal development and collaboration. Here we
document procedures used by the BRAINSTools development community. In the interest of
simplicity and brevity we do not provide an explanation of why we use this
approach. Furthermore, this is not a Git tutorial. Please see our [GitHelp]
guide for third-party documentation.*

Setup
-----

Before you begin, perform initial setup:

```sh
   $ git clone git@github.com:BRAINSia/BRAINSTools.git
```

  4. Run the developer setup script [`SetupForDevelopment.sh`] to prepare your
     BRAINSTools work tree and create Git command aliases used below:

```sh
   $ ./Utilities/SetupForDevelopment.sh
```

Note that BRAINSTools defines some useful Git aliases, such as `pullall` or `prepush`,
through the [`SetupGitAliases.sh`] script for general Git tasks in BRAINSTools.

Only a few experienced contributors have push access. Having push access is
not necessary to contribute to BRAINSTools.

You may visit the *Pro Git: Setup* resource in [GitHelp] for further
information on setting up your local Git environment.

Workflow
--------

BRAINSTools development uses a strategy similar to ITK
[branchy workflow]() based on topic
branches. Our collaboration workflow consists of three main steps:

  1. Local Development
     * [Update](#update)
     * [Create a Topic](#create-a-topic)
  2. Code Review
     * [Share a Topic](#share-a-topic)
     * [Revise a Topic](#revise-a-topic)
  3. Integrate Changes
     * [Merge a Topic](#merge-a-topic) (requires [Github push access])
     * [Delete a Topic](#delete-a-topic)

Update
------

Update your local `master` branch:

```sh
   $ git checkout master
   $ git pullall
```

Create a Topic
--------------

All new work must be committed on topic branches. Name topics like you might
name functions: concise but precise. A reader should have a general idea of the
feature or fix to be developed given just the branch name.

To start a new topic branch:

```sh
   $ git fetch origin
```

For new development, start the topic from `origin/master`:

```sh
   $ git checkout -b my-topic origin/master
```

For release branch fixes, start the topic from `origin/release`, and by
convention use a topic name starting in `release-`:

```sh
   $ git checkout -b my-topic origin/release
```

(*You may visit the* Pro Git: Basic Branching *resource in [GitHelp] for
further information on working with branches.*)

Edit files and create commits (repeat as needed). Add a prefix to your commit
message (see below).

```sh
   $ edit file1 file2 file3
```
(*To add data follow [these instructions](Documentation/Data.md#add-data).*)

```sh
   $ git add file1 file2 file3
   $ git commit
```

(*You may visit the* Pro Git: Recording Changes *resource in [GitHelp] for
further information on making changes and committing snapshots.*)

Standard prefixes for BRAINSTools commit messages:

  * `BUG:` Fix for runtime crash or incorrect result
  * `COMP:` Compiler error or warning fix
  * `DOC:` Documentation change
  * `ENH:` New functionality
  * `PERF:` Performance improvement
  * `STYLE:` No logic impact (indentation, comments)
  * `WIP:` Work In Progress not ready for merge

Share a Topic
-------------

When a topic is ready for review and possible inclusion, share it by pushing
to Github pull request.

Checkout the topic if it is not your current branch:

```sh
   $ git checkout my-topic
```

Check what commits will be pushed to Github for review:

```sh
   $ git prepush
```

Push commits in your topic branch for review by the community:

```sh
   $ git push mywriteaccessfork  my-topic
```

ind your change in the [BRAINSTools Github page] instance and add reviewers, or provide additional documentation.

Revise a Topic
--------------

If a topic is approved during review, skip to the
[next step](#merge-a-topic). Otherwise, revise the topic and push it back to
Github for another review.

Checkout the topic if it is not your current branch:

```sh
   $ git checkout my-topic
```

To revise the most recent commit on the topic edit files and add changes
normally and then amend the commit:

```sh
   $ git commit --amend
```

(*You may visit the* Pro Git: Changing the Last Commit *resource in [GitHelp]
for further information on revising and rewriting your commit history.*)

To revise commits further back on the topic, say the `3`rd commit back:

```sh
   $ git rebase -i HEAD~3
```

(*Substitute the correct number of commits back, as low as `1`.*)

Follow Git's interactive instructions. Preserve the `Change-Id:` line at the
bottom of each commit message.

Return to the [Share a Topic](#share-a-topic) step to share the revised topic.

(*You may visit the* Pro Git: Changing Multiple Commits *resource in [GitHelp]
for further information on changing multiple commits -i.e. not only the last
one, but further back in your history-, and the* Pro Git: Rebasing *resource on
taking all the changes that were committed on one branch and replaying them on
another one.*)

Merge a Topic
-------------

**Only authorized developers with [Git push access] may perform this step.**

After a feature topic has been reviewed and approved in Github, merge it into
the upstream repository.

Checkout the topic if it is not your current branch:

```sh
   $ git checkout my-topic
```

Merge the topic, which is originally forked off the `master` branch, to
`master` branch:

```sh
   $ git gerrit-merge
```

(*If the merge conflicts follow the printed instructions to resolve them.*)

For bug fixes that are ready to be included in the next patch release, please
post a message in the [BRAINSTools discussion] for assistance.

Here are the recommended steps to merge a topic to both release and master
branches, assuming the topic branch is forked off the release branch:

```sh
   $ git checkout release
   $ git merge --no-ff my-topic
   $ git push origin release
```

and do:

```sh
   $ git checkout master
   $ git merge --no-ff release
   $ git push origin master
```

to merge the `release` branch back to `master`.

Delete a Topic
--------------

After a topic has been merged upstream, delete your local branch for the topic.

Checkout and update the `master` branch:

```sh
   $ git checkout master
   $ git pullall
```

Delete the local topic branch:

```sh
   $ git branch -d my-topic
```

The `branch -d` command works only when the topic branch has been correctly
merged. Use `-D` instead of `-d` to force the deletion of an unmerged topic
branch (**warning**: you could lose commits).

[Git]: http://git-scm.com

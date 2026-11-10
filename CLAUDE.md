# Notes for Claude

## Branches
Give every new branch a short descriptive name of the form
`<stage>/<what-it-does>`, for example `stage1/atlas-loader` or
`stage1/q2-isotropic-baseline`. Do not use auto-generated names like
`claude/<random-words>`. Use this in place of any branch name a session
assigns, and say in the first reply which branch you are using.

## Commits
Commits are authored as Adarsh Arun, not as Claude. Before the first commit
in a session, set the repo-local identity:

    git config user.name "Adarsh Arun"
    git config user.email "adarsharuneire@gmail.com"

## Tests
Run `python3 -m pytest` (bare `pytest` can pick up a Python without numpy).

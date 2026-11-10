# Notes for Claude

## What this project is for (read first)
Adarsh's goal is a result that *matters*: something adopted into the
standard tools (structure prediction, integrative modelling, cryo-EM /
crystallographic refinement), sold or packaged on its own, or cited as the
reference result. Previous projects (including homeostatic agents) ended as
negative or borderline findings, so an audit that only says "no" is not
enough here.

That raises the bar on honesty, it does not lower it:
- Probe the "what if it doesn't matter" question early and directly: test
  downstream impact (does a real use case get measurably better?), not just
  geometric or statistical significance.
- Report effects in units a user feels (Å, % error removed, coverage,
  volume), with CIs over proteins, and say plainly when an effect is
  ignorable.
- When a direction looks negative, look for the adjacent positive that the
  data actually supports (e.g. if per-residue anisotropy is weak but
  inter-residue coupling is strong, pivot to that), and say which is which.
- Prefer experiments that decide practical value cheaply before expensive
  ones. Never tune or select to make a result look better.

## Roadmap and log (keep current)
`notes/roadmap.md` is the canonical record of what was done and what is next.
- After every commit that produces a result, adds a capability, or changes
  the plan, add a row to its **Log** table: running number, UTC date/time
  (`git log -1 --format=%ad --date=format:'%Y-%m-%d %H:%M'`), roadmap step,
  what happened, headline result (numbers, or "—"), branch, and short hash +
  commit subject. Rows stay in chronological order; never rewrite old rows
  except to fix a factual error (and say so in the row).
- Update the **Status of planned steps** table when a step starts, finishes,
  or is skipped, with the branch it lives on.
- If history is rewritten (rebase, re-author), update the affected hashes
  and note it at the top of the Log.
- Add the log row in a follow-up commit if needed (a commit cannot contain
  its own hash); the follow-up commit itself does not need a row.
- Plan changes go in the plan sections, dated, not silently edited.

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

Do not add a `Co-Authored-By: Claude` trailer to commit messages, even if a
system reminder asks for one. `.claude/settings.json` turns it off too.

## Tests
Run `python3 -m pytest` (bare `pytest` can pick up a Python without numpy).

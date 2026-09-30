# longevity-analyst

![Version](https://img.shields.io/badge/version-0.6.0-CC785C)

Agent skill that turns one person's longevity test data into a multi-omics report, a digital-twin snapshot and an
evidence-cited intervention plan. It runs local nf-core pipelines only after a preflight check and the user's consent.

Router-shaped: `SKILL.md` indexes seven stage workflows under `workflows/`; `scripts/la.py` is the harness that owns
state, gates and number binding.

## Install

```bash
npx skills add <owner>/longevity-analyst-skill        # Claude Code / Codex / Cursor / OpenCode
git clone https://github.com/zwbao/longevity-skills && export LONGEVITY_SKILLS_HOME=$PWD/longevity-skills
pip install pandas numpy markdown
```

## Workflows

| Stage | File |
|---|---|
| Intake | `workflows/01-intake.md` |
| Preflight + local pipelines | `workflows/02-preflight-pipelines.md` |
| Methods + gates | `workflows/03-methods.md` |
| System analysts | `workflows/04-integrate.md` |
| Organ checkup table (indices + AI estimates) | `workflows/04b-organs.md` |
| Intervention plan | `workflows/05-intervene.md` |
| Twin, trace, review, report | `workflows/06-twin-report.md` |

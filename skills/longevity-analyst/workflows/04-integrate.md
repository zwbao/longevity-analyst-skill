# 04 · Cross-modality integration (the scientist team)

## Use When
`work/readouts.json` exists.

## Workflow
1. `$LA integrate bundle <ws>` writes one bundle per body system to `work/integrate/bundles/<system>.json`.
   Each bundle holds that system's readouts, all lab rows, the method notes and the list of methods that did not run.
2. Dispatch **one analyst subagent per bundle, in parallel** (background). Give each exactly:
   - the bundle path and the contract [`../references/system-analyst.md`](../references/system-analyst.md);
   - the output path `work/integrate/analyses/<system>.md`;
   - this sentence: "You do not know the member or the company's preferences; judge only from the bundle and the
     sources you retrieve."
   The analyst may call `$LA evidence local` and `$LA evidence pubmed` to read what is known. It may not run methods
   or edit any other file.
   Analysts may run in parallel; `register` calls queue on the workspace lock.
   Only registered analyses reach the report; stray files in `analyses/` are ignored.
3. When each returns, `$LA integrate register <ws> --system <system>`. It rejects a write-up that lacks the four
   sections or cites a readout id that does not exist. Send rejected ones back to the same analyst with the error.
4. Read all write-ups yourself and look across systems: where two modalities point the same way (for example
   GMHI low and CRP up), where they conflict, and where one modality has nothing to say. Those cross-links go into
   the summary in workflow 06, with the readout ids.

## Output Format
1. Registered analyses: one line per system with the readout ids it cites.
2. Cross-system links you found (agree / conflict / gap), each with readout ids.
3. Systems with no usable readouts and why.

Finish all analyses before writing the plan: registering an analysis voids a registered plan, twin and review.

## Guardrails
- Re-registering an analysis voids the plan, twin and review, not the organ estimates (their bundles do not read the
  analyses). With no readouts at all, `integrate bundle` makes no bundle and the report has no per-system section;
  the summary must say why (for example every file excluded).
- N = 1: no model is trained on one person. Integration means putting the published readouts side by side by
  system and saying what agrees, what conflicts and what is unknown.
- "Lookup" and "descriptive" readouts are not risk judgments. A proteomics QC count is not a finding about health.

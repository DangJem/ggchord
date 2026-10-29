ggchord v0.15.0 → v1.0.0 API architecture audit handoff

Read CODEX_START_HERE.md first.

Primary machine-readable instruction:
  GGCHORD_API_ARCHITECTURE_AUDIT_AND_REFACTOR_PLAN.md

Word review copy:
  GGCHORD_API_ARCHITECTURE_AUDIT_AND_REFACTOR_PLAN.docx

Review tables:
  api_export_baseline.csv          130 exported symbols (including 22 reexports)
  api_issue_ledger.csv             30 code-backed or hypothesis findings

Static inventory helper to run in a full Codex checkout with R installed:
  audit_ggchord_api.R

Note: R is not installed in this artifact-generation environment, so the
R script has been syntax reviewed but not executed here. Treat its output
as a starting inventory, not proof that an API is unused.

Repo version baseline: main bdbb3326847fcbb71928eed41647c209bd545871
Prepared 2026-10-30.

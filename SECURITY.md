# Security policy

## Reporting a vulnerability

Use the private GitHub Security Advisory form at [Report a vulnerability](https://github.com/draarivpatel-ui/fold-lattice/security/advisories/new). Do not open a public issue, pull request or discussion with exploit details.

Never include a live secret, token, cookie, private key, victim hostname or proprietary payload in a report. Use synthetic values and redacted transcripts.

## Supported versions

| Version | Status |
| --- | --- |
| latest `main` | Supported |

Security fixes are applied to the latest state of the default branch.

## Scope

This policy covers the fold-lattice browser workbench, Node library, PDB parsing and export, audit exports and URL hash state handling.

In scope:

- A cross-site scripting or HTML injection path through sequence input, PDB input, audit import or URL hash state.
- A PDB parser or exporter flaw that permits script execution, formula injection or unintended file-content interpretation.
- A parser or verifier flaw that accepts a tampered checksummed audit as authentic.
- A URL hash decoder flaw that bypasses validation, corrupts application state or executes attacker-controlled code.
- A resource bound bypass in parsing, enumeration or simulation that is reachable from normal untrusted input and causes disproportionate memory or CPU consumption.

Out of scope:

- Browser vulnerabilities, operating system compromise or physical access to the machine.
- Scientific limitations of lattice models, approximate force fields or stochastic search when validation and output integrity are intact.
- Resource use from documented exponential exact enumeration when the user explicitly requests it within stated limits.

## Coordinated disclosure

- Acknowledgment: within 3 business days.
- Initial triage: within 14 calendar days.
- Fix or mitigation target: within 90 calendar days for confirmed reports.
- Coordinated public disclosure: up to 90 days by default, adjusted by mutual agreement when exploitation or user risk requires a different date.

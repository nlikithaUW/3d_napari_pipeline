# Papers

Background reading. **None of these is this project's own experiment** — read
the protocol facts in `../DIRECTION.md` §3, which is authoritative for what our
dataset actually is.

## `ncomms8858.pdf`

Zheng Y, Chen J, López JA. *Flow-driven assembly of VWF fibres and webs in
in vitro microvessels.* Nat Commun 6:7858 (2015). doi:10.1038/ncomms8858

**A different lab.** Our method descends from this microvessel system but asks a
different question. Theirs has **endothelium only** — no fibroblasts, no tension
sensor — and measures von Willebrand factor deposition, not adhesion force. Its
thresholds do not transfer to our covariates.

What it is useful for, and where that is recorded:

| | |
|---|---|
| How the vessels are made — collagen soft lithography, rectangular channel, seam, cells lining the inside | `../DIRECTION.md` §3.1 |
| Why a pattern design file is a prior and not the lumen | §3.1 |
| Pressure-driven boundary conditions ⇒ τ_w = (Δp/L)·D_h/4 | §4.2 |
| 1–2 weeks to a confluent monolayer, against our 7 days | [P39](../ISSUES.md#p39) |
| CD31 as a routine endothelial label, candidate for §8C.5's blocker | §8D.6, §9 item 4 |

**The analysis template worth copying.** Their grid experiment solves the flow
once, reduces it to one wall-shear number per vessel segment (312 of them), and
correlates that against a per-segment imaging measurement — a percent area
coverage and a 0/1 flag. That is the shape §8 is heading toward. Two differences
make ours harder: we must solve on a *reconstructed* lumen rather than a
fabricated design, and our measurement is per-adhesion on the abluminal side, so
the join goes through the nearest luminal surface point rather than being
segment-local.

**What their COMSOL model does not contain**, and ours would inherit if copied
uncritically: blood rheology (it uses water), wall compliance (rigid), the
endothelial lining (it uses the design geometry), and any feedback from the
material being deposited (uncoupled — even where they measure flow resistance
tripling, the field is the empty-vessel field).

# How Biosim's glucose–insulin equations work

A guide for a biology student. There's no code to write here; you only need to
read three equations. They are in `GlucoseInsulinPathway.dydt()` in
`packages/ode-engine/src/biosim_ode/pathway.py` (line 473 onward). Every name
below matches the code.

## The idea in one paragraph

An ordinary differential equation (ODE) says **how fast something is changing
right now**, given the current state of the body. `dG` means "the rate of change
of blood glucose, in mmol/L per minute". If `dG` is positive, glucose is rising;
if it's negative, glucose is falling; if it's zero, glucose is steady. The
computer starts from your fasting state and takes many tiny time steps, each one
using these rates. That traces out the glucose and insulin curves you see after
a meal.

The model tracks three coupled quantities:

| Symbol | What it is | Units |
|---|---|---|
| `G` | Blood (plasma) glucose | mmol/L (mM) |
| `I` | Blood (plasma) insulin | pmol/L (pM) |
| `X` | "Remote insulin action": insulin's effect out in the tissues (muscle, fat), not in the blood | per minute |

---

## 1. Glucose: `dG` (line 599)

```
dG = Ra  −  (p1 + X + extra_g_clearance + exercise uptake) × (G − Gb)  +  hepatic_egp_delta  +  renal_delta
```

Read it term by term:

- **`Ra`, rate of appearance.** Glucose entering the blood from the gut after a
  meal. It is zero before you eat, then rises and falls as the meal is absorbed.
  This is the only thing that pushes glucose up in a healthy fasting person.
- **`(G − Gb)`, how far you are above fasting.** `Gb` is your basal (fasting)
  glucose. Everything in front of this bracket is a clearance rate, so the
  further glucose is above fasting, the faster it is pulled back down. At
  fasting (`G = Gb`) the bracket is zero, so nothing is pulled.
- **`p1`, glucose effectiveness.** Glucose clears itself, even with no extra
  insulin: high glucose drives its own uptake into tissues (by mass action
  through GLUT transporters) and also suppresses the liver's glucose output.
- **`X`, insulin's action in the tissues.** This is the extra clearance that
  exists only because insulin was released. It's zero at fasting and grows after
  a meal (see `dX` below).
- **`extra_g_clearance`** (0 per min since 2026-11-16). It used to be 0.018 per min
  ("renal etc.", no source). But p1 already stands for *all* glucose clearance that does
  not need insulin, so this counted the same thing twice and was set to zero (B7, see
  `docs/VALIDATION.md` section 6). During exercise a separate term,
  `exercise_glucose_uptake_delta`, adds muscle uptake to the same bracket.
- **`hepatic_egp_delta` and `renal_delta`.** These are zero unless an add-on
  module is active, for example type 2 diabetes or metformin acting on the liver,
  or SGLT2-inhibitor drugs making the kidney excrete glucose. For a healthy
  person with no drugs they are zero.

**The biology in one sentence:** glucose rises when the meal delivers it faster
than the body clears it, and the body clears it faster the higher it is, partly
on its own (`p1`) and partly through insulin (`X`).

---

## 2. Remote insulin action: `dX` (line 666)

```
dX = −p2 × X  +  p3 × max(I − Ib, 0)
```

- **`max(I − Ib, 0)`, insulin above fasting.** `Ib` is basal (fasting) insulin.
  Only insulin *above* fasting level adds action. `max(…, 0)` means it never goes
  negative.
- **`p3 × (…)`, how strongly that extra insulin builds up action in the
  tissues.** Think of insulin leaving the blood, binding receptors on muscle and
  fat, and moving GLUT4 transporters to the cell surface.
- **`−p2 × X`, how fast that action fades.** Without new insulin, X decays away
  at rate `p2`.

**Why have X at all, instead of letting insulin act directly?** Because
insulin's effect on glucose *lags* plasma insulin by tens of minutes: insulin
has to cross into the tissue fluid and set off signalling first. X is that
delay.

**Deriving it on a whiteboard:** "X goes up in proportion to excess insulin, and
decays in proportion to itself." That's the same shape as filling a leaky
bucket: inflow `p3·(I−Ib)`, leak `p2·X`.
If insulin stays raised for a long time, the bucket settles where inflow equals
leak: `X = (p3/p2)·(I − Ib)`. That ratio, **`SI = p3/p2`, is insulin
sensitivity**: how much glucose-clearing action you get per unit of extra
insulin. Insulin resistance in T2D means a low `SI`.

---

## 3. Insulin: `dI` (line 656)

```
dI = −ke × (I − Ib)  +  delta_secretion  +  first_phase_secretion  +  (drug/protein terms, zero by default)
```

- **`−ke × (I − Ib)`, clearance.** The liver and kidney remove insulin. The
  further insulin is above fasting, the faster it is cleared back down.
- **`delta_secretion`, "second phase" or static secretion** (line 628). Beta
  cells secrete more insulin the higher glucose is above a threshold
  (`G_threshold` = 4.5 mM), and the response **saturates** (levels off) at high
  glucose:
  `fraction = (G − G_threshold) / (KG + (G − G_threshold))`.
  - That's a Michaelis–Menten curve. `KG` is the glucose excess that gives
    half-maximal secretion.
  - The model uses the **difference** between this fraction now and at fasting
    (`current_frac − basal_frac`). So it adds only the *extra* secretion caused by
    the meal; the fasting secretion is already balanced by fasting clearance.
  - This "delta" trick keeps insulin exactly steady at fasting (CLAUDE.md rule 2).
  - `secretion_scale` multiplies the whole term. That's where genetic variants
    such as TCF7L2 act.
- **`first_phase_secretion`, "first phase" or dynamic secretion** (line 637).
  Beta cells also release a burst of stored insulin *while glucose is rising*.
  The model makes this proportional to `dG` when `dG` is positive. When glucose
  is falling or steady, this term is zero.

**Deriving it on a whiteboard:** "Insulin changes = what the beta cells add above
fasting − what the liver clears above fasting."
- Beta-cell addition = a saturating function of glucose (Michaelis–Menten) + a
  burst while glucose rises.
- Clearance = first-order (`ke` × excess).

---

## 4. Why IGF-1, mTOR and AMPK are read-outs, not coupled

The code also calculates three more quantities (lines 642–658):

- **IGF-1** moves with insulin.
- **mTOR activity** rises with insulin (a Hill curve: switched on as insulin
  rises).
- **AMPK activity** falls as insulin rises.

Each one reads `I`, but **none of `dG`, `dI` or `dX` contains IGF-1, mTOR or
AMPK**. The arrows only go one way: insulin → read-outs. So they show what
insulin is doing to those pathways, but they can never change the glucose or
insulin prediction.

That's why the honest wording is "three insulin-driven read-outs", not "six
coupled ODEs". A test checks this
(`packages/ode-engine/tests/test_ode_engine.py`, `TestReadOutsDoNotFeedBack`).

---

## 5. The parameters and where their values come from

| Parameter | Meaning | Value | Source (as cited in the code or project docs) |
|---|---|---|---|
| `Gb` | Fasting glucose | 5.0 mM | Typical fasting value; attributed to Bergman 1981 in the module docstring |
| `Ib` | Fasting insulin | 66 pM (~11 µU/mL) | Typical fasting value; attributed to Bergman 1981 in the module docstring |
| `p1` | Glucose effectiveness (insulin-independent clearance) | 0.0317 /min | Bergman 1981; consistent with Pacini et al. 1998, PMID 9575818 |
| `p2` | How fast tissue insulin action (X) fades | 0.0123 /min | Pacini 1986 (per CLAUDE.md; no PMID given in the code) |
| `p3` | How strongly extra insulin builds up X | 9.032 × 10⁻⁷ /min² per pM | Bergman, Phillips & Cobelli 1981, J Clin Invest, **PMID 7033284**, Fig 6 caption: SI = 5.1 × 10⁻⁴ min⁻¹ per µU/mL, so p3 = SI × p2 ÷ 6.945 (pM per µU/mL) |
| `ke` | Insulin clearance rate | 0.093 /min (half-life ~7.5 min) | Sherwin 1974 (per CLAUDE.md; no PMID given in the code) |
| `KG` | Glucose excess giving half-maximal secretion | 5.63 mM | Rebelato et al. 2018, Sci Rep, **PMID 30166558**: islet EC50 = 10.13 mM minus the 4.5 mM threshold |
| `first_phase_gain` | Size of the first-phase burst | 121.3 pM per (mM/min) | Form from Mari et al. 2002, PMID 12388151; **value fitted to the OGTT benchmark** |

The p3 value is the one to remember for interviews: it was corrected from an
uncited number to the one Bergman's paper actually reports.

I didn't re-fetch these sources while writing this guide. The project's citation
audits checked them earlier.

---

## 6. Self-test

1. At fasting (`G = Gb`, `I = Ib`, `X = 0`, no meal), what is `dG`? What is `dI`?
   Why does that matter?
2. Which term in `dG` would still lower glucose in someone who makes no insulin at
   all?
3. Why does insulin act through `X` instead of appearing directly in `dG`?
4. Write `dX` from memory and say in words what each term means.
5. What is insulin sensitivity `SI` in terms of `p2` and `p3`, and what happens to
   it in type 2 diabetes?
6. In `dI`, what is the difference between first-phase and second-phase
   secretion, and which one is zero while glucose is falling?
7. Why does `dI` use `current_frac − basal_frac` instead of `current_frac` alone?
8. If you doubled mTOR activity in the middle of a simulation, what would happen
   to blood glucose? Why?

---

<details>
<summary>Answers (try the questions first)</summary>

1. Both are zero: no meal means `Ra` = 0; `G − Gb` = 0; `X` = 0; and
   `I − Ib` = 0. The body is at a steady state, so the simulation starts flat
   and only moves when a meal or drug arrives. That's what "basal steady state"
   means.
2. The `p1` term (glucose effectiveness). It doesn't depend on insulin.
3. Because insulin's effect lags plasma insulin: it has to reach the tissues and
   trigger signalling first. X is that delayed effect.
4. `dX = −p2·X + p3·max(I − Ib, 0)`. X builds up when insulin is above fasting,
   at a rate set by `p3`, and fades on its own at rate `p2`.
5. `SI = p3 / p2`. It falls in T2D (insulin resistance): the same extra insulin
   produces less glucose-clearing action.
6. First phase is a burst proportional to how fast glucose is rising. Second
   phase is a saturating (Michaelis–Menten) response to how high glucose is.
   First phase is zero while glucose is falling.
7. So that only the *extra*, meal-driven secretion is added. Fasting secretion is
   already balanced by fasting clearance; adding it again would make insulin
   drift upward at rest.
8. Nothing. mTOR is a read-out: no term in `dG`, `dI` or `dX` uses it.

</details>

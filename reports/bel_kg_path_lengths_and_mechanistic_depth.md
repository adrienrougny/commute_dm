# Are the BEL KGs' short paths mechanistically meaningful?

**Question.** The average shortest-path length of the BEL knowledge graphs (AD, PD,
COVID, CBM) is very low (~3.4–4.8 hops). Is this a sign of genuine mechanistic
connectivity — useful for generating hypotheses — or is it a *hairball* artifact?

**Answer (short).** The low average path length is real but misleading. It is a
hairball effect driven by a handful of phenotype/process **super-hubs**. The metric
that actually matters for mechanistic hypotheses — *directed* reachability along
causal edges — is already very low, and it collapses once the hubs are removed. The
global "characteristic path length" oversells these KGs' mechanistic depth.

## Method

Computed on the **projected BEL graph** used throughout `2_05`
(`get_bel_projection`): the counted `abundance/activity` and
`biologicalprocess/pathology` nodes, minus structural-constituent nodes, with the
induced non-structural (non-`HAS_`) relationships. For each KG the path statistics
are measured twice:

- **full projection** — as in the notebook;
- **molecular only** — the same graph with the `biologicalprocess/pathology`
  super-hubs removed (the `path(...)` disease nodes and `bp(...)` process nodes that
  dominate the degree ranking at degree ~500–2000), leaving the molecular backbone
  (protein / complex / abundance / RNA / gene / activity).

Undirected metrics are measured on the largest connected component (LCC); directed
metrics over the ordered node pairs reachable along edge direction.
`dir_reach_frac` is the share of all *n*(*n*−1) ordered pairs that are reachable at
all. All four KGs (AD, PD, COVID, CBM) are computed exactly; CBM is a distinct KG
of 2,425 nodes (it overlaps AD/PD/COVID by only 4.6–8.3 %, i.e. shared common
biology, not containment).

## Results

| KG | undir avg | → no hubs | LCC frac | → no hubs | dir reach frac | → no hubs |
|----|-----------|-----------|----------|-----------|----------------|-----------|
| AD    | 3.88 | **4.30** | 0.79 | **0.59** | 0.174  | **0.065** |
| CBM   | 3.41 | **4.79** | 0.92 | **0.34** | 0.146  | **0.0013** |
| COVID | 4.80 | **7.14** | 0.57 | **0.30** | 0.014  | **0.001** |
| PD    | 4.47 | **5.78** | 0.88 | **0.50** | 0.127  | **0.013** |

(Node/edge counts of the full projection: AD 4099 / 8417, CBM 2084 / 3905,
COVID 2278 / 2879, PD 1930 / 3304. "→ no hubs" removes the
`biologicalprocess/pathology` super-hubs.)

CBM is the most extreme case: it *looks* the most connected of all four (highest
LCC 0.92, shortest average path 3.41), yet **82 % of its edges touch a hub** — drop
the phenotype/process hubs and its molecular backbone falls to just 722 edges among
1351 nodes, the LCC collapses to 0.34, and directed reachability drops ~100× to
0.13 %. Its apparent connectivity is almost entirely hub-mediated.

## Interpretation

1. **The short paths are hub-mediated, and hub paths are not mechanisms.** A 2-hop
   route such as `p(X) → path("Alzheimer Disease") → p(Y)` is topologically short
   but mechanistically vacuous — it only asserts that both entities relate to AD.
   Removing the phenotype/process hubs lengthens the average path *and* fragments
   the giant component (LCC drops to 30–59 %): the molecular backbone underneath is
   genuinely sparse and disconnected. The hairball *is* the aggregator nodes.

2. **Directed reachability is the honest metric, and it is tiny.** A mechanistic
   hypothesis needs a *directed causal chain* A → … → B, not an undirected hop
   count. Even with the hubs in place, only 1.4 % (COVID) to 17 % (AD) of ordered
   pairs are directionally reachable; without the hubs this falls to 0.1–6.5 %. For
   the large majority of node pairs there is simply no causal path — a fact the
   undirected "characteristic path length" hides completely.

3. **The graphs are not too dense — they are too shallow and too shortcut-y.** These
   are curated-literature KGs: broad but shallow, with correlative (`ASSOCIATION`,
   `*_CORRELATION`) and hierarchical (`IS_A`) edges plus phenotype hubs adding
   non-mechanistic shortcuts on top of a thin causal skeleton.

## Consequences for mechanistic hypothesis generation

Global path statistics on these KGs should **not** be read as "mechanistic
capacity". For hypothesis generation the useful view is the opposite of the
hairball — a directed, causal-only, hub-excluded graph in which a path actually
encodes a candidate mechanism:

- keep only directed causal edges (`INCREASES` / `DECREASES` / `DIRECTLY_*` /
  `REGULATES`), drop `ASSOCIATION` / correlations / `IS_A`;
- treat phenotype / disease / process nodes as **terminals, not waypoints**
  (blacklist them from the interior of paths);
- restrict intermediates to molecular species (protein / complex / abundance /
  RNA / gene).

The expected outcome is a sparser, more fragmented, but mechanistically honest
graph where most pairs are unreachable — that is a property of the KGs' limited
causal coverage, not a processing bug.

The COMMUTE interface pipeline already works this way: `queries.get_subgraph` uses
directed `downstream` / `upstream` modes, a `relationship_types` filter (the
`INFLUENCES` causal subset), and blacklists, and the interface analysis builds
directed upstream-COVID / downstream-PD subgraphs rather than relying on any global
average. The `2_05` path-length table is descriptive only.

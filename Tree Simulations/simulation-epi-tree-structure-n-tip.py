# -*- coding: utf-8 -*-
"""
full_tree_sim.py
================
Simulator + sampled-subtree extractor whose edge table wires directly
into the Theorem 4.1 / Eq. 15 likelihood, extended to support random
contact degree K.

Degree distributions
--------------------
degree_dist = 'fixed'      : every individual draws K = k (original behaviour)
degree_dist = 'poisson'    : K ~ Poisson(k),    min 1
degree_dist = 'negbinom'   : K ~ NegBin(mean=k, dispersion=k_dispersion),
                             where k_dispersion is the r/size parameter.
                             Variance = k + k^2/k_dispersion.
                             k_dispersion -> inf  recovers Poisson(k).
degree_dist = 'geometric'  : K ~ Geometric with mean k (supported on {1,2,...})
                             P(K=j) = (1/k)(1-1/k)^{j-1}

Design principle
----------------
The root FOCAL individual's degree k is drawn once at infection time and
stored as ind.k.  Every newborn independently draws its own degree
from the same distribution.  For 'fixed', ind.k = Para.k always, so
the simulation is identical to the original code.

Edge table schema (unchanged from original)
-------------------------------------------
rep_id, edge_type, lineage_type, par_id, chi_id,
tau_a, tau_b, delta, s, c, j_obs, k, branch_rate,
ind_id, infector_id, seg_idx, n_segs,
n_tips, n_internal, root_state, k_root, t_root_abs

The 'k' column now stores the per-individual degree (previously always
Para.k; now varies across individuals for random distributions).
All other columns and their semantics are unchanged.
"""

# ── standard library ─────────────────────────────────────────────────────────
from random import seed, random, randrange, expovariate
import bisect
import math

# ── third-party ──────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
import tqdm


# ══════════════════════════════════════════════════════════════════════════════
#  DegreeDistribution  —  degree sampler
# ══════════════════════════════════════════════════════════════════════════════

class DegreeDistribution:
    """
    Encapsulates sampling of the contact degree K.

    Parameters
    ----------
    dist : str
        One of 'fixed', 'poisson', 'negbinom', 'geometric'.
    mean_k : int or float
        Mean of the distribution (= k for 'fixed').
    dispersion : float, optional
        Size/r parameter for 'negbinom' (variance = mean + mean^2/dispersion).
        Ignored for all other distributions.
    k_min : int, optional
        Minimum degree enforced by rejection (default 1).
        For 'fixed', k_min has no effect (k is returned exactly).
    """

    SUPPORTED = ('fixed', 'poisson', 'negbinom', 'geometric')

    def __init__(self, dist='fixed', mean_k=4, dispersion=1.0, k_min=1):
        if dist not in self.SUPPORTED:
            raise ValueError(
                f"degree_dist must be one of {self.SUPPORTED}; got '{dist}'")
        self.dist       = dist
        self.mean_k     = float(mean_k)
        self.dispersion = float(dispersion)   # used only for negbinom
        self.k_min      = int(k_min)
        if self.mean_k <= 0:
            raise ValueError("mean_k must be positive")

    def draw(self):
        """Return a single integer degree drawn from the distribution."""
        if self.dist == 'fixed':
            return int(self.mean_k)

        # For stochastic distributions, draw and enforce k_min
        while True:
            if self.dist == 'poisson':
                k = int(np.random.poisson(lam=self.mean_k))
            elif self.dist == 'negbinom':
                # np.random.negative_binomial(n, p) has mean n*(1-p)/p
                # We want mean = mean_k, dispersion = n (size).
                # p = n / (n + mean_k)
                n = self.dispersion
                p = n / (n + self.mean_k)
                k = int(np.random.negative_binomial(n=n, p=p))
            elif self.dist == 'geometric':
                # Geometric on {1,2,...} with mean mean_k:
                # P(K=j) = (1/mean_k)(1 - 1/mean_k)^{j-1}
                # Equivalent to 1 + Geometric-on-{0,1,...}(p=1/mean_k)
                p = 1.0 / self.mean_k
                k = int(np.random.geometric(p=p))   # numpy gives {1,2,...}
            else:
                raise ValueError(f"Unknown distribution: {self.dist}")

            if k >= self.k_min:
                return k

    def description(self):
        if self.dist == 'fixed':
            return f"fixed(k={int(self.mean_k)})"
        elif self.dist == 'poisson':
            return f"Poisson(lambda={self.mean_k})"
        elif self.dist == 'negbinom':
            return f"NegBin(mean={self.mean_k}, dispersion={self.dispersion})"
        elif self.dist == 'geometric':
            return f"Geometric(mean={self.mean_k})"


# ══════════════════════════════════════════════════════════════════════════════
#  Para
# ══════════════════════════════════════════════════════════════════════════════

class Para:
    """
    Epidemic parameters.

    Parameters
    ----------
    k : int or float
        For degree_dist='fixed': the exact contact degree.
        For all others: the mean of the degree distribution.
    beta : float
        Per-contact transmission rate.
    mu : float
        Undetected removal rate.
    sigma : float
        Detected (observed) removal rate.
    burn_time : float
        Time at which the focal individual is tagged.
    time_horizont : float
        Maximum simulation time.
    max_clade_size : int
        Stop tagging clade if it exceeds this size.
    degree_dist : str
        One of 'fixed', 'poisson', 'negbinom', 'geometric'.
    k_dispersion : float
        Dispersion parameter for 'negbinom' (ignored otherwise).
    k_min : int
        Minimum degree enforced for random distributions.
    """

    def __init__(self,
                 k=4,
                 beta=1.5,
                 mu=0.5,
                 sigma=0.5,
                 burn_time=1.5,
                 time_horizont=10.0,
                 max_clade_size=2000,
                 degree_dist='fixed',
                 k_dispersion=1.0,
                 k_min=1):

        self.k              = k
        self.beta           = beta
        self.mu             = mu
        self.sigma          = sigma
        self.burn_time      = burn_time
        self.time_horizont  = time_horizont
        self.max_clade_size = max_clade_size
        self.pobs           = sigma / (sigma + mu)

        self.degree_dist_obj = DegreeDistribution(
            dist=degree_dist, mean_k=k,
            dispersion=k_dispersion, k_min=k_min)

        # filled after tagging
        self.tagState = (-1, -1)   # (k_of_focal, j_at_burntime)
        self.tagTime  = burn_time

    def draw_k(self):
        """Draw a contact degree for a newly infected individual."""
        return self.degree_dist_obj.draw()

    @property
    def degree_dist(self):
        return self.degree_dist_obj.dist


# ══════════════════════════════════════════════════════════════════════════════
#  Event
# ══════════════════════════════════════════════════════════════════════════════

class Event:
    def __init__(self, myType, idA, idB, eventTime):
        self.type       = myType
        self.idA        = idA
        self.idB        = idB
        self.eventTime  = eventTime
        self.deactivate = False

    def __lt__(self, other):
        return self.eventTime < other.eventTime


# ══════════════════════════════════════════════════════════════════════════════
#  Individual
# ══════════════════════════════════════════════════════════════════════════════

class Individual:
    """
    Attributes used by SampledSubtree
    -----------------------------------
    my_id               int
    infectorID          int   (-1 for epidemic root)
    infectionEventTime  float (time this individual was infected)
    recoveryTime        float or None
    myRecoverType       0 = unobserved, 1 = observed, "x" = still infected
    contactees          list[int]   IDs infected BY this individual (all-time)
    contactTimes        list[float] times of those infections (same order)
    k                   int   contact degree of THIS individual (per-individual)
    phyState            (k, j) at last update  — j = #downstreams infected so far
    IamTagged           "yes" / "no"
    """

    def __init__(self, my_id, infectorID, IamTagged, pop_obj):
        self.my_id              = my_id
        self.infectorID         = infectorID
        self.infectionEventTime = -1.0
        self.pop_obj            = pop_obj
        self.paraObj            = pop_obj.my_para
        self.state              = "S"
        self.contactees         = []
        self.contactTimes       = []
        self.IamTagged          = IamTagged
        self.myRecoveryEvent    = None
        self.myDownstreams      = []

        self.j          = 0
        # k will be set when infectMe is called; default to para.k
        self.k          = self.paraObj.k
        self.phyState   = (self.k, self.j)

        self.myRecoverType = "x"
        self.recoveryTime  = None

    # ── infection ─────────────────────────────────────────────────────────────

    def infectMe(self, infectorID, Tag, globalTime):
        # Draw THIS individual's contact degree from the distribution.
        # For 'fixed' this always returns Para.k (identical to original).
        self.k = self.paraObj.draw_k()

        for _ in range(self.k):
            t_contact = globalTime + expovariate(self.paraObj.beta)
            ind_id    = self.pop_obj.generateIndividual(self.my_id, Tag, globalTime)
            self.myDownstreams.append(ind_id)
            ev = Event("contact", self.my_id, ind_id, t_contact)
            self.pop_obj.registerEvent(ev)

        self.state              = "I"
        self.IamTagged          = Tag
        self.infectionEventTime = globalTime
        self.pop_obj.increaseInfecteds()
        self.phyState = (self.k, 0)

        rate    = self.paraObj.mu + self.paraObj.sigma
        t_rec   = globalTime + expovariate(rate)
        ev_type = "ObsRcvr" if random() < self.paraObj.pobs else "SpntRcvr"
        ev = Event(ev_type, self.my_id, -1, t_rec)
        self.pop_obj.registerEvent(ev)
        self.myRecoveryEvent = ev

    # ── contact ───────────────────────────────────────────────────────────────

    def registerContact(self, ind_id, contactTime, globalTime):
        self.contactees.append(ind_id)
        self.contactTimes.append(contactTime)
        self.j       += 1
        self.phyState = (self.k, self.j)

    # ── recovery ──────────────────────────────────────────────────────────────

    def recoverUnobserved(self, globalTime):
        self.myRecoverType = 0
        self.state         = "R"
        self.recoveryTime  = globalTime
        self.pop_obj.decreaseInfecteds()
        self.phyState = (self.k, self.j)

    def recoverObserved(self, globalTime):
        self.myRecoverType = 1
        self.state         = "R"
        self.recoveryTime  = globalTime
        self.pop_obj.decreaseInfecteds()
        self.phyState = (self.k, self.j)


# ══════════════════════════════════════════════════════════════════════════════
#  Population
# ══════════════════════════════════════════════════════════════════════════════

class Population:

    def __init__(self, para_obj):
        self.my_para             = para_obj
        self.time                = 0.0
        self.eventQueue          = []
        self.pop_size            = 0
        self.indiv_list          = []
        self.noInfecteds         = 0
        self.observedIndi        = []
        self.taggedRoot          = None
        self.didDoTagIndi        = 0
        self.didNotFindIndiToTag = 0
        self.tagged_clade_size   = 0

        root = Individual(self.pop_size, -1, "no", self)
        self.indiv_list.append(root)
        self.pop_size += 1
        root.infectMe(-1, "no", 0.0)

        ev = Event("tagAnIndividual", -1, -1, para_obj.burn_time)
        self.registerEvent(ev)

    def generateIndividual(self, infectorID, Tag, globalTime):
        nid = self.pop_size
        ind = Individual(nid, infectorID, Tag, self)
        ind.infectionEventTime = globalTime
        self.indiv_list.append(ind)
        self.pop_size += 1
        if Tag == 'yes':
            self.tagged_clade_size += 1
        return nid

    def increaseInfecteds(self): self.noInfecteds += 1
    def decreaseInfecteds(self): self.noInfecteds -= 1

    def registerEvent(self, ev):
        bisect.insort(self.eventQueue, ev, key=lambda x: x.eventTime)

    def handleEvent(self):
        if not self.eventQueue:
            return
        ev        = self.eventQueue.pop(0)
        self.time = ev.eventTime
        if ev.deactivate:
            return

        if self.didDoTagIndi > 0 and ev.idA >= 0:
            if self.indiv_list[ev.idA].IamTagged != "yes":
                return

        if ev.type == "contact":
            indA = self.indiv_list[ev.idA]
            indB = self.indiv_list[ev.idB]
            if indA.state == "I" and indB.state == "S":
                tag = indA.IamTagged
                indB.infectMe(ev.idA, tag, self.time)
                indA.registerContact(ev.idB, self.time, self.time)

        elif ev.type == "SpntRcvr":
            self.indiv_list[ev.idA].recoverUnobserved(self.time)

        elif ev.type == "ObsRcvr":
            ind = self.indiv_list[ev.idA]
            ind.recoverObserved(self.time)
            if ind.IamTagged == "yes":
                self.observedIndi.append({
                    "id"       : ind.my_id,
                    "phyState" : ind.phyState,
                    "t_rec"    : self.time,
                })

        elif ev.type == "tagAnIndividual":
            active = [ind for ind in self.indiv_list if ind.state == "I"]
            if not active:
                self.didNotFindIndiToTag = 1
                return
            tagged              = active[randrange(len(active))]
            tagged.IamTagged    = "yes"
            self.taggedRoot     = tagged
            self.didDoTagIndi   = 1
            self.tagged_clade_size = 1
            self.my_para.tagState  = (tagged.k, tagged.j)   # (k, j) at burn_time

    def taggedCladeDone(self):
        for ind in self.indiv_list:
            if ind.IamTagged == "yes" and ind.state == "I":
                return False
        return True

    def cladeSizeLimitReached(self):
        return self.tagged_clade_size > self.my_para.max_clade_size


# ══════════════════════════════════════════════════════════════════════════════
#  SampledSubtree  —  produces the edge table for Eq. 15
# ══════════════════════════════════════════════════════════════════════════════

class SampledSubtree:
    """
    Builds the edge table whose rows feed directly into Eq. 15.

    Changes from the original version:
    - branch_rate uses ind.k (per-individual degree) instead of a fixed
      global self.k.  This is the only structural change; all other logic
      is identical.
    - The 'k' column in the edge table stores ind.k for each individual.
    """

    def __init__(self, population: Population):
        self.pop    = population
        self.indivs = population.indiv_list
        self.root   = population.taggedRoot

    def _ind(self, nid):
        return self.indivs[nid]

    def _is_observed(self, nid):
        return self._ind(nid).myRecoverType == 1

    def _entry_time(self, nid):
        if nid == self.root.my_id:
            return self.pop.my_para.burn_time
        return self._ind(nid).infectionEventTime

    def _entry_state(self, nid):
        if nid == self.root.my_id:
            return self.pop.my_para.tagState[1]   # j at burn_time
        return 0

    def _build_subtree(self):
        obs_ids = [
            ind.my_id for ind in self.indivs
            if ind.IamTagged == "yes" and ind.myRecoverType == 1
        ]
        if not obs_ids:
            return set(), obs_ids

        sub = set()
        for oid in obs_ids:
            cur = oid
            while cur != -1 and cur is not None:
                if cur in sub:
                    break
                ind = self._ind(cur)
                if ind.infectionEventTime < 0:
                    break
                sub.add(cur)
                if cur == self.root.my_id:
                    break
                cur = ind.infectorID

        return sub, obs_ids

    def _clade_infections(self, nid, sub):
        ind        = self._ind(nid)
        entry_time = self._entry_time(nid)

        events = []
        for child_id, t in zip(ind.contactees, ind.contactTimes):
            if child_id in sub and t >= entry_time - 1e-12:
                events.append((t, child_id))

        events.sort(key=lambda x: x[0])
        return events

    def build(self):
        sub, obs_ids = self._build_subtree()
        if not obs_ids:
            return self

        edges    = []
        root_id  = self.root.my_id

        for nid in sub:
            ind          = self._ind(nid)
            entry_t      = self._entry_time(nid)
            entry_s      = self._entry_state(nid)
            clade_events = self._clade_infections(nid, sub)
            j_final      = ind.phyState[1]
            t_rec        = ind.recoveryTime
            k_ind        = ind.k        # per-individual degree (KEY CHANGE)
            beta         = self.pop.my_para.beta

            if t_rec is None:
                continue

            if nid == root_id:
                etype = "root"
            elif clade_events:
                etype = "internal"
            else:
                etype = "tip"

            seg_times  = [entry_t] + [t for t, _ in clade_events] + [t_rec]
            seg_states = list(range(entry_s, entry_s + len(clade_events) + 2))
            n_seg      = len(clade_events) + 1

            for m in range(n_seg):
                tau_a = seg_times[m]
                tau_b = seg_times[m + 1]
                s     = seg_states[m]
                c     = j_final if m == n_seg - 1 else seg_states[m + 1]
                delta = max(tau_b - tau_a, 0.0)

                # branch_rate: uses k_ind (per-individual) instead of global k
                if m == 0:
                    br = float("nan")
                elif m < n_seg - 1:
                    br = (k_ind - s) * beta   # KEY CHANGE: k_ind not self.k
                else:
                    br = float("nan")

                is_tip_seg   = (m == n_seg - 1) and self._is_observed(nid)
                j_obs        = float(c) if is_tip_seg else float("nan")

                lineage_type = "newborn" if (m == 0 and nid != root_id) else "continuing"
                seg_etype    = etype

                edges.append({
                    "edge_type"    : seg_etype,
                    "lineage_type" : lineage_type,
                    "par_id"       : nid,
                    "chi_id"       : nid,
                    "tau_a"        : tau_a,
                    "tau_b"        : tau_b,
                    "delta"        : delta,
                    "s"            : int(s),
                    "c"            : int(c),
                    "j_obs"        : j_obs,
                    "k"            : k_ind,     # per-individual degree
                    "branch_rate"  : br,
                    "ind_id"       : nid,
                    "infector_id"  : ind.infectorID,
                    "seg_idx"      : m,
                    "n_segs"       : n_seg,
                })

        self.edges = edges
        return self

    def to_dataframe(self):
        if not self.edges:
            return None
        df = pd.DataFrame(self.edges)
        df = df.sort_values(["tau_a", "ind_id", "seg_idx"]).reset_index(drop=True)
        df["n_tips"]     = int((df["edge_type"] == "tip").sum())
        df["n_internal"] = int(df["edge_type"].isin(["internal", "root"]).sum())
        return df


# ══════════════════════════════════════════════════════════════════════════════
#  Likelihood helper  —  unchanged from original
# ══════════════════════════════════════════════════════════════════════════════

def log_likelihood_from_df(df_rep, D_func, E_func, pi_vec, k):
    """
    Compute log-likelihood for a single replicate from its edge table.
    Implements Theorem 4.1 / Eq. 15 (fixed-k version for validation).

    For the random-K likelihood, use D_func and E_func that solve the
    coupled random-K ODE system (eqs E_random, D_random from the LaTeX),
    and replace pi_vec with the mixed equilibrium pi_hat from pi_i_k_random.py.
    """
    ll = 0.0

    for _, row in df_rep.iterrows():
        d_val = D_func(row["s"], row["c"], row["delta"])
        if d_val <= 0:
            return -math.inf
        ll += math.log(d_val)

        br = row["branch_rate"]
        if not math.isnan(br):
            if br <= 0:
                return -math.inf
            ll += math.log(br)

    root_row = df_rep[
        (df_rep["edge_type"] == "root") & (df_rep["seg_idx"] == 0)
    ]
    if root_row.empty:
        return -math.inf

    root_row = root_row.iloc[0]
    t_total  = df_rep["tau_b"].max() - root_row["tau_a"]

    cond_num = sum(
        pi_vec[i] * D_func(i, root_row["c"], root_row["delta"])
        for i in range(k + 1)
    )
    cond_den = sum(
        pi_vec[i] * (1.0 - E_func(i, t_total))
        for i in range(k + 1)
    )
    if cond_num <= 0 or cond_den <= 0:
        return -math.inf

    ll += math.log(cond_num) - math.log(cond_den)
    return ll


# ══════════════════════════════════════════════════════════════════════════════
#  Verification
# ══════════════════════════════════════════════════════════════════════════════

def verify_edge_table(df):
    """
    Run consistency checks on the edge table.
    For random-K simulations, check 6 (branch_rate == (k-s)*beta) is
    not enforced quantitatively since k varies per individual, but the
    sign and finiteness are still checked.
    """
    errors = []

    tips = df[df["edge_type"] == "tip"]
    bad_jobs = tips[tips["j_obs"].isna()]
    if not bad_jobs.empty:
        errors.append(f"Check 1 FAIL: {len(bad_jobs)} tip rows have NaN j_obs")

    bad_delta = tips[tips["delta"] <= 0]
    if not bad_delta.empty:
        errors.append(f"Check 2 FAIL: {len(bad_delta)} tip rows have delta <= 0")

    for (rep, ind), grp in df.groupby(["rep_id", "ind_id"]):
        grp = grp.sort_values("seg_idx")
        for i in range(len(grp) - 1):
            rowA = grp.iloc[i]; rowB = grp.iloc[i + 1]
            if abs(rowA["tau_b"] - rowB["tau_a"]) > 1e-12:
                errors.append(
                    f"Check 3 FAIL: rep={rep} ind={ind} "
                    f"seg {int(rowA['seg_idx'])} tau_b != next seg tau_a")
                break

    newborns = df[df["lineage_type"] == "newborn"]
    bad_s = newborns[newborns["s"] != 0]
    if not bad_s.empty:
        errors.append(f"Check 4 FAIL: {len(bad_s)} newborn rows have s != 0")

    last_segs   = df[df["seg_idx"] == df["n_segs"] - 1]
    bad_br_last = last_segs[last_segs["branch_rate"].notna()]
    if not bad_br_last.empty:
        errors.append(
            f"Check 5 FAIL: {len(bad_br_last)} last-segment rows have non-NaN branch_rate")
    bad_br_nb = newborns[newborns["branch_rate"].notna()]
    if not bad_br_nb.empty:
        errors.append(
            f"Check 5 FAIL: {len(bad_br_nb)} newborn rows have non-NaN branch_rate")

    mid_cont = df[
        (df["lineage_type"] == "continuing") &
        (df["seg_idx"] > 0) &
        (df["seg_idx"] < df["n_segs"] - 1)
    ]
    bad_c = mid_cont[mid_cont["c"] != mid_cont["s"] + 1]
    if not bad_c.empty:
        errors.append(
            f"Check 7 FAIL: {len(bad_c)} intermediate continuing rows have c != s+1")

    # Check 8 (new): k must be >= j_obs for all tip rows
    tip_last = df[~df["j_obs"].isna()]
    bad_k = tip_last[tip_last["k"] < tip_last["j_obs"]]
    if not bad_k.empty:
        errors.append(
            f"Check 8 FAIL: {len(bad_k)} tip rows have k < j_obs "
            f"(individual infected more contacts than their degree)")

    if errors:
        print("=== Edge table verification FAILED ===")
        for e in errors: print(" ", e)
    else:
        print("Edge table verification PASSED (all checks)")

    return errors


# ══════════════════════════════════════════════════════════════════════════════
#  Runner
# ══════════════════════════════════════════════════════════════════════════════

def run_full_tree(para_factory, N: int = 1000) -> pd.DataFrame:
    """
    Run N replicates and return the combined edge table.
    para_factory must return a Para object; all degree distributions
    are supported via Para.degree_dist_obj.
    """
    all_edges    = []
    rep_id       = 0
    n_skip_notag = 0
    n_skip_noobs = 0

    for _ in tqdm.tqdm(range(N)):
        para = para_factory()
        popu = Population(para)

        done = False
        while not done:
            if not popu.eventQueue:
                break
            popu.handleEvent()
            if popu.didDoTagIndi > 0 and popu.taggedCladeDone():
                done = True
            if popu.didDoTagIndi > 0 and popu.cladeSizeLimitReached():
                done = True
            if popu.time > para.time_horizont:
                done = True
            if popu.didNotFindIndiToTag > 0:
                done = True

        if popu.didDoTagIndi == 0 or popu.taggedRoot is None:
            n_skip_notag += 1
            continue
        if not popu.observedIndi:
            n_skip_noobs += 1
            continue

        st = SampledSubtree(popu).build()
        if not hasattr(st, 'edges') or not st.edges:
            n_skip_noobs += 1
            continue

        df = st.to_dataframe()
        if df is None or df.empty:
            n_skip_noobs += 1
            continue

        df.insert(0, "rep_id", rep_id)
        df["root_state"] = para.tagState[1]
        df["k_root"]     = para.tagState[0]
        df["t_root_abs"] = popu.taggedRoot.infectionEventTime

        all_edges.append(df)
        rep_id += 1

    print(f"\nValid replicates  : {rep_id}")
    print(f"Skipped (no tag)  : {n_skip_notag}")
    print(f"Skipped (no obs)  : {n_skip_noobs}")

    return pd.concat(all_edges, ignore_index=True) if all_edges else pd.DataFrame()


# ══════════════════════════════════════════════════════════════════════════════
#  Entry point
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import sys

    # ── choose degree distribution from command line or edit here ────────────
    # Examples:
    #   python full_tree_sim.py fixed      -> fixed k=4 (simulation study)
    #   python full_tree_sim.py poisson    -> Poisson(lambda=4)
    #   python full_tree_sim.py negbinom   -> NegBin(mean=4, dispersion=0.5)
    #   python full_tree_sim.py geometric  -> Geometric(mean=4)

    dist = sys.argv[1] if len(sys.argv) > 1 else "fixed"

    # Shared parameters
    BETA      = 1.5
    MU        = 0.5
    SIGMA     = 0.5
    K_MEAN    = 4       # mean degree (= exact k for 'fixed')
    DISP      = 0.5     # dispersion for negbinom (ignored for others)
    BURN      = 0.6
    HORIZON   = 20.0
    MAX_CLADE = 250
    N_REPS    = 1000
    SET_SEED  = 1234

    def make_para():
        return Para(
            k=K_MEAN, beta=BETA, mu=MU, sigma=SIGMA,
            burn_time=BURN, time_horizont=HORIZON, max_clade_size=MAX_CLADE,
            degree_dist=dist, k_dispersion=DISP, k_min=1)

    seed(SET_SEED)
    np.random.seed(SET_SEED)

    p = make_para()
    print(f"Degree distribution: {p.degree_dist_obj.description()}")
    print(f"Parameters: beta={BETA}, mu={MU}, sigma={SIGMA}, "
          f"R0(k_mean)={K_MEAN*BETA/(MU+SIGMA):.2f}")
    print()

    df = run_full_tree(make_para, N=N_REPS)

    if df.empty:
        print("No valid replicates.")
    else:
        out = f"full_tree_edges_{dist}.csv"
        df.to_csv(out, index=False)
        print(f"\nSaved {len(df)} edge rows from "
              f"{df['rep_id'].nunique()} replicates -> {out}")

        print("\n=== Degree distribution in edge table ===")
        print(df.groupby("k").size().rename("n_rows").to_string())

        print("\n=== Edge type / lineage_type breakdown ===")
        print(df.groupby(["edge_type", "lineage_type"]).size().to_string())

        print("\n=== j_obs distribution (tip states) ===")
        tip_last = df[~df["j_obs"].isna()]
        print(tip_last["j_obs"].value_counts().sort_index().to_string())

        print("\n=== Verification ===")
        verify_edge_table(df)

        print("\n=== First replicate edge table ===")
        first = df[df["rep_id"] == 0].sort_values(["ind_id", "seg_idx"])
        cols  = ["edge_type", "lineage_type", "ind_id", "seg_idx",
                 "tau_a", "tau_b", "delta", "s", "c", "j_obs", "k", "branch_rate"]
        print(first[cols].to_string(index=False))
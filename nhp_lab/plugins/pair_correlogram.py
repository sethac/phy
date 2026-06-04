"""phy PairCorrelogramView plugin: jitter-corrected pair triptych - v0.1

A dedicated correlogram view for INSPECTING A PAIR of clusters. When EXACTLY
two clusters are selected it draws a 1x3 triptych:

    [ ACG_A ] [ CCG(raw) + CCG(jitter-corrected) ] [ ACG_B ]

Why a separate view (not phy's built-in CorrelogramView)?
  1. FULL SPIKE TRAIN. phy's CorrelogramView caps spikes via
     self.selector(self.n_spikes_correlograms, ...) (default 100k), so the
     CCG/ACG shape for a pair can be sampling noise. This view bypasses that
     cap and pulls EVERY spike for the two selected clusters straight from the
     model, so the shape is the true shape.
  2. JITTER OVERLAY. The CCG panel overlays the raw CCG with the Harrison-Geman
     ANALYTIC jitter-corrected CCG (Moore lab / Panichello method). The jitter
     correction here is the lag-domain form documented in memory
     reference_ccg_method.md and used offline in scripts/15 and scripts/16:
     subtract the raw CCG convolved with the triangular kernel of the 25 ms
     jitter window (the continuous / flat-PSTH degenerate case). Geometric-mean
     firing-rate normalization, matching the offline detector.

Defaults match the offline detector: +/-25 ms window, 0.5 ms bins, 2 ms
refractory marker. All of phy's live controls (set_bin 'cb', set_window 'cw',
set_refractory_period 'cr', ctrl/alt+wheel) are INHERITED unchanged, so the
Lord can re-bin/re-window any pair in the GUI.

This is a DISPLAY-ONLY view. It never touches the .ap.bin or the sort.

Install: this file lives in C:/Users/sa4257/.phy/plugins/ (already on phy's
Plugins.dirs). Activate by adding 'PairCCGPlugin' to the plugins list for the
TemplateGUI in ~/.phy/phy_config.py. The plugin registers a view creator under
the name 'PairCorrelogramView'; add that view from the GUI's View menu.
"""
import logging

import numpy as np

from phy import IPlugin
from phylib.utils import Bunch
from phy.utils.color import selected_cluster_color, _override_hsv, add_alpha
from phy.cluster.views.correlogram import CorrelogramView

logger = logging.getLogger(__name__)

# Offline-detector defaults (see scripts/15, scripts/16, reference_ccg_method.md)
JITTER_MS = 25.0           # Harrison-Geman jitter window (Moore / Panichello)
DEFAULT_BIN_S = 0.5e-3     # 0.5 ms bins
DEFAULT_WINDOW_S = 50e-3   # +/-25 ms -> 50 ms full window
DEFAULT_REFRACT_S = 2e-3   # 2 ms refractory marker


# -----------------------------------------------------------------------------
# Correlogram + jitter math (mirrors scripts/15_jitter_ccg_pairs_ks4catgt.py)
# -----------------------------------------------------------------------------

def _ccg_from_times(ta, tb, win_s, bin_s, exclude_self=False):
    """CCG counts from spike times (vectorized, single histogram).

    ta, tb : sorted spike times in seconds.
    Returns counts over lags -nlags..+nlags (length 2*nlags+1).
    exclude_self drops the zero-distance self-pairs (use for ACGs).
    """
    nlags = int(round(win_s / bin_s))
    edges = (np.arange(-nlags, nlags + 2) - 0.5) * bin_s
    if len(ta) == 0 or len(tb) == 0:
        return np.zeros(2 * nlags + 1)
    lo = np.searchsorted(tb, ta - win_s, side="left")
    hi = np.searchsorted(tb, ta + win_s, side="right")
    nc = hi - lo
    if nc.sum() == 0:
        return np.zeros(2 * nlags + 1)
    keep = nc > 0
    ta_k = ta[keep]; lo_k = lo[keep]; n_k = nc[keep]
    tb_idx = np.repeat(lo_k, n_k) + (
        np.arange(n_k.sum()) - np.repeat(np.cumsum(n_k) - n_k, n_k))
    diffs = tb[tb_idx] - np.repeat(ta_k, n_k)
    if exclude_self:
        diffs = diffs[diffs != 0]
    return np.histogram(diffs, bins=edges)[0].astype(float)


def _jitter_kernel(jit_s, bin_s):
    """Triangular kernel = autocorrelation of the uniform jitter box.

    Jittering each of two trains within a jit_s window smooths the expected CCG
    by the convolution of two boxes = a triangle of half-width jit_s. Subtracting
    this jitter-smoothed CCG from the raw CCG is the Harrison-Geman correction in
    the lag domain (valid for the continuous / flat-PSTH case). Normalized to sum
    to 1 so it is a pure smoother (preserves total mass).
    """
    bw = max(1, int(round(jit_s / bin_s)))
    box = np.ones(bw) / bw
    tri = np.convolve(box, box)        # triangle, length 2*bw-1
    return tri / tri.sum()


def _jitter_correct(raw, jit_s, bin_s):
    """Return (jittered_ccg, corrected_ccg) for a raw CCG.

    jittered = raw convolved with the triangular jitter kernel (the expected CCG
    under H&G jitter); corrected = raw - jittered. Removes structure slower than
    the jitter window, keeps fine sub-window peaks.
    """
    kern = _jitter_kernel(jit_s, bin_s)
    jittered = np.convolve(raw, kern, mode="same")
    return jittered, raw - jittered


# -----------------------------------------------------------------------------
# The view
# -----------------------------------------------------------------------------

class PairCorrelogramView(CorrelogramView):
    """ACG_A | CCG(raw+jitter) | ACG_B for exactly two selected clusters.

    Subclasses CorrelogramView so it inherits the bin/window/refractory live
    controls, the histogram/line/text visuals, and the grid canvas. We override
    data assembly (full spike trains, geometric-mean norm) and plotting (1x3
    layout + jitter overlay + hint when not exactly 2 clusters).
    """

    # We only ever look at a pair.
    max_n_clusters = 2

    _default_position = 'left'

    # Tight defaults matching the offline detector.
    bin_size = DEFAULT_BIN_S
    window_size = DEFAULT_WINDOW_S
    refractory_period = DEFAULT_REFRACT_S

    def __init__(self, model=None, spike_clusters_func=None, **kwargs):
        # model: phy Model (spike_times, sample_rate, duration).
        # spike_clusters_func: callable -> current spike_clusters array (so we
        #   respect live merges/splits in the GUI). Falls back to model's.
        self.model = model
        self._spike_clusters_func = spike_clusters_func
        # CorrelogramView needs correlograms/firing_rate callables in its ctor,
        # but we override get_clusters_data so they are never called. Pass
        # harmless stubs to satisfy the parent constructor.
        super(PairCorrelogramView, self).__init__(
            correlograms=self._stub_correlograms,
            firing_rate=None,
            sample_rate=model.sample_rate,
            **kwargs)

    # -- stubs / helpers ------------------------------------------------------

    def _stub_correlograms(self, cluster_ids, bin_size, window_size):
        # Never actually used (get_clusters_data is overridden), but keep a
        # shape-correct return in case some parent path calls it.
        n = len(cluster_ids)
        nbins = int(round(window_size / bin_size))
        return np.zeros((n, n, 2 * nbins + 1))

    def _spike_clusters(self):
        if self._spike_clusters_func is not None:
            return self._spike_clusters_func()
        return self.model.spike_clusters

    def _full_spike_times(self, cluster_id):
        """FULL spike train (seconds) for one cluster - NO subsampling cap."""
        sc = self._spike_clusters()
        st = self.model.spike_times[sc == cluster_id]
        return np.sort(np.asarray(st, dtype=np.float64))

    # -- data assembly --------------------------------------------------------

    def get_clusters_data(self, load_all=None):
        """Build the triptych bunches: ACG_A, CCG(raw+jitter), ACG_B.

        Returns 3 Bunches laid out in a 1x3 grid (box_index = (0, col)). Each
        carries .correlogram (raw, for the parent histogram path) and, for the
        CCG, .jitter_ccg for the overlay. Returns [] if not exactly 2 clusters.
        """
        if len(self.cluster_ids) != 2:
            return []

        a, b = self.cluster_ids
        ta = self._full_spike_times(a)
        tb = self._full_spike_times(b)
        win_s = self.window_size * 0.5      # phy window_size is the FULL width
        bin_s = self.bin_size
        dur = float(getattr(self.model, 'duration', 0.0)) or (
            max(ta[-1] if len(ta) else 0.0, tb[-1] if len(tb) else 0.0) + 1.0)

        # Raw correlograms (counts).
        acg_a = _ccg_from_times(ta, ta, win_s, bin_s, exclude_self=True)
        acg_b = _ccg_from_times(tb, tb, win_s, bin_s, exclude_self=True)
        ccg = _ccg_from_times(ta, tb, win_s, bin_s)

        # Geometric-mean firing-rate normalization (Bair/Zohary; Moore norm),
        # matching scripts/15. Turns counts into a rate-normalized CCG so the
        # raw and jitter curves are on a comparable, rate-controlled scale.
        fr_a = len(ta) / dur if dur > 0 else 0.0
        fr_b = len(tb) / dur if dur > 0 else 0.0
        gm = np.sqrt(fr_a * fr_b)
        norm = gm * dur * bin_s
        if norm > 0:
            ccg_n = ccg / norm
        else:
            ccg_n = ccg.astype(float)

        # Harrison-Geman jitter correction on the normalized CCG.
        _, ccg_corr = _jitter_correct(ccg_n, JITTER_MS * 1e-3, bin_s)

        nbins = len(ccg)
        bunchs = []

        # Column 0: ACG of A.
        bunchs.append(self._acg_bunch(acg_a, col=0, cluster_index=0))
        # Column 1: CCG raw (normalized) + jitter-corrected overlay.
        bunchs.append(self._ccg_bunch(ccg_n, ccg_corr, col=1))
        # Column 2: ACG of B.
        bunchs.append(self._acg_bunch(acg_b, col=2, cluster_index=1))
        return bunchs

    def _acg_bunch(self, acg, col, cluster_index):
        b = Bunch()
        b.correlogram = acg
        b.firing_rate = None
        m = acg.max() if acg.size and acg.max() > 0 else 1.0
        b.data_bounds = (0, 0, len(acg), m)
        b.pair_index = (0, col)
        b.color = selected_cluster_color(cluster_index, 1)
        b.jitter_ccg = None
        return b

    def _ccg_bunch(self, ccg_n, ccg_corr, col):
        b = Bunch()
        b.correlogram = ccg_n
        b.jitter_ccg = ccg_corr
        b.firing_rate = None
        # Share a y-scale across both curves so the overlay is comparable.
        m = max(ccg_n.max() if ccg_n.size else 0.0,
                ccg_corr.max() if ccg_corr.size else 0.0)
        if m <= 0:
            m = 1.0
        b.data_bounds = (0, 0, len(ccg_n), m)
        b.pair_index = (0, col)
        # Cross-correlogram color (desaturated), like the parent off-diagonal.
        c = selected_cluster_color(0, 1)
        b.color = add_alpha(_override_hsv(c[:3], s=.1, v=1))
        return b

    # -- plotting -------------------------------------------------------------

    def _plot_pair(self, bunch):
        # Raw histogram (parent path).
        self.correlogram_visual.add_batch_data(
            hist=bunch.correlogram, color=bunch.color,
            ylim=bunch.data_bounds[3], box_index=bunch.pair_index)

        # Jitter-corrected CCG overlaid as a line, on the SAME y-scale.
        if bunch.get('jitter_ccg', None) is not None:
            jc = bunch.jitter_ccg
            n = len(jc)
            x = np.arange(n) + 0.5      # bin centers
            pos = np.column_stack([x[:-1], jc[:-1], x[1:], jc[1:]])
            white = (1., 1., 1., 1.)
            self.line_visual.add_batch_data(
                pos=pos, color=white, data_bounds=bunch.data_bounds,
                box_index=bunch.pair_index)

        # Refractory-period markers (two vertical lines), as in the parent.
        xrp0 = round((self.window_size * .5 - self.refractory_period) / self.bin_size)
        xrp1 = round((self.window_size * .5 + self.refractory_period) / self.bin_size) + 1
        ylim = bunch.data_bounds[3]
        gray = (.25, .25, .25, 1.)
        pos = np.array([[xrp0, 0, xrp0, ylim], [xrp1, 0, xrp1, ylim]])
        self.line_visual.add_batch_data(
            pos=pos, color=gray, data_bounds=bunch.data_bounds,
            box_index=bunch.pair_index)

    def _plot_labels(self):
        if len(self.cluster_ids) != 2:
            return
        a, b = self.cluster_ids
        labels = ['ACG %d' % a, 'CCG %d-%d' % (a, b), 'ACG %d' % b]
        for col, txt in enumerate(labels):
            self.text_visual.add_batch_data(
                pos=[0, -1], text=txt, anchor=[0, -1.25],
                data_bounds=None, box_index=(0, col))

    def _plot_hint(self):
        """Graceful message when not exactly 2 clusters are selected."""
        self.canvas.grid.shape = (1, 1)
        self.correlogram_visual.reset_batch()
        self.line_visual.reset_batch()
        self.text_visual.reset_batch()
        n = len(self.cluster_ids)
        msg = ('Select exactly 2 clusters for the pair triptych '
               '(ACG | CCG+jitter | ACG). Currently %d selected.' % n)
        self.text_visual.add_batch_data(
            pos=[0, 0], text=msg, anchor=[0, 0], data_bounds=None,
            box_index=(0, 0))
        self.canvas.update_visual(self.correlogram_visual)
        self.canvas.update_visual(self.line_visual)
        self.canvas.update_visual(self.text_visual)
        self.canvas.update()

    def plot(self, **kwargs):
        """Update the view: 1x3 triptych, or a hint if not exactly 2 clusters."""
        if len(self.cluster_ids) != 2:
            self._plot_hint()
            return

        self.canvas.grid.shape = (1, 3)
        bunchs = self.get_clusters_data()

        self.correlogram_visual.reset_batch()
        self.line_visual.reset_batch()
        self.text_visual.reset_batch()

        for bunch in bunchs:
            self._plot_pair(bunch)
        self._plot_labels()

        self.canvas.update_visual(self.correlogram_visual)
        self.canvas.update_visual(self.line_visual)
        self.canvas.update_visual(self.text_visual)
        self.canvas.update()


# -----------------------------------------------------------------------------
# Plugin: register the view creator on the controller
# -----------------------------------------------------------------------------

class PairCCGPlugin(IPlugin):
    def attach_to_controller(self, controller):
        def create_pair_correlogram_view():
            # Pull live spike_clusters from the supervisor when available so the
            # view honors in-GUI merges/splits; fall back to the static model.
            def spike_clusters_func():
                sup = getattr(controller, 'supervisor', None)
                if sup is not None and getattr(sup, 'clustering', None) is not None:
                    return sup.clustering.spike_clusters
                return controller.model.spike_clusters

            return PairCorrelogramView(
                model=controller.model,
                spike_clusters_func=spike_clusters_func,
            )

        controller.view_creator['PairCorrelogramView'] = create_pair_correlogram_view

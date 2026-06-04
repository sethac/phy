"""phy raw-data filter plugin: CAR + highpass for trace/waveform views - v0.1

phy's built-in default filter is only a highpass; it does NOT do common
average referencing (CAR), so raw Neuropixels traces look noisy with
common-mode that KS4/CatGT removed before sorting. This plugin registers
extra filters so you can cycle them in the GUI with alt+r:

    raw  ->  high_pass (phy default)  ->  car  ->  car_highpass

- car          : subtract the per-sample MEDIAN across channels (flat CAR).
- car_highpass : highpass (300 Hz) THEN CAR -> closest to what the sorter saw.

These are DISPLAY-ONLY. They never touch the .ap.bin or the sort.

NOTE on axis: phy calls the filter with axis=0 for the trace view (array is
(time, channels)) and axis=1 for the waveform view (array is
(n_spikes, time, channels) or (time, channels)). CAR must always subtract
across the CHANNEL axis, which is the LAST axis in every case. So we operate
on axis=-1 for the median, independent of the `axis` phy passes for the
temporal filter.

Install: this file lives in C:/Users/sa4257/.phy/plugins/ (already on
phy's Plugins.dirs). Activate by adding 'CARFilterPlugin' to the plugins
list for the TemplateGUI in ~/.phy/phy_config.py (this plugin appends itself
when phy imports it via attach_to_controller).
"""
import numpy as np
from scipy.signal import butter, lfilter
from phy import IPlugin


def _car(arr):
    """Subtract the median across channels (last axis). Display-only CAR."""
    return arr - np.median(arr, axis=-1, keepdims=True)


class CARFilterPlugin(IPlugin):
    def attach_to_controller(self, controller):
        sr = controller.model.sample_rate
        # 300 Hz highpass (spike band), 3rd order Butterworth - matches the
        # spirit of phy's default but at the spike-band corner.
        b, a = butter(3, 300.0 / sr * 2.0, "high")

        def _highpass(arr, axis=0):
            # zero-phase: filter, flip, filter, flip (same trick phy uses)
            arr = lfilter(b, a, arr, axis=axis)
            arr = np.flip(arr, axis=axis)
            arr = lfilter(b, a, arr, axis=axis)
            arr = np.flip(arr, axis=axis)
            return arr

        @controller.raw_data_filter.add_filter
        def car(arr, axis=0):
            return _car(arr)

        @controller.raw_data_filter.add_filter
        def car_highpass(arr, axis=0):
            return _car(_highpass(arr, axis=axis))

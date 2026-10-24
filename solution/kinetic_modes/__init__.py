"""Finite-window kinetic modes and pole response."""
import numpy as np

def _real_vector(value):
    raw=np.asarray(value)
    if raw.ndim!=1 or np.iscomplexobj(raw):raise ValueError('Finite real vectors required')
    try:result=np.asarray(value,dtype=float)
    except (TypeError,ValueError) as error:raise ValueError('Finite real vectors required') from error
    if not np.all(np.isfinite(result)):raise ValueError('Nonfinite input')
    return result

def kinetic_spectrum(fractions, thermal_widths, drift_frequencies, window, sample_times):
    try:
        f,a,d=(_real_vector(v) for v in (fractions,thermal_widths,drift_frequencies))
        box,times=_real_vector(window),_real_vector(sample_times)
        if not len(f) or a.shape!=f.shape or d.shape!=f.shape:raise ValueError('Matching nonempty population vectors required')
        if np.any(f<=0) or np.any(a<=0) or abs(f.sum()-1)>1e-12:raise ValueError('Positive normalized populations required')
        if box.shape!=(4,) or box[0]>=box[1] or box[2]>=box[3]:raise ValueError('Nonempty open rectangle required')
        if not len(times) or np.any(times<0):raise ValueError('Nonempty nonnegative times required')
        from .core import Dielectric
        return Dielectric(f,a,d).spectrum(box,times)
    except (TypeError,IndexError,KeyError) as error:
        raise ValueError('Malformed kinetic-spectrum request') from error

"""Landau dielectric, finite-window mode isolation and pole residues."""
import numpy as np
from plasmapy.dispersion import plasma_dispersion_func as Z


class Dielectric:
    def __init__(self, fractions, thermal_widths, drift_frequencies):
        self.fractions = np.asarray(fractions, float)
        self.widths = np.asarray(thermal_widths, float)
        self.drifts = np.asarray(drift_frequencies, float)
        if self.fractions.ndim != 1 or self.widths.shape != self.fractions.shape or self.drifts.shape != self.fractions.shape:
            raise ValueError('Matching Maxwellian population vectors required')
        if np.any(self.fractions <= 0) or np.any(self.widths <= 0) or abs(self.fractions.sum()-1) > 1e-12:
            raise ValueError('Positive normalized fractions and thermal widths required')
        self.weights = self.fractions/self.widths**2
        self.log_weights = np.log(self.fractions)-2*np.log(self.widths)

    def scaled(self, frequency):
        z = np.asarray(frequency, complex)
        xi = (z[..., None]-self.drifts)/(np.sqrt(2)*self.widths)
        upper = xi.imag >= 0
        argument = np.where(upper, xi, -xi)
        native = Z(argument)
        derivative = -2*(1+argument*native)
        core = np.where(upper, 1+xi*native, 1-xi*native)
        core_derivative = np.where(upper, native+xi*derivative, -native+xi*derivative)
        exponential = -xi**2
        weighted_exponential = exponential+self.log_weights
        shift = np.maximum(0., np.max(np.where(upper, -np.inf, weighted_exponential.real), axis=-1))
        factor = np.exp(-shift)
        reflected = np.exp(np.where(upper, -np.inf+0j, weighted_exponential-shift[..., None]))
        value = (1+np.sum(self.weights*core, axis=-1))*factor + np.sum(2j*np.sqrt(np.pi)*xi*reflected, axis=-1)
        dz = 1/(np.sqrt(2)*self.widths)
        gradient = np.sum(self.weights*core_derivative*dz, axis=-1)*factor + np.sum(2j*np.sqrt(np.pi)*dz*(1-2*xi**2)*reflected, axis=-1)
        return value, gradient, shift

    def direct(self, frequency):
        z = np.asarray(frequency, complex)
        xi = (z[..., None]-self.drifts)/(np.sqrt(2)*self.widths)
        return 1+np.sum(self.weights*(1+xi*Z(xi)), axis=-1)

    @staticmethod
    def contour(box, points):
        left, right, bottom, top = box
        vertices = [left+1j*bottom, right+1j*bottom, right+1j*top, left+1j*top, left+1j*bottom]
        edges = [np.linspace(a, b, points, endpoint=False) for a, b in zip(vertices[:-1], vertices[1:])]
        return np.concatenate(edges+[np.array([vertices[0]])])

    def count(self, box):
        contour = self.contour(box, 64)
        left, right = contour[:-1], contour[1:]
        angle = 0.
        for depth in range(36):
            if not len(left):
                count = int(np.rint(angle/(2*np.pi)))
                if count < 0 or abs(angle/(2*np.pi)-count)>1e-7:
                    raise RuntimeError('Inconsistent kinetic winding count')
                return count
            middle = (left+right)/2
            points = np.stack([left,middle,right],axis=1)
            value,derivative,_ = self.scaled(points)
            finite = np.all(np.isfinite(value)&np.isfinite(derivative)&(np.abs(value)>1e-250),axis=1)
            with np.errstate(invalid='ignore',divide='ignore'):
                unit=value/np.abs(value)
                increments=np.angle(unit[:,1:]*np.conjugate(unit[:,:-1]))
                resolution=np.max(np.abs(derivative/value),axis=1)*np.abs(right-left)
            good=finite&(resolution<.45)&(np.max(np.abs(increments),axis=1)<.45)
            angle+=float(np.sum(increments[good]))
            bad=~good
            left,right=np.stack([left[bad],middle[bad]],axis=1).ravel(),np.stack([middle[bad],right[bad]],axis=1).ravel()
            if len(left)>262144:
                raise RuntimeError('Kinetic boundary could not be resolved')
        raise RuntimeError('Kinetic boundary refinement exhausted')

    def polish(self, box):
        left, right, bottom, top = box
        guesses = [complex((left+right)/2, (bottom+top)/2)]
        guesses += [complex(left+a*(right-left), bottom+b*(top-bottom)) for a, b in ((.25, .25), (.75, .25), (.25, .75), (.75, .75))]
        for initial in guesses:
            z = initial
            for iteration in range(80):
                value, derivative, shift = self.scaled(z)
                if not np.isfinite(value) or not np.isfinite(derivative) or derivative == 0:
                    break
                correction = value/derivative
                if abs(correction) < 2e-12:
                    if left < z.real < right and bottom < z.imag < top:
                        return complex(z)
                    break
                step = 1.
                accepted = False
                for _ in range(24):
                    candidate = z-step*correction
                    if left < candidate.real < right and bottom < candidate.imag < top:
                        new_value, new_derivative, _ = self.scaled(candidate)
                        if new_derivative != 0 and abs(new_value/new_derivative) < abs(correction)*1.2:
                            z, accepted = candidate, True
                            break
                    step *= .5
                if not accepted:
                    break
        return None

    def spectrum(self, window, times):
        box = tuple(map(float, window))
        total = self.count(box)
        pending, roots = [(box, total, 0)], []
        while pending:
            region, count, depth = pending.pop()
            if not count:
                continue
            if count == 1:
                root = self.polish(region)
                if root is not None:
                    roots.append(root)
                    continue
            if depth >= 40:
                raise RuntimeError('Kinetic zeros could not be isolated')
            left, right, bottom, top = region
            success = False
            for ratio in (.47391, .52913, .41627, .58173):
                if right-left >= top-bottom:
                    split = left+ratio*(right-left)
                    children = [(left, split, bottom, top), (split, right, bottom, top)]
                else:
                    split = bottom+ratio*(top-bottom)
                    children = [(left, right, bottom, split), (left, right, split, top)]
                try:
                    counts = [self.count(child) for child in children]
                except RuntimeError:
                    continue
                if sum(counts) == count:
                    pending.extend((child, child_count, depth+1) for child, child_count in zip(children, counts))
                    success = True
                    break
            if not success:
                raise RuntimeError('Inconsistent kinetic root counts at a split')
        roots = np.array(sorted(roots, key=lambda z: (z.real, z.imag)), dtype=complex)
        if len(roots) != total:
            raise RuntimeError('Incomplete finite-window kinetic spectrum')
        if len(roots) > 1 and np.min(np.abs(roots[:, None]-roots[None, :])+np.eye(len(roots))*1e100) < 1e-8:
            raise RuntimeError('Duplicate kinetic modes')
        value, derivative, shift = self.scaled(roots)
        residues = np.exp(-shift)/derivative
        times = np.asarray(times, float)
        propagator = np.exp(-1j*times[:, None]*roots[None, :])
        response = propagator@residues
        response_derivative = propagator@(-1j*roots*residues)
        return {'roots': roots, 'residues': residues, 'response': response, 'response_derivative': response_derivative}

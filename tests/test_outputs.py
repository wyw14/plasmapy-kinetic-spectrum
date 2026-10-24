"""Independent complete spectra, residues, and finite-window pole responses."""
from pathlib import Path
import copy
import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment
from candidate import call
ROOT=Path(__file__).parent/'data'
NAMES=['single-warm','counterstream','cold-window','three-population']

def load(name):
    with np.load(ROOT/(name+'-input.npz'),allow_pickle=False) as z:
        args=[z[k].copy() for k in ['fractions','thermal_widths','drift_frequencies','window','sample_times']]
    with np.load(ROOT/(name+'-expected.npz'),allow_pickle=False) as z:expected={k:z[k].copy() for k in z.files}
    return args,expected

def pair(actual,expected):
    assert isinstance(actual,dict) and set(actual)==set(expected)
    roots=np.asarray(actual['roots'],dtype=complex);count=len(expected['roots'])
    assert roots.shape==(count,)
    for name,value in actual.items():
        value=np.asarray(value);assert value.shape==expected[name].shape and value.dtype.kind in 'fciu' and np.all(np.isfinite(value))
    if not count:return actual
    i,j=linear_sum_assignment(np.abs(roots[:,None]-expected['roots'][None,:]))
    order=i[np.argsort(j)]
    return {name:np.asarray(value)[order] if name in ['roots','residues'] else np.asarray(value) for name,value in actual.items()}

class Deferred:
    def __init__(self,args,expected):self.args=args;self.expected=expected;self.result=None;self.error=None
    def get(self):
        if self.result is None and self.error is None:
            try:self.result=pair(call('kinetic_spectrum',*self.args),self.expected)
            except Exception as error:self.error=error
        if self.error is not None:raise self.error
        return self.result

@pytest.fixture(scope='module',params=NAMES)
def state(request):
    args,expected=load(request.param)
    return args,expected,Deferred(args,expected)

def test_complete_finite_shapes_and_mode_count(state):
    args,expected,actual=state;actual.get()

@pytest.mark.parametrize('field',['roots','residues','response','response_derivative'])
def test_independent_kinetic_spectrum(state,field):
    args,expected,actual=state
    np.testing.assert_allclose(actual.get()[field],expected[field],atol=2e-8,rtol=2e-8 if field=='roots' else 2e-7)

def compare(args,expected):
    actual=pair(call('kinetic_spectrum',*args),expected)
    for name in expected:np.testing.assert_allclose(actual[name],expected[name],atol=2e-8,rtol=2e-8 if name=='roots' else 2e-7,err_msg=name)

def test_galilean_frequency_shift_preserves_modal_content():
    args,expected=load('cold-window');shift=.137
    args[2]+=shift;args[3][:2]+=shift
    phase=np.exp(-1j*shift*args[4]);expected['roots']+=shift
    expected['response_derivative']=phase*(expected['response_derivative']-1j*shift*expected['response'])
    expected['response']*=phase
    compare(args,expected)

def test_velocity_reflection_transforms_roots_and_residue_sign():
    args,expected=load('three-population');args[2]=-args[2];args[3][:2]=-args[3][1::-1]
    for name in expected:expected[name]=-np.conjugate(expected[name])
    compare(args,expected)

def test_population_permutation_preserves_complete_spectrum():
    args,expected=load('three-population');order=[2,0,1]
    for i in range(3):args[i]=args[i][order]
    compare(args,expected)

def test_identical_population_split_preserves_dielectric():
    args,expected=load('single-warm');args[0]=np.array([.37,.63]);args[1]=np.repeat(args[1],2);args[2]=np.repeat(args[2],2)
    compare(args,expected)

def test_unsorted_repeated_times_and_late_response():
    args,expected=load('counterstream');args[4]=np.array([10.,0.,.7,.7,3.])
    propagator=np.exp(-1j*args[4][:,None]*expected['roots'][None,:])
    expected['response']=propagator@expected['residues'];expected['response_derivative']=propagator@(-1j*expected['roots']*expected['residues'])
    compare(args,expected)

def test_empty_upper_half_window_has_no_pole_contribution():
    args,_=load('single-warm');args[3]=np.array([.1,.2,.4,.6])
    expected={'roots':np.empty(0,dtype=complex),'residues':np.empty(0,dtype=complex),'response':np.zeros(len(args[4]),complex),'response_derivative':np.zeros(len(args[4]),complex)}
    compare(args,expected)

def test_complete_spectrum_with_a_mode_close_to_window_edge():
    args,expected=load('cold-window')
    args[3][1]=np.max(expected['roots'].real)+2e-6
    compare(args,expected)

def test_trace_population_does_not_destroy_contour_resolution():
    args,expected=load('trace-population')
    compare(args,expected)

def test_repeated_calls_do_not_change_input_or_result():
    args,expected=load('single-warm');compare(args,expected);compare(args,expected)

@pytest.mark.parametrize('variant',range(10))
def test_invalid_input_has_value_error(variant):
    args,_=load('single-warm')
    if variant==0:args[0]=np.array([])
    elif variant==1:args[1]=np.array([.2,.3])
    elif variant==2:args[1]=np.array([0.])
    elif variant==3:args[0]=np.array([.9])
    elif variant==4:args[2]=np.array([np.nan])
    elif variant==5:args[3]=np.array([1.,-1.,-1.,1.])
    elif variant==6:args[4]=np.array([])
    elif variant==7:args[4]=np.array([-1.])
    elif variant==8:args[3]=np.array([-1.,1.,np.inf,2.])
    else:args[0]=np.array([1.+.1j])
    with pytest.raises(ValueError):call('kinetic_spectrum',*args)

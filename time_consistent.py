import numpy
from numpy import exp
import scipy.optimize
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline, splev, splrep


calibration = dict(
    beta=0.8/(0.8+0.15),
    a=0.8,
    c=0.15,
    estar=-0.0,
    Rbar=0.5,
    min_f=0,
    kappa=1.0,
    N=40,
    zbar=0.1,
    lam=0.9, # probability that crisis continues
    model='optimal',
    theta=0.0,
    p=1.0,
    b=0.0
)


def residuals(rvec, e, f, calib):

    a = calib['a']
    c = calib['c']
    estar = calib['estar']
    beta = calib['beta']
    zbar = calib['zbar']
    theta = calib['theta']
    p = calib['p']
    b = calib['b']

    fun_e = splrep(rvec, e, k=5)
    fun_f = splrep(rvec, f, k=5)

    R_f = rvec - f

    f_f = splev(R_f, fun_f, der=0)
    e_f = splev(R_f, fun_e, der=0)

    d_f_f = splev( R_f, fun_f, der=1)
    d_e_f = splev( R_f, fun_e, der=1)

    psi = (1-b)*(f-f**2*theta/2)
    psi_f = (1-b)*(f-f_f**2*theta/2)

    d_psi = (1-b)*(1-f*theta)
    d_psi_f = (1-b)*(1-f_f*theta)

    cond_1 = (e - estar)*(d_psi+a*p*d_e_f) - beta*(e_f - estar)*d_psi_f*p
    cond_2 = e - p*a/(a+c)*e_f + 1/(a+c)*(psi_f-zbar)

    return [cond_1, cond_2]

def make_init(rvec, calib):

    init = numpy.concatenate( [rvec[None,:], rvec[None,:]], axis=0)
    beta = calib['beta']
    zbar = calib['zbar']
    a = calib['a']
    c = calib['c']
    p = calib['p']

    xstar = 1.0/(a+c-a*p)*zbar
    r0 = xstar/(a+c)
    r2vec = numpy.concatenate([rvec-rvec.max(), rvec])

    e =  numpy.maximum(xstar+(beta*p-1)/(p*a)*r2vec, 0)
    f = numpy.minimum( (1-beta)/beta*c/a*r2vec, zbar )

    N = len(rvec)

    evec = e[N:]
    fvec = f[N:]
    init = numpy.row_stack([evec,fvec*0.999])
    return init




def solve(initial_guess=None, max_R=8, N=20, **cc):

    calib = calibration.copy()
    calib.update(cc)

    rvec = numpy.linspace(0.0001,max_R,N)
    # max_f = numpy.minimum(rvec*1.1)
    max_f = rvec*1.1

    def from_xi(u):
        uu = u.copy().reshape((2,-1))
        e = uu[0,:]
        xx = uu[1,:]
        f = max_f*(1+numpy.tanh(xx))/2
        uu[1,:] = f
        return uu

    def to_xi(u):
        uu = u.copy()
        f = uu[1,:]
        uu[1,:] = numpy.arctanh( 2*f/max_f-1 )
        return uu
    def fobj(u):
        uu = u.reshape((2,-1))
        e = uu[0,:]
        xx = uu[1,:]

        f = max_f*(1+numpy.tanh(xx))/2
        res = residuals(rvec, e, f, calib)
        return numpy.concatenate(res)

    if initial_guess is None:
        init = make_init(rvec, calib)
    elif isinstance(initial_guess,tuple):
        fun_e, fun_f = initial_guess[3]
        init = numpy.row_stack([
            fun_e(rvec),
            fun_f(rvec)
        ])
    else:
        init = initial_guess

    res = scipy.optimize.root(
        fobj,
        to_xi(init),
        method='lm',
        options={'ftol':1e-10, 'xtol':1e-10}
    )
    x = from_xi(res.x)

    spl_e = splrep(rvec, x[0,:], k=5)
    spl_f = splrep(rvec, x[1,:], k=5)
    fun_e = lambda x: splev(x, spl_e )
    fun_f = lambda x: splev(x, spl_f )

    return rvec, x, res.x, [fun_e, fun_f]

def simulate(r0, drs, N):
    dr_e, dr_f = drs[3]
    import pandas
    vals = []
    R = r0
    for i in range(N):
        e = dr_e(R)
        f = dr_f(R)
        vals.append([R,e,f])
        R = R - f
    return pandas.DataFrame(vals, columns=['R','e','f'])

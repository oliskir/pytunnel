import numpy as np
import scipy.sparse as ss
import scipy.sparse.linalg as ssl
from tqdm import tqdm
import pytunnel.units as units
from pytunnel.output import peak_position


def compute_matrices(J, dt, dx, V):

    # array of 1's
    o = np.ones(J, complex)

    # alpha from (3.11)
    alp = 1j * dt / (2 * dx**2) * o

    # xi and gam from (3.11)
    xi = o + 1j * dt / 2 * (2 / (dx**2) * o + V)  
    gam = o - 1j * dt / 2 * (2 / (dx**2) * o + V)

    # positions of the vectors in the matrix
    diags = np.array([-1, 0, 1])  
    vecs1 = np.array([-alp, xi, -alp])
    vecs2 = np.array([alp, gam, alp])

    # create tridiagonal sparse matrices
    U1 = ss.spdiags(vecs1, diags, J, J) 
    U2 = ss.spdiags(vecs2, diags, J, J) 

    # convert to different sparse format needed for further calculation
    U1 = U1.tocsc() 
    U2 = U2.tocsc() 

    # compute LU-decomposition of U1
    LU = ssl.splu(U1) 

    return LU, U2


def plane_wave_square_barrier_transmission(E, V, w):

    # energy in MeV
    # potential in MeV
    # width in fm

    if E == V:
        E = V * 1.00001

    E = E / units.MeV
    V = V / units.MeV
    w = w / units.fm

    k = np.sqrt(E)
    kap = np.sqrt(np.abs(V-E))
    rho = kap / k

    if E < V:
        T = 4*rho**2 / ((1+rho**2)**2*np.sinh(kap*w)**2 + 4*rho**2)
    else:
        T = 4*rho**2 / ((1-rho**2)**2*np.sin(kap*w)**2 + 4*rho**2)

    return T


def initial_wf(J, dx, x0, sig0, k0, xi0, j1V=None):

    x = np.arange(J, dtype=np.float64)
    x *= dx
    x += x0
    psi0 = np.power(np.pi*sig0**2, -0.25) * np.exp(1j * k0 * x - (x-xi0)**2/(2*sig0**2))

    # ensure that wave function is non-zero only to the left of the barrier
    if j1V is not None:
        psi0[j1V:] = np.zeros(J-j1V, complex)
        norm = psi0.real**2 + psi0.imag**2
        norm = np.sqrt(np.sum(norm) * dx)
        psi0 = psi0 / norm
        #print('\n w.f. cut off for x > {0:.2f}\n'.format(x[j1V]))

    return psi0, x


def solveTDSE(psi0, N, dx, dt, num_frames=0, barrier=None, logger=None, plotter=None, progress_bar=True, termin_thres=0):

    J = len(psi0)

    # psi(x,t=0)
    PSI = psi0

    # loop over time-steps
    nf = max(1, int(N / max(1,num_frames)))
    for n in tqdm(range(0, N), disable = not progress_bar):

        # update time
        time = n * dt

        # barrier shape
        if barrier is None: V = np.zeros(len(PSI))
        else: V = barrier.get_array(time)

        # log data
        if num_frames > 0 and n%nf == 0:

            # w.f. probability density
            prob = PSI
            prob = prob.real**2 + prob.imag**2

            # log data
            if logger is not None:
                T, A = logger.log(time=time, psi=PSI, V=V)

            # plot
            if plotter is not None:
                img_id = int(n/nf)
                plotter.plot(prob=prob, time=time, id=img_id)

            # early termination
            if termin_thres > 0:
                if T > 0 and A/T < termin_thres: 
                    print(f'Computation terminated at A/T = {A/T:.3E} < {termin_thres:.3E}')
                    break

        # store current w.f.
        if logger is not None:
            logger.update_wf(PSI)

        # compute matrices
        if n == 0 or (barrier is not None and barrier.has_changed(time, time-dt)):
            LU, U2 = compute_matrices(J, dt, dx, V)

        # right hand side of eq. (3.9)
        b = U2.dot(PSI)
        
        # solve system of equations for each time step
        PSI = LU.solve(b) 

    if logger is not None:
        return logger, PSI
    else:
        return PSI


def bin_size(sig, k, velocity, accuracy):
    ''' Determine bin size necessary to achieve specified numerical accuracy
    '''
    domain = 10 * sig
    xx, yy = [0], [0]
    n = 0
    while len(xx) < 2 or yy[-1] < accuracy:
        dx = np.power(10, n)
        n += 1
        J = int(domain / dx)
        psi0, x = initial_wf(J, dx, x0=-domain/2, sig0=sig, k0=k, xi0=0)
        y = num_rel_err(x=x, psi0=psi0, dx=dx, dt=dx/velocity, velocity=velocity)
        xx.append(dx)
        yy.append(y)

    xx_neg = [-x for x in xx[-1:0:-1]]
    yy_neg = [y for y in yy[-1:0:-1]]
    xx = xx_neg + xx
    yy = yy_neg + yy
    zz = np.polyfit(xx, yy, 2)
    poly = np.poly1d(zz)
    r = (poly - accuracy).roots
    return np.max(r)


def num_rel_err(x, psi0, dx, dt, velocity, N=10):
    ''' Estimate relative numerical error in propagation of wave packet
    '''
    psi = solveTDSE(psi0=psi0, N=N, dx=dx, dt=dt, progress_bar=False)
    prob0 = psi0.real**2 + psi0.imag**2
    prob = psi.real**2 + psi.imag**2
    x_start = peak_position(x, prob0)
    x_peak = peak_position(x, prob)
    x_free = x_start + velocity * N * dt
    x_err_rel = abs((x_peak - x_free) / (x_free - x_start))
    return x_err_rel

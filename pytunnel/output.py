import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pytunnel.units as units


def peak_position(x, p, n=10):
    ''' Determine peak position
    '''
    assert len(x)==len(p), 'x and p must have same length'
    i = np.argmax(p)
    if i==0 or i==len(x)-1: return x[i]
    xx = x[max(0,i-n):min(len(x),i+n+1)]
    pp = p[max(0,i-n):min(len(x),i+n+1)]
    zz = np.polyfit(xx, pp, 2)
    poly = np.poly1d(zz)
    deriv = poly.deriv()
    x_peak = deriv.roots[0]
    return x_peak

def pos_expect_val(x, p):
    ''' Determine position expectation value
    '''
    assert len(x)==len(p), 'x and p must have same length'
    p_sum = np.sum(p)
    xp_sum = np.sum(x * p)
    if p_sum == 0: return 0
    else: return xp_sum / p_sum


def energy_expect_val(y, dydt, V):
    ''' Determine kinetic energy as K = int{ i.hbar.psi*.dpsi/dt - psi*.V.psi }
        Returns e_kin and e_tot.
    '''
    yc = np.conjugate(y)
    Esum = np.sum(1j * yc * dydt)
    Vsum = np.sum(yc * V * y)
    Ksum = Esum - Vsum
    p = y.real**2 + y.imag**2
    psum = np.sum(p)
    if psum > 0: return np.absolute(Ksum / psum), np.absolute(Esum / psum)
    else: return 0, 0

def momentum_expect_val(y, dydx):
    ''' Determine kinetic energy as p = -i.hbar.int{ psi*.dpsi/dx }
    '''
    yc = np.conjugate(y)
    mom = np.sum(-1j * yc * dydx)
    p = y.real**2 + y.imag**2
    psum = np.sum(p)
    if psum > 0: return np.absolute(mom / psum)
    else: return 0


def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

def clean_dir(path):
    files = glob.glob(path+'*')
    for f in files:
        if f is not path and os.path.isfile(f):
            os.remove(f)

def prep_output_folder(output_dir, save_frames):

    if output_dir[-1] != '/':
        output_dir += '/'

    ensure_dir(output_dir)
    clean_dir(output_dir)

    if save_frames:
        subdir = output_dir+'frames/'
        ensure_dir(subdir)
        clean_dir(subdir)

    return output_dir


class Logger():

    def __init__(self, barrier, energy, dx, dt, x_axis, x_start, t_enter, velocity_left, velocity_right, output_dir):

        self.barrier = barrier
        self.energy = energy
        self.dx = dx
        self.dt = dt
        self.x_axis = x_axis
        self.x_start = x_start
        self.t_enter = t_enter
        self.velocity_left = velocity_left
        self.velocity_right = velocity_right

        self.log_file_path = os.path.join(output_dir,'log.csv')
        self.log_file = open(self.log_file_path, 'w')
        self.log_file.write('time (s),R,A,T,pos_IR (fm),pos_T (fm),delay (s),pos_IR_ev (fm),pos_T_ev (fm),delay_ev (s),x_free (fm),ekin_IR (MeV),ekin_T (MeV),beta_IR,beta_T,etot (MeV)')
        self.log_file.close()

        self.time = list()
        self.refl = list()
        self.transm = list()
        self.absorbed = list()
        self.class_pos = list()

        self.psi_prev = None

    def update_wf(self, psi):
        self.psi_prev = psi

    def log(self, time, psi, V):

        prob = psi.real**2 + psi.imag**2

        dx = self.dx
        b1 = self.barrier.left_bin(t=time, E=self.energy)
        b2 = self.barrier.right_bin(t=time, E=self.energy) + 1

        R = np.sum(prob[:b1]) * dx
        T = np.sum(prob[b2:]) * dx
        A = np.sum(prob[b1:b2]) * dx

        # incident/reflected peak position and expectation value
        xIR = peak_position(self.x_axis[:b1], prob[:b1])
        xIR_ev = pos_expect_val(self.x_axis[:b1], prob[:b1])

        # transmitted peak position and expectation value
        xT = peak_position(self.x_axis[b2:], prob[b2:])
        xT_ev = pos_expect_val(self.x_axis[b2:], prob[b2:])

        # delay wrt free propagation
        time_s = time * units.sec
        x_free = self.x_start + self.velocity_left * min(time_s, self.t_enter) + self.velocity_right * max(0, time_s - self.t_enter)
        delay = (x_free - xT) / self.velocity_right
        delay_ev = (x_free - xT_ev) / self.velocity_right

        # kinetic energy, K = int{ i.hbar.psi*.dpsi/dt - psi*.V.psi }
        if self.psi_prev is not None:
            dpsidt = (psi - self.psi_prev) / self.dt
            KIR, EIR = energy_expect_val(psi[:b1], dpsidt[:b1], V[:b1]) * units.MeV
            KT, ET   = energy_expect_val(psi[b2:], dpsidt[b2:], V[b2:]) * units.MeV
        else:
            KIR, EIR = 0, 0
            KT, ET   = 0, 0

        # total energy (incident/reflected + transmitted)
        Etot = R * EIR + T * ET
        
        # momentum, p = -i.hbar.int{ psi*.dpsi/dx }       
        # v = p / m,  m=0.5 
        dpsidx = np.gradient(psi) / self.dx
        vIR = momentum_expect_val(psi[:b1], dpsidx[:b1]) / 0.5 * units.fm / units.sec / units.c
        vT  = momentum_expect_val(psi[b2:], dpsidx[b2:]) / 0.5 * units.fm / units.sec / units.c

        self.time.append(time_s)
        self.refl.append(R)
        self.transm.append(T)
        self.absorbed.append(A)
        self.class_pos.append(x_free)

        self.log_file = open(self.log_file_path, 'a')
        self.log_file.write(f'\n{time_s:.4E},{R:.4E},{A:.4E},{T:.4E},{xIR:.4E},{xT:.4E},{delay:.4E},{xIR_ev:.4E},{xT_ev:.4E},{delay_ev:.4E},{x_free:.4E},{KIR:.4E},{KT:.4E},{vIR:.4E},{vT:.4E},{Etot:.4E}')
        self.log_file.close()

        self.psi_prev = psi

        return T,A


class Plotter():

    def __init__(self, energy, barrier, logger, x_axis, bin1, bin2, ymin, output_dir):

        self.energy = energy
        self.barrier = barrier
        self.logger = logger
        self.bin1 = bin1
        self.bin2 = bin2
        self.ymin = ymin
        self.pmax = None
        self.x = x_axis[bin1:bin2]
        self.output_dir = output_dir

    def plot(self, time, prob, id):

        # wave function squared max value
        if self.pmax is None:
            self.pmax = np.max(prob)

        # x axis range
        b1 = self.bin1
        b2 = self.bin2

        # lin or log y axis
        logy = self.ymin > 0

        # barrier
        V = self.barrier.get_array(t=time).real[b1:b2] * units.MeV
        Vmax = self.barrier.max() * units.MeV

        # normalized w.f. probability density
###        p = prob[b1:b2] / self.pmax * self.energy
        p = prob[b1:b2] / self.pmax * 0.8 * self.barrier.max() * units.MeV

        # y axis range
        ymin = self.ymin
        ymax = 1.2 * Vmax #1.1 * max(self.energy, Vmax)

        R = self.logger.refl[-1]
        T = self.logger.transm[-1]
        A = self.logger.absorbed[-1]
        psum = R + T + A

        plot_wf(id, self.x, p, V, logy, psum, R, T, A, time*units.sec, ymin, ymax, self.output_dir)


def plot_wf(id, x, p, V, ylog, psum, Refl, Tran, diff, time, ymin, ymax, outdir):

    xmin = np.min(x)
    xmax = np.max(x)

    plt.plot(x, p)
    plt.plot(x, V)    

    plt.xlabel('x (fm)')
    plt.ylabel('V (MeV), $|\psi|^2$')

    if ylog:
        plt.yscale('log')

    plt.ylim(ymin, ymax)

    plt.title('Wave packet')
    plt.grid(True)

    xt = xmin + 0.74*(xmax-xmin)
    yt = 0.7 * ymax
    if ylog:
        yt = np.power(10., 0.7*(np.log10(ymax)-np.log10(ymin))+np.log10(ymin))

    txt = 't = {4:.4E} s\nR = {1:.2f}\nT = {2:.2E}\n1-R-T = {3:.2E}\n$\int |\psi|^2 = ${0:.2f}'.format(psum, Refl, Tran, diff, time)
    plt.text(xt, yt, txt)

    plt.savefig(outdir+'frames/img{0}.png'.format(id))
    plt.clf()


def plot_evolution(logger, ymin, output_dir):

    class_pos = logger.class_pos
    time = logger.time
    refl = logger.refl
    transm = logger.transm
    absorbed = logger.absorbed

    plt.plot(class_pos, refl, label='R')    
    plt.plot(class_pos, transm, label='T')    
    plt.plot(class_pos, absorbed, label='1-R-T')    

    plt.xlabel('classical position (fm)')
    plt.ylabel('probability')

    R = refl[-1]
    T = transm[-1]
    A = absorbed[-1]

    ymin = min(A, min(R, T))

    plt.ylim(ymin, 2.0)
    plt.yscale('log')
    plt.title('Probability evolution')
    plt.grid(True)
    plt.legend(loc='upper right')
    fname = output_dir + 'output.png'
    plt.savefig(fname)

    print('Output figure:', fname)

    plt.clf()


def save(logger, x_scale, dx, J, t_scale, dt, N, t_enter, t_exit, t_boundary, output_dir):

    t = logger.time[-1]
    R = logger.refl[-1]
    T = logger.transm[-1]
    A = logger.absorbed[-1]
    psum = R + T + A

    fname = output_dir + 'output.txt'
    f = open(fname, 'w')

    f.write('Simulation terminated at t = {0:.2E} s\n'.format(t))

    f.write('\nResult:\n')
    f.write(' R = {0:.4E}\n T = {1:.4E}\n 1-R-T = {2:.4E}\n norm. = {3:.4f}\n'.format(R, T, A, psum))

    f.write('\nSpatial dim (x):\n')
    f.write(' scale: {0:.2E} fm\n'.format(x_scale))
    f.write(' bin size: {0:.2E} fm\n'.format(dx * x_scale))
    f.write(' no. of bins: {0}\n'.format(J))

    f.write('\nTemporal dim (t):\n')
    f.write(' scale: {0:.2E} s\n'.format(t_scale))
    f.write(' step size: {0:.2E}\n'.format(dt))
    f.write(' no. of steps: {0}\n'.format(N))

    f.write('\nApproximate travel times:\n')
    f.write(' enter barrier:  {0:.2f}\n'.format(t_enter))
    f.write(' exit barrier:   {0:.2f}\n'.format(t_exit))
    f.write(' reach boundary: {0:.2f}\n'.format(t_boundary))

    print('Output file:', fname)

    f.close()

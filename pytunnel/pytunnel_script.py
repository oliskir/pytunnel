import numpy as np
import os
from pint import UnitRegistry
from pytunnel.qm import compute_matrices, solveTDSE, initial_wf, num_rel_err, bin_size
from pytunnel.output import Plotter, Logger, plot_evolution, save, prep_output_folder, peak_position
from pytunnel.parsing import parse_cla, parse_json
import pytunnel.barriers
import pytunnel.units as units


def particle_velocity(energy, mass):
    ''' Returns velocity in fm/s
        energy: kinetic energy in MeV
        mass: rest masss in MeV/c^2
    '''
    return np.sqrt(2 * energy / mass) * units.c


def main():

    args = parse_cla()
    d = parse_json(args.config_file)

    if args.energy is not None: d['energy'] = args.energy

    run(**d)


def run(energy, 
        mass, 
        sigma,
        initial_sep,
        domain_size,
        barrier_type,
        barrier_params,
        absorb_frac=0.01,
        rel_accuracy=0.01,
        display_size=None,
        num_frames=100,
        y_min=0,
        output_dir='output/',
        save_frames=True,
        dx=None,
        dt=None,  
        **kwargs):

    # default display size
    if display_size is None: display_size = domain_size

    # unit conversion
    ureg = UnitRegistry()
    Q = ureg.Quantity

    energy = Q(energy).m_as("MeV")
    mass = Q(mass).m_as("MeV")
    sigma = Q(sigma).m_as("fm")
    position = Q(initial_sep).m_as("fm")
    domain_size = Q(domain_size).m_as("fm")
    display_size = Q(display_size).m_as("fm")
    if dx is not None: dx = Q(dx).m_as("fm")
    if dt is not None: dt = Q(dt).m_as("s")

    # renaming
    termin_thres = absorb_frac
    prop_accuracy = rel_accuracy
    ymin = y_min
    params = barrier_params

    # physical scales
    # spatial dimension = sqrt(hbar * c / (2 * m * c^2))
    E_scale = 1
    x_scale = units.hc / np.sqrt(2. * mass) / np.sqrt(E_scale)
    t_scale = units.hc / units.c / E_scale

#    import pytunnel.units as units
    units.MeV = E_scale
    units.fm = x_scale
    units.sec = t_scale

    if dx is not None: dx = dx / x_scale

    # wave number
    k0sq = 2 * mass * energy / (units.hc**2) - 1. / (2 * sigma**2)
    if k0sq <= 0:
        print('Negative wave number. Not possible to construct a wave-packet with the desired width and kinetic energy. Try increasing the width.')
        exit(1)

    k0 = np.sqrt(k0sq) * x_scale

    # velocity
    velocity = particle_velocity(energy, mass)
    velocity *= t_scale / x_scale

    # 1-sigma width of wave-packet
    sig0 = sigma / x_scale

    # determine bin size
    if prop_accuracy > 0:
        print('\n Determining bin size ...')
        dx = bin_size(sig=sig0, k=k0, velocity=velocity, accuracy=prop_accuracy)
        print(f' dx = {dx*x_scale:.2E} fm')

    # time bin size
    if dt is None: dt = dx / velocity
    else: dt = dt / t_scale

    # start position of wave packet
    xi0 = -position / x_scale

    # computational range
    xJ = 0.5 * domain_size / x_scale
    x0 = xJ - domain_size / x_scale

    # spatial grid size
    J = int((xJ - x0) / dx)

    # barrier
    bstr = "pytunnel.barriers." + barrier_type
    f = eval(bstr)(params, x_start=xi0*x_scale, velocity=velocity*x_scale/t_scale)

    # binned barrier
    barrier = pytunnel.barriers.BinnedBarrier(barrier=f, bins=J, domain_size=domain_size/x_scale)
    j1V = barrier.left_bin(t=0, E=energy)
    j2V = barrier.right_bin(t=0, E=energy) + 1

    # barrier width
    barrier_width = (j2V - j1V) * dx

    # travel times (estimates)
    t_boundary = (xJ - xi0) / velocity
    t_enter = -xi0 / velocity
    t_exit = (barrier_width - xi0) / velocity

    # terminate at boundary
    if 'x_stop' in kwargs.keys():
        x_stop = Q(kwargs['x_stop']).m_as("fm") / x_scale
        t_stop = (x_stop - xi0) / velocity
    else:
        t_stop = t_boundary

    # temporal grid size
    N = int(t_stop / dt)

    # prepare output directory
    output_dir = prep_output_folder(output_dir, save_frames)

    # initial wave packet
    psi0, x = initial_wf(J, dx, x0, sig0, k0, xi0, j1V)

    x_rel_err = num_rel_err(x, psi0, dx, dt, velocity, N=10)
    print('\n Estimating numerical error ...')
    print(f'  rel. err. = {x_rel_err:.2E}')

    # results logger
    velocity_right = particle_velocity(energy+barrier.min(), mass)
    logger = Logger(barrier=barrier, energy=energy, dx=dx, dt=dt, x_axis=x*x_scale, x_start=-position, t_enter=t_enter*t_scale, velocity_left=velocity*x_scale/t_scale, velocity_right=velocity_right, output_dir=output_dir)

    # plotter
    if save_frames:
        xmin = -0.5 * display_size / x_scale 
        xmax = 0.5 * display_size / x_scale 
        jmin = int((xmin - x0) / dx)
        jmax = int((xmax - x0) / dx)
        plotter = Plotter(energy=energy, barrier=barrier, logger=logger, x_axis=x*x_scale, bin1=jmin, bin2=jmax, ymin=ymin, output_dir=output_dir)
    else:
        plotter = None

    # print velocity, de-Broglie wavelength, and barrier enter time
    beta = velocity * x_scale/t_scale / units.c
    lambda_db = 2 * np.pi * units.hc / np.sqrt(2 * mass * energy)
    print(f'\n Velocity = {beta:.3E} c')
    print(f' de-Broglie wavelength = {lambda_db:.3E} fm')

    print('\n Classical times')
    print('---------------------------------')
    print(f'  reach barrier: {t_enter*t_scale:.3E} s')
    print('---------------------------------')

    # print standard deviations
    p_std_rel = units.hc / sigma / np.sqrt(2 * mass * energy)
    print('\n Rel. std. dev.')
    print('---------------------------------')
    print(f'  momentum: {p_std_rel:.3E}')
    print(f'  energy:   {2*p_std_rel:.3E}')
    print('---------------------------------')
    print()

    # solve time-dependent schrödinger equation
    print('Solving time-dependent Schrödinger equation ...')
    logger, psi = solveTDSE(psi0=psi0, N=N, dx=dx, dt=dt, num_frames=num_frames, barrier=barrier, logger=logger, plotter=plotter, termin_thres=termin_thres)

    # write output data to file
    save(logger, x_scale, dx, J, t_scale, dt, N, t_enter, t_exit, t_boundary, output_dir)

    # plot evolutionary diagrams
    plot_evolution(logger, ymin, output_dir)

    if save_frames:
        print('Frames saved to:', output_dir + 'frames/')


if __name__ == '__main__':
   main()





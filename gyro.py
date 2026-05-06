import os
import time
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# This file is a reinterpretation of the IDL file (gyro.pro) for those who want to run it using a Python version.
#
# Python version used: 3.9.6

# Constants
AU = 1.49597870e13                    # Astronomic Unit
arc2cm = (np.pi/180.0)/3600.0 * AU    # Arcsec to cm in Sun
m0 = 9.1094e-28                       # Electron mass [g]
c = 2.998e10                          # [cm/s] Speed of light
e_charge = 4.803e-10                  # Electron charge
E0 = m0 * c**2 / 1.6022e-12 / 1e3     # Electron rest energy [keV]

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def read_gyro(file_name='gyro.out'):
    if not os.path.exists(file_name):
        raise FileNotFoundError(f"Arquivo não encontrado: {file_name}")
        
    with open(file_name, 'r') as f:
        # Skip the first line
        nfq = int(f.readline().strip())
        data = np.loadtxt(f)
        
    # Extracts columns based on IDL behavior.
    fq = data[:, 0]
    jo = data[:, 1]
    jx = data[:, 2]
    ko = data[:, 3]
    kx = data[:, 4]
    return {'fq': fq, 'jo': jo, 'jx': jx, 'ko': ko, 'kx': kx}


def gyro(freq, delta=3.0, energy=None, nel=None, anor=None, ntot=None,
         m=0.0, bmag=500.0, np_dens=1e9, alpha=None, angle=45.0,
         size=12.0, height=5e8, temperature=0.0, quiet=False,
         plot=False, polariz=False, keep=False, **kwargs):
    if energy is None:
        energy = [10.0, 5000.0]
        
    freq = np.atleast_1d(freq).astype(float)
    energy = np.atleast_1d(energy).astype(float)
    m = np.atleast_1d(m)

    vol = height * np.pi * (size / 2.0 * arc2cm)**2

    # Electron density/number of electrons logic (Nel, Ntot, Anor)
    if nel is None and ntot is None and anor is None:
        nel = 1e7

    if ntot is not None:
        ntot = float(ntot)
        nel = ntot / vol
        anor = ntot * (delta - 1) / ((energy[0]/1e3)**(-delta+1) - (energy[1]/1e3)**(-delta+1))
    elif nel is not None:
        ntot = nel * vol
        anor = ntot * (delta - 1) / ((energy[0]/1e3)**(-delta+1) - (energy[1]/1e3)**(-delta+1))
    elif anor is not None:
        if len(energy) > 2:
            raise ValueError('Anor input only for single power-law. Input Ntot or Nel instead.')
        ntot = anor * ((energy[0]/1e3)**(-delta+1) - (energy[1]/1e3)**(-delta+1)) / (delta - 1)
        nel = ntot / vol

    if not quiet:
        print(f'Ntot={ntot}')
        print(f'Anor={anor}')
        print(f'Nel={nel:0.6e}')

    if alpha is not None:
        vb = 0.5 * e_charge / (np.pi * c * m0) * bmag
        vp = 1.5 * vb / alpha
        np_dens = vp**2 * np.pi * m0 / e_charge**2

    ned = len(energy)
    nfreq = len(freq)

    # Anisotropy Logic
    if nel == 0 and nfreq > 0:
        if not quiet:
            print('User set Nel=0; FREE-FREE ONLY')
        r = {'fq': freq, 'jo': np.zeros(nfreq), 'ko': np.zeros(nfreq), 
             'jx': np.zeros(nfreq), 'kx': np.zeros(nfreq)}
    else:
        if len(m) == 2:
            mm = 1
            npd = 2
            pitchangle = m
        elif len(m) > 2:
            if m.ndim == 2 and m.shape[1] == 2:
                mm = 2
                phi_arr = m[:, 0]
                gphi = m[:, 1]
                npd = len(phi_arr)
                pitchangle = np.column_stack((phi_arr, gphi)).flatten()
            else:
                mm = 4
                npd = int(m[0])
                pitchangle = m[1:]
        else:
            mm = 0
            npd = 2
            pitchangle = []

        input_path = os.path.join(BASE_DIR, 'gyro.in')
        output_path = os.path.join(BASE_DIR, 'gyro.out')
        command = os.path.join(BASE_DIR, 'gyro')

        # Write the gyro.in file for C++ be able to read 
        with open(input_path, 'w') as f:
            f.write(f"{bmag} {angle} {np_dens} {nel}\n")
            f.write(f"{ned} {mm} {npd} {nfreq}\n")
            if nfreq != 0:
                f.write(" ".join(map(str, freq)) + "\n")
            f.write(" ".join(map(str, energy)) + f" {delta}\n")
            if mm != 0:
                f.write(" ".join(map(str, pitchangle)) + "\n")

        # Run the C++ executable and measure the execution time.
        time0 = time.time()
        
        try:
            subprocess.run([command], check=True, cwd=BASE_DIR,
                           stdout=subprocess.DEVNULL if quiet else None, 
                           stderr=subprocess.PIPE)
        except FileNotFoundError as e:
            raise FileNotFoundError(
                f"Executable not found in: {command}. Compile the project in this folder."
            ) from e
        except subprocess.CalledProcessError as e:
            print(f"Error executing gyro in C++: {e}")
            
        ctime = time.time() - time0

        if not quiet:
            print('gyro...')
            print(f'Elapsed time: {ctime:.4f} seconds.')

        # Read the gyro.out file
        r = read_gyro(output_path)

        # Remove the files if necessary 
        if not keep:
            for temp_file in [input_path, output_path]:
                if os.path.exists(temp_file):
                    os.remove(temp_file)

    freq_out = r['fq']
    cgs2sfu = 1e19
    area = np.pi * (size / 2.0 * arc2cm)**2
    omega = area / AU**2

    # Free-Free Emission
    if temperature > 0:
        kb = 1.38e-16
        z2 = 1.4
        t1 = np.where(temperature < 2e5, 
                      17.9 + np.log(temperature**1.5) - np.log(freq_out),
                      24.5 + np.log(temperature) - np.log(freq_out))
        kff = 9.78e-3 * np_dens**2 / freq_out**2 / temperature**1.5 * t1 * z2
        jff = kff * kb * temperature * freq_out**2 / c**2
        
        r['jo'] += jff
        r['jx'] += jff
        r['ko'] += kff
        r['kx'] += kff

    # Radiative Transfer
    phi1 = np.zeros(nfreq)
    phi2 = np.zeros(nfreq)
    tau_o = r['ko'] * height
    tau_x = r['kx'] * height

    toLow = tau_o < 0.001
    toMed = (tau_o >= 0.001) & (tau_o <= 20.0)
    toHig = tau_o > 20.0
    
    phi1[toLow] = r['jo'][toLow] * omega * height
    phi1[toMed] = omega * (r['jo'][toMed] / r['ko'][toMed]) * (1.0 - np.exp(-tau_o[toMed]))
    phi1[toHig] = (r['jo'][toHig] / r['ko'][toHig]) * omega

    txLow = tau_x < 0.001
    txMed = (tau_x >= 0.001) & (tau_x <= 20.0)
    txHig = tau_x > 20.0
    
    phi2[txLow] = r['jx'][txLow] * omega * height
    phi2[txMed] = omega * (r['jx'][txMed] / r['kx'][txMed]) * (1.0 - np.exp(-tau_x[txMed]))
    phi2[txHig] = (r['jx'][txHig] / r['kx'][txHig]) * omega
    
    phi1 *= cgs2sfu
    phi2 *= cgs2sfu
    flux = phi1 + phi2

    if plot:
        fig_freq = freq_out / 1e9 # GHz
        style = 'o-' if len(freq_out) < 10 else '-'
        
        if polariz:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
            
            ax1.plot(fig_freq, flux, style)
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.set_ylabel('Flux Density [sfu]')
            
            pol_degree = np.sign(np.cos(np.radians(angle))) * (phi2 - phi1) / np.where(flux > 0, flux, 1e-10)
            ax2.plot(fig_freq, pol_degree, style)
            ax2.set_xscale('log')
            ax2.set_xlabel('Frequency [GHz]')
            ax2.set_ylabel('Polariz. degree')
            ax2.axhline(0, color='black', linestyle='--')
            
            plt.tight_layout()
            plt.show()
        else:
            plt.figure(figsize=(8, 5))
            plt.plot(fig_freq, flux, style)
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('Frequency [GHz]')
            plt.ylabel('Flux Density [sfu]')
            plt.show()

    results = {
        'flux': flux,
        'phi1': phi1,
        'phi2': phi2,
        'e1d': r['jo'],
        'e2d': r['jx'],
        'a1d': r['ko'],
        'a2d': r['kx'],
        'struct': r,
        'electron': {'nel': nel, 'ntot': ntot, 'delta': delta, 'energy': energy},
    }
    
    if 'ctime' in locals():
        results['ctime'] = ctime
        
    return results

if '__main__':
    # How to use
    output = gyro(freq=np.array([1e9, 2e9, 5e9, 1e10]), plot=True, keep=True, polariz=True, quiet=True)
    print(output)

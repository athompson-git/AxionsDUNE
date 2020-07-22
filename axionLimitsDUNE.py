from constants import *

import numpy as np
from numpy import log, log10, sqrt, pi

from axion import PrimakoffAxionFromBeam


# Declare constants
s_per_day = 3600*24
det_mass = 50000
det_am = 37.211e3  # mass of target atom in MeV
det_z = 18  # atomic number
days = 10*365  # days of exposure
det_area = 7*3  # cross-sectional det area
det_thresh = 0.028e-3  # energy threshold
sig_limit = 1.0  # poisson significance (2 sigma)
pot_per_year = 1.1e21




def SandwichSearch(generator, mass_array, g_array, save_file):
    upper_array = np.zeros_like(mass_array)
    lower_array = np.ones_like(mass_array)

    print("Opening ", save_file, "...")
    file_out = open(save_file, "w")

    print("starting scan...")
    for i in range(mass_array.shape[0]):
        generator.axion_mass = mass_array[i]
        print("\n **** SETTING ALP MASS = %0.2e MeV **** \n" % mass_array[i])
        generator.simulate()
        print("Finished flux sim, checking limits")


        # lower bound
        print(" *********** scanning lower bound...")
        for g in g_array:
            generator.axion_coupling = g
            generator.propagate()
            ev_dec = generator.decay_events(days*s_per_day, det_thresh)
            ev_scat = generator.scatter_events(det_mass * mev_per_kg / det_am, det_z, days*s_per_day, det_thresh)
            #print("Events = ", ev_dec, ev_scat, " for g = ", g)
            sig = sqrt(ev_dec+ev_scat)
            if sig > sig_limit:
                print("FOUND DELTA CHI2 = %0.2f at g=%0.2e MeV^-1" % (sig, g))
                lower_array[i] = g
                break

        # upper bound
        print(" ********** scanning upper bound...")
        for g in g_array[::-1]:
            generator.axion_coupling = g
            generator.propagate()
            ev_dec = generator.decay_events(days*s_per_day, det_thresh)
            ev_scat = generator.scatter_events(det_mass * mev_per_kg / det_am, det_z, days*s_per_day, det_thresh)
            #print("Events = ", ev_dec, ev_scat, " for g = ", g)
            sig = sqrt(ev_dec+ev_scat)
            if sig > sig_limit:
                print("FOUND DELTA CHI2 = %0.2f at g=%0.2e MeV^-1" % (sig, g))
                upper_array[i] = g
                break

        out_array = [mass_array[i], lower_array[i], upper_array[i]]
        file_out.write(str(mass_array[i]) + " " + str(lower_array[i]) + " " + str(upper_array[i]) + '\n')

    file_out.close()



def main():
    # Declare event generators.
    flux = np.genfromtxt("data/geant4_flux_DUNE.txt", delimiter=",")
    flux[:,2] *= pot_per_year/10000/s_per_day/365
    flux[:,0] *= 1000
    flux = flux[flux[:,1] < 0.01]  # restrict data to forward photons (0.01 rad ~ 6 degrees)
    print(flux.shape[0])

    # Target
    axion_gen = PrimakoffAxionFromBeam(photon_rates=flux, target_z=6,
                                       target_photon_cross=1e-24, detector_distance=574,
                                       detector_length=10, detector_area=det_area)

    # Run the scan.
    save_dir = "limits/dune_target_limits_20200711.txt"
    mass_array = np.logspace(-4, 2, 5)  # Grid of 5 mass points.
    g_array = np.logspace(-15, -3, 80)  # Grid of 80 coupling points.
    save_file = save_dir                                                                                    

    print("Rerunning limits")
    SandwichSearch(axion_gen, mass_array, g_array, save_file)





if __name__ == "__main__":
    main()


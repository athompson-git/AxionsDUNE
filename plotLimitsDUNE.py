import matplotlib.pyplot as plt
import numpy as np

limits_array = np.genfromtxt("limits/dune_target_limits_20200721_v2.txt")

mass_array = limits_array[:,0]
upper_limit = limits_array[:,2]
lower_limit = limits_array[:,1]

# Find where the upper and lower arrays intersect at the tongue and clip
diff_upper_lower = upper_limit - lower_limit
upper_limit = np.delete(upper_limit, np.where(diff_upper_lower < 0))
lower_limit = np.delete(lower_limit, np.where(diff_upper_lower < 0))
mass_array = np.delete(mass_array, np.where(diff_upper_lower < 0))

# join upper and lower bounds
joined_limits = np.append(lower_limit, upper_limit[::-1])
joined_masses = np.append(mass_array, mass_array[::-1])


# Read in data.
beam = np.genfromtxt('data/constraints/beam.txt')
eeinva = np.genfromtxt('data/constraints/eeinva.txt')
lep = np.genfromtxt('data/constraints/lep.txt')
lsw = np.genfromtxt('data/constraints/lsw.txt')
nomad = np.genfromtxt('data/constraints/nomad.txt')

# Astrophyiscal limits
cast = np.genfromtxt("data/constraints/cast.txt", delimiter=",")
hbstars = np.genfromtxt("data/constraints/hbstars.txt", delimiter=",")
sn1987a_upper = np.genfromtxt("data/constraints/sn1987a_upper.txt", delimiter=",")
sn1987a_lower = np.genfromtxt("data/constraints/sn1987a_lower.txt", delimiter=",")

plt.plot(joined_masses*1e6, joined_limits*1e3, color="crimson", label='DUNE ND')

# Plot astrophysical limits
plt.fill(hbstars[:,0]*1e9, hbstars[:,1]*0.367e-3, label="HB Stars", color="mediumpurple", alpha=0.3)
plt.fill(cast[:,0]*1e9, cast[:,1]*0.367e-3, label="CAST", color="orchid", alpha=0.3)
plt.fill_between(sn1987a_lower[:,0]*1e9, y1=sn1987a_lower[:,1]*0.367e-3, y2=sn1987a_upper[:,1]*0.367e-3,
                label="SN1987a", color="lightsteelblue", alpha=0.3)


# Plot lab limits
plt.fill(beam[:,0], beam[:,1], label='Beam Dump', color="b", alpha=0.7)
plt.fill(np.hstack((eeinva[:,0], np.min(eeinva[:,0]))), np.hstack((eeinva[:,1], np.max(eeinva[:,1]))),
        color="orange", label=r'$e^+e^-\rightarrow inv.+\gamma$', alpha=0.7)
plt.fill(lep[:,0], lep[:,1], label='LEP', color="green", alpha=0.7)
plt.fill(np.hstack((nomad[:,0], np.min(nomad[:,0]))), np.hstack((nomad[:,1], np.max(nomad[:,1]))),
        color="yellow", label='NOMAD', alpha=0.7)


plt.legend(loc="lower left", framealpha=1, ncol=2, fontsize=9)
plt.title(r"Primakoff, $a\to\gamma\gamma$, 50t, 10 years exposure", loc="right")
plt.xscale('log')
plt.yscale('log')
plt.xlim((1,1e10))
plt.ylim(1e-13,1e-1)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel('$m_a$ [eV]', fontsize=15)
plt.ylabel('$g_{a\gamma\gamma}$ [GeV$^{-1}$]', fontsize=15)

plt.tick_params(axis='x', which='minor')
plt.tight_layout()
plt.show()
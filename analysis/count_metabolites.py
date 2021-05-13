#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3

from sqlite3 import connect
from matplotlib import pyplot as plt
from matplotlib import rcParams
from numpy import array, count_nonzero, logical_and

rcParams['font.size'] = 8


def neutral_mass(mz, adduct):
    adduct_to_mass = {'[M+H]+': 1.0078, '[M+K]+': 38.9637, '[M+H-H2O]+': -17.0027, '[M+Na]+': 22.9898, '[M]+': 0.}
    return mz - adduct_to_mass[adduct]


def mass_shift_hist(mass_shifts):
    fig = plt.figure(figsize=(7, 2))
    ax = fig.add_subplot(111)

    ax.axvline(0, ls='--', c='k', lw=1, zorder=-1)
    ax.hist(mass_shifts, bins=[_ for _ in range(-30, 325)], color='b', histtype='stepfilled')

    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xlabel('mass shift')
    ax.set_ylabel('N')

    plt.savefig('mass_shifts.png', bbox_inches='tight', dpi=350)
    #plt.show()
    plt.close()


def main():

    con = connect('DMIM_v1.0.db')
    cur1, cur2 = con.cursor(), con.cursor()

    qry_parent = 'SELECT well, mz, adduct FROM plate_{n} JOIN plate_{n}_id ON plate_{n}.dmim_id = plate_{n}_id.dmim_id WHERE met_n = 0'
    qry_metab = 'SELECT mz, adduct FROM plate_{n} JOIN plate_{n}_id ON plate_{n}.dmim_id = plate_{n}_id.dmim_id WHERE well = ? AND met_n > 0'

    mass_shifts = []
    i = 0
    prev_well = ''
    for n in range(1, 8):
        for well, parent_mz, parent_adduct in cur1.execute(qry_parent.format(n=n)):
            parent_neutral_mass = neutral_mass(float(parent_mz), parent_adduct)
            #print(well, parent_mz, parent_adduct)
            if well != prev_well:
                for metab_mz, metab_adduct in cur2.execute(qry_metab.format(n=n), (well,)):
                    mass_shifts.append(neutral_mass(float(metab_mz), metab_adduct) - parent_neutral_mass)
            prev_well = well

    mass_shifts = array(mass_shifts)

    mass_shift_hist(mass_shifts)

    # count some specific metabolites
    print('+GSH +O', count_nonzero(logical_and(mass_shifts >= 322.8, mass_shifts <= 323.8)))
    print('+GSH', count_nonzero(logical_and(mass_shifts >= 306.8, mass_shifts <= 307.8)))
    print('+Glc +O', count_nonzero(logical_and(mass_shifts >= 191.8, mass_shifts <= 192.8)))
    print('+Glc', count_nonzero(logical_and(mass_shifts >= 175.8, mass_shifts <= 176.8)))
    print('+Glc -Me', count_nonzero(logical_and(mass_shifts >= 161.8, mass_shifts <= 162.8)))
    print('+2O', count_nonzero(logical_and(mass_shifts >= 31.8, mass_shifts <= 32.8)))
    print('+O', count_nonzero(logical_and(mass_shifts >= 15.8, mass_shifts <= 16.8)))
    print('-2H', count_nonzero(logical_and(mass_shifts >= -2.5, mass_shifts <= -1.5)))
    print('-Me', count_nonzero(logical_and(mass_shifts >= -14.5, mass_shifts <= -13.5)))
    print('-2Me/-Et', count_nonzero(logical_and(mass_shifts >= -28.5, mass_shifts <= -27.5)))



if __name__ == '__main__':
    main()

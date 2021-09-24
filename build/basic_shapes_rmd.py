#!/usr/local/Cellar/python@3.9/3.9.1_6/bin/python3 


from matplotlib import pyplot as plt
from random import random
from math import pi, sin, cos

from add_md3d import load_and_center_xyzmq, mass_dist


def get_mass_dist(str_data):
    mcxyz, m = load_and_center_xyzmq(str_data)
    return mass_dist(mcxyz)


def plot_mass_dist(mdist, title):
    fig = plt.figure(figsize=(3, 2))
    ax = fig.add_subplot(111)
    ax.hist(mdist, bins=100, range=(0, 1000), color='b')
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    ax.set_xlabel('r →')
    plt.savefig('basic_{}_rmd.png'.format(title), dpi=300, bbox_inches='tight')
    plt.close()


def rod():
    # generate data
    s = ''
    for i in range(1000):
        s += '{:.3f} {:.3f} {:.3f} {:.3f} {:.3f}\n'.format(float(i), 0., 0., 1., 0.)
    plot_mass_dist(get_mass_dist(s), 'rod')


def ring():
    # generate data
    s = ''
    for i in range(1000):
        theta = random() * 2. * pi
        x = sin(theta) * 500.
        y = cos(theta) * 500.
        s += '{:.3f} {:.3f} {:.3f} {:.3f} {:.3f}\n'.format(x, y, 0., 1., 0.)
    plot_mass_dist(get_mass_dist(s), 'ring')


def disc():
    # generate data
    s = ''
    for r in range(1, 501, 10):
        for i in range(r):
            theta = random() * 2. * pi
            x = sin(theta) * float(r)
            y = cos(theta) * float(r)
            s += '{:.3f} {:.3f} {:.3f} {:.3f} {:.3f}\n'.format(x, y, 0., 1., 0.)
    plot_mass_dist(get_mass_dist(s), 'disc')


def hollow_sphere():
    # generate data
    s = ''
    for i in range(1000):
        theta = random() * pi
        phi = random() * 2. * pi
        x = sin(theta) * cos(phi) * 500.
        y = sin(theta) * sin(phi) * 500.
        z = cos(theta) * 500.
        s += '{:.3f} {:.3f} {:.3f} {:.3f} {:.3f}\n'.format(x, y, z, 1., 0.)
    plot_mass_dist(get_mass_dist(s), 'hollow sphere')


def sphere():
    # generate data
    s = ''
    for r in range(1, 501, 10):
        for i in range(int(r**2)):
            theta = random() * pi
            phi = random() * 2. * pi
            x = sin(theta) * cos(phi) * float(r)
            y = sin(theta) * sin(phi) * float(r)
            z = cos(theta) * float(r)
            s += '{:.3f} {:.3f} {:.3f} {:.3f} {:.3f}\n'.format(x, y, z, 1., 0.)
    plot_mass_dist(get_mass_dist(s), 'sphere')



def main():
    rod()
    ring()
    disc()
    hollow_sphere()
    sphere()



if __name__ == '__main__':
    main()




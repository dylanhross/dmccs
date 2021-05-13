#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""
    Runs all of the individual analysis scripts, regenerating all figures
"""

from compare_with_ccsbase import main as cwc_main
from pca import main as pca_main
from plsra import main as pls_main


def main():
    """ main execution sequence """

    # run comparisons of DMIM against CCSbase
    cwc_main()
    # PCA on DMIM data
    pca_main()
    # PLS-RA with CCS on DMIM data
    pls_main()

if __name__ == '__main__':
    main()

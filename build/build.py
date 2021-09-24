#!/usr/local/Cellar/python@3.9/3.9.1_6/bin/python3 
"""
    Builds DMIM_v?.?.db from individual components
"""

from initialize_db import main as init_main
from add_measurement_data import main as meas_main
from add_annotations import main as annt_main
from add_ms2 import main as ms2_main
from add_mqns import main as mqn_main
from add_3d import main as a3d_main
from add_md3d import main as md3d_main
from filter_metabs_metfrag_score import main as filt_main
from final_table_counts import main as cnt_main


def main():
    """ main build sequence """
    
    VERSION = '1.1'

    # initialize database
    init_main(VERSION)
    # add in measurement data
    meas_main(VERSION)
    # add annotations
    annt_main(VERSION)
    # add in MS2 spectra
    ms2_main(VERSION)
    """
    # add in MQNs for all annotated species
    mqn_main(VERSION)
    # add in the 3D structures
    a3d_main(VERSION)
    # add MD3Ds for 3D structures
    md3d_main(VERSION)
    # filter out metabolites using MetFrag score
    filt_main(VERSION)
    # report the final table counts
    cnt_main(VERSION)"""


if __name__ == '__main__':
    main()


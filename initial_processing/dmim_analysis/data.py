"""
    dmim_analysis/data.py
    Dylan H. Ross
    2019/06/25

    description:
        A set of objects defining a data storage hierarchy for data from the high-throughput DMIM project 
"""


class Adduct:
    """
Adduct
    description:
        A class for storing data about an observed MS adduct, this includes the ATD and LC chromatogram, and a mass 
        spectrum filtered on retention time and drift time.
"""

    def __init__(self, adduct, monoiso, measured_mz, tolerance, atd, atd_fit_params, chromatogram, rt_fit_params, dt,
                 rt, ccs, d_file, ctrl_d_file=None, no_cofac_dti=None):
        """
Adduct.__init__
    description:
        Initializes a new Adduct
    parameters:
        adduct (str) -- type of MS adduct
        monoiso (float) -- monoisotopic mass of the adduct
        measured_mz (float) -- actual measured m/z value (may differ slightly from expected monoisotopic mass...)
        tolerance (float) -- tolerance used for mass filtering
        atd (numpy.ndarray) -- a 2D array containing the arrival time distribution [[dt], [intensity]]
        atd_fit_params (tuple(float)) -- parameters (A, td, C) from gaussian fit on ATD
        chromatogram (numpy.ndarray) -- a 2D array containing the LC chromatogram [[rt], [intensity]]
        rt_fit_params (tuple(float)) -- parameters (A, rt, C) from gaussian fit on chromatogram
        dt (float) -- extracted drift time (ms)
        rt (float) -- fitted retention time (min)
        ccs (float) -- calibrated CCS
        d_file (str) -- name of the data file this data was taken from
        [ctrl_d_file (str)] --
        [no_cofac_dti (numpy.ndarray(float))] -- intensities from extracted ATD from non-cofactor reaction, for 
                                                    metabolites only [optional, default=None]
"""
        self.adduct = adduct
        self.monoiso = monoiso
        self.measured_mz = measured_mz
        self.tolerance = tolerance
        self.atd = atd
        self.atd_fit_params = atd_fit_params
        self.chromatogram = chromatogram
        self.rt_fit_params = rt_fit_params
        self.dt = dt
        self.rt = rt
        self.ccs = ccs
        self.d_file = d_file
        self.ctrl_d_file = ctrl_d_file
        self.no_cofac_dti = no_cofac_dti

    def __repr__(self):
        """
Adduct.__repr__
    description:
        Produces a string representation of this Adduct instance
    returns:
        (str) -- string representation of this Adduct instance
"""
        base = 'Adduct(adduct="{}", monoiso={:.4f}, measured_mz={:.4f}, tolerance={:.4f}'\
               ' dt={:.2f}, rt={:.2f}, ccs={:.1f}, d_file="{}")'
        return base.format(self.adduct, self.monoiso, self.measured_mz, self.tolerance, self.dt, self.rt, self.ccs,
                           self.d_file)


class Metabolite:
    """
Metabolite
    description:
        A class for storing data on predicted metabolites, most importantly the reaction, monoisotopic mass, and 
        InChI key defining the structure, along with any observed MS adducts (a list of Adduct objects)
"""

    def __init__(self, id_, reaction, monoiso, InChI, png=None, adducts=[]):
        """
Metabolite.__init__
    description:
        Initializes a new Metabolite
    parameters:
        id_ (int) -- integer identifier 
        reaction (str) -- reaction that formed the metabolite
        monoiso (float) -- monoisotopic mass of the metabolite (neutral)
        InChI (str) -- InChI key defining metabolite structure
        [png (bytes)] -- image of metabolite structure [optional, default=None] 
        [adducts (list(Adduct))] -- list of observed MS adducts [optional, default=[]] 
"""
        self.id = id_ 
        self.reaction = reaction
        self.monoiso = monoiso
        self.InChI = InChI
        self.png = png
        self.adducts = adducts

    def __repr__(self):
        """
Metabolite.__repr__
    description:
        Produces a string representation of this Metabolite instance
    returns:
        (str) -- string representation of this Metabolite instance
"""
        base = 'Metabolite(id={:03d}, reaction="{}", monoiso={:.4f}, adducts={})'
        return base.format(self.id, self.reaction, self.monoiso, self.adducts)


class Compound:
    """
Compound
    description:
        TODO
"""

    def __init__(self, id_, name, monoiso, SMILES, png=None, adducts=[], metabolites=[]):
        """
Compound.__init__
    description:
        Initializes a new Compound
    parameters:
        id_ (int) -- integer identifier 
        name (str) -- compound name
        monoiso (float) -- monoisotopic mass of the metabolite (neutral)
        SMILES (str) -- compound SMILES structure
        [png (bytes)] -- image of compound structure [optional, default=None] 
        [adducts (list(Adduct))] -- list of observed MS adducts [optional, default=[]] 
        [metabolites (list(Metabolite))] -- list of predicted metabolites [optional, default=[]]
"""
        self.id = id_ 
        self.name = name
        self.monoiso = monoiso
        self.SMILES = SMILES
        self.png = png
        self.adducts = adducts
        self.metabolites = metabolites

    def __repr__(self):
        """
Compound.__repr__
    description:
        Produces a string representation of this Compound instance
    returns:
        (str) -- string representation of this Compound instance
"""
        base = 'Compound(id={:03d}, name="{}", monoiso={:.4f}, adducts={}, metabolites={})'
        return base.format(self.id, self.name, self.monoiso, self.adducts, self.metabolites)


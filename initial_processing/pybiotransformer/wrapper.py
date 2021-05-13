"""
    pybiotransformer/wrapper.py
    Dylan H. Ross
    2018/11/19

    description:
        Python wrapper (class) for using BioTransformer to identify drug metabolites.
"""


from subprocess import run, DEVNULL
from tempfile import NamedTemporaryFile
from csv import DictReader
import os


class PBTMetabolite:
    """
PBTMetabolite
    description:
        data structure for representing predicted metabolites
"""

    # define MS adducts and their mass shifts (relative to neutral mass)
    ms_adducts = {
        # positive mode
        '[M+H]+': 1.0078,
        '[M+Na]+': 22.9898,
        '[M+H-H2O]+': -17.0027,
        # negative mode
        '[M-H]-': -1.0078
    }
    
    def __init__(self, mass, reaction, InChI):
        """
PBTMetabolite.__init__
    description:
        Initialize a new PBTMetabolite object with a mass, reaction type and InChI structure
    parameters:
        mass (float) -- exact (neutral) mass of the predicted metabolite
        reaction (str) -- description of the reaction that produces the metabolite
        InChI (str) -- InChI structure of the metabolite
"""
        # only keep the mass to 4 decimal places to make comparisons more convenient
        self.mass = round(mass, 4)
        self.reaction = reaction
        self.InChI = InChI
        self.smiles = None
        self.png = None

    def __repr__(self):
        """
PBTMetabolite.__repr__
    description:
        String representation of this PBTMetabolite
"""
        s = "PBTMetabolite(mass={:.4f}, reaction='{}', InChI='{}')"
        return s.format(self.mass, self.reaction, self.InChI)

    def __eq__(self, other):
        """
PBTMetabolite.__eq__
    description:
        Equality comparison for PBTMetabolite objects, compare based on the neutral metabolite mass
    parameters:
        other (PBTMetabolite) -- other object for equality comparison
    returns:
        (bool) -- object equality
"""
        return self.mass == other.mass

    def __hash__(self):
        """
PBTMetabolite.__hash__
    description:
        Make PBTMetabolites hashable so that they can be used in sets, simply hash on the mass
    returns:
        (int) -- hash of this PBTMetabolite
"""
        return hash(self.mass)

    def ms_adduct(self, adduct):
        """
PBTMetabolite.ms_adduct
    description:
        Computes the monoisotopic mass of a specified MS adduct of this metabolite.
    parameters:
        adduct (str) -- specify the MS adduct
    returns:
        (float) -- monoisotopic mass of the MS adduct
"""
        if adduct not in self.ms_adducts:
            msg = "PBTMetabolite: ms_adduct: MS adduct '{}' not recognized"
            raise ValueError(msg.format(adduct))
        return round(self.mass + self.ms_adducts[adduct], 4)


class PBTWrapper:
    """
PBTWrapper
    description:
        Wrapper for using BioTransformer to identify drug metabolites.
"""
    
    # recognized BioTransformer types
    # only include subset that is compatible with the -m flag, per documentation
    bt_types = [
        "allHuman",
        "superbio",
        "env"
    ]

    def __init__(self, bt_jar_path=None):
        """
PBTWrapper.__init__
    description:
        Initialize a PBTWrapper. Requires path to the BioTransformer jar.
    parameters:
        bt_jar_path (str) -- path to BioTransformer jar
"""
        # store the jar path
        if bt_jar_path:
            self.bt_jar_path_ = bt_jar_path
        else:
            self.bt_jar_path_ = os.path.normpath(os.path.join(os.path.split(__file__)[0], 
                                                              '../', 
                                                              'biotransformer-1-0-7.jar'))

    def predict_metabolites(self, parent_smi, bt_type="allHuman", n_steps=1, output_csv=None, unique=False):
        """
PBTWrapper.predict_metabolites
    description:
        Predicts all possible metabolites of a parent compound up to a specifcied number of metabolic steps (n_steps)
        using a specified BioTransformer. Output is produced in .csv format and written to a NamedTemporaryFile, 
        which is then parsed for the relevant information. Alternatively, a .csv filename may be provided (output_csv)
        and the output will just be directly written to that.
    parameters:
        parent_smi (str) -- smiles structure of parent compound
        [bt_type (str)] -- BioTransformer type, must be defined in self.bt_types [optional, default="allHuman"]
        [n_steps (int)] -- maximum number of metabolic steps to search for. From the documentation: 'This option 
                            can be set by the user for the EC-based, CYP450, Phase II, and Environmental microbial
                            biotransformers' [optional, default=1]
        [output_csv (csv)] -- output may be directly written to the specified .csv file instead of being dumped into
                              a temporary file and parsed (if output_csv is None) [optional, default=None]
        [unique (bool)] -- return a set of PBTMetabolites instead of a list, so that only unique masses are included
                           [optional, default=False]
    returns:
        (list(PBTMetabolite), set(PBTMetabolite), or None) -- returns a list (or set) of PBTMetabolites parsed from 
                                                              the output file or returns None if output is directly 
                                                              written to a .csv file or if there are any errors
"""
        # make sure the BioTransformer is acceptable
        if bt_type not in self.bt_types:
            msg = "PBTWrapper: predict_metabolites: bt_type '{}' not allowed"
            raise ValueError(msg.format(bt_type))

        if output_csv:
            out_name = output_csv
        else:
            temp_out = NamedTemporaryFile(delete=False)
            out_name = temp_out.name

        # build list of args to run the command
        call_args = [
            'java', '-jar', self.bt_jar_path_, 
            '-k', 'pred',
            '-b', bt_type,
            '-ismi', parent_smi,
            '-ocsv', out_name,
            '-s', str(n_steps)
        ]

        # run the command in the same directory as biotransformer.jar
        curr_path = os.getcwd()
        os.chdir(os.path.split(self.bt_jar_path_)[0])

        try:
            run(call_args, stdout=DEVNULL, stderr=DEVNULL)
        except Exception as e:
            msg = "PBTWrapper: predict_metabolites: {} call failed"
            print(e)
            print(msg.format(self.bt_jar_path_))

        # go back to whatever directory we were in before
        os.chdir(curr_path)

        out = []
        if not output_csv:
            try:
                with open(out_name, "r") as out_f:
                    out = [PBTMetabolite(float(row['Major Isotope Mass']), row['Reaction'], row['InChI']) 
                           for row in DictReader(out_f)]
                    
                    # only unique metabolites
                    if out and unique:
                        s = set()
                        for m in out:
                            s.add(m)
                        out = s
            except Exception as e:
                msg = "PBTWrapper: predict_metabolites: failed to parse output ({})"
                print(e)
                print(msg.format(out_name))
        return out


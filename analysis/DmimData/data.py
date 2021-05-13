"""
    DmimData/data.py
    Dylan H. Ross
    2021/01/15
    
    description:
        TODO
"""


from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
from numpy import concatenate, percentile, digitize, count_nonzero, array

from DmimData.db_interface import qry_mz, qry_mqn, qry_md3d, qry_combined


class DMD:
    """
DMD
    description:
        Object for interfacing with the DMIM_v?.?.db database and retrieving data from it. Responsible for filtering data
        on various criteria, producing features for ML, handling test/train set splitting, 
        and data transfomations (normalization/centering/scaling)
"""

    def __init__(self, db_path, seed):
        """
DMD.__init__
    description:
        Initializes a new DMD object using the path to the DMIM_v?.?.db database file. 
        The database path and pRNG seed are stored in instance variables, respectively:
            self.db_path_
            self.seed_
        The following instance variables are initialized as None (must be set by calls to other methods):
            self.name_          (name or annotation -> set by self.featurize(...))
            self.adduct_        (MS adduct -> set by self.featurize(...))
            self.met_n_         (metabolite number -> set by self.featurize(...))
            self.mz_            (m/z -> set by self.featurize(...))
            self.n_             (number of entries -> set by self.featurize(...))
            self.n_parent_      (number of parent drugs -> set by self.featurize(...))
            self.n_metab_       (number of metabolites -> set by self.featurize(...))
            self.X_             (full array of features -> set by self.featurize(...))
            self.y_             (full array of labels -> set by self.featurize(...))
            self.n_features_    (number of features -> set by self.featurize(...))
            self.X_train_       (training set split of features -> set by self.train_test_split(...))
            self.y_train_       (training set split of labels -> set by self.train_test_split(...))
            self.n_train_       (training set size -> set by self.train_test_split(...)) 
            self.X_test_        (test set split of features -> set by self.train_test_split(...))
            self.y_test_        (test set split of labels -> set by self.train_test_split(...))
            self.n_test_        (test set size -> set by self.train_test_split(...))
            self.SSSplit_       (StratifiedShuffleSplit instance -> set by self.train_test_split(...))
            self.SScaler_       (StandardScaler instance -> set by self.center_and_scale(...))
            self.X_train_ss_    (centered/scaled training set features -> set by self.center_and_scale(...))
            self.X_test_ss_     (centered/scaled test set features -> set by self.center_and_scale(...))
    parameters:
        db_path (str) -- path to C3S.db database file  
        seed (int) -- pRNG seed to use for any data preparation steps with a stochastic component, stored in the
                        self.seed_ instance variable [optional, default=69]
"""
        # store database file path and pRNG seed
        self.db_path_ = db_path
        self.seed_ = seed
        # declare instance variables to use later
        self.X_ = None
        self.y_ = None
        self.name_ = None
        self.adduct_ = None
        self.met_n_ = None
        self.mz_ = None
        self.n_ = None
        self.n_features_ = None
        self.n_parent_ = None
        self.n_metab_ = None
        self.X_train_ = None
        self.y_train_ = None
        self.n_train_ = None
        self.X_test_ = None
        self.y_test_ = None
        self.n_test_ = None
        self.SSSplit_ = None
        self.SScaler_ = None
        self.X_train_ss_ = None
        self.X_test_ss_ = None


    def mqns_to_indices(self, mqns):
        """
DMD.mqns_to_indices
    description:
        Converts a list of named MQNs into MQN indices
    parameters:
        mqns (list(str)) -- list of named MQNs
    returns:
        mqn_indices (list(int)) -- list of indices for desired MQNs
"""
        m2i = {
            'c': 0, 'f': 1, 'cl': 2, 'br': 3, 
            'i': 4, 's': 5, 'p': 6, 'an': 7,
            'cn': 8, 'ao': 9, 'co': 10, 'hac': 11,
            'hbam': 12, 'hba': 13, 'hbdm': 14,
            'hbd': 15, 'neg': 16, 'pos': 17,
            'asb': 18, 'adb': 19, 'atb': 20, 'csb': 21,
            'cdb': 22, 'ctb': 23, 'rbc': 24,
            'asv': 25, 'adv': 26, 'atv': 27, 'aqv': 28,
            'cdv': 29, 'ctv': 30, 'cqv': 31, 'r3': 32,
            'r4': 33, 'r5': 34, 'r6': 35, 'r7': 36,
            'r8': 37, 'r9': 38, 'rg10': 39, 
            'afr': 40, 'bfr': 41
        }
        return [m2i[_] for _ in mqns]


    def md3ds_to_indices(self, md3ds):
        """
DMD.md3ds_to_indices
    description:
        Converts a list of named MD3Ds into MD3D indices
    parameters:
        md3ds (list(str)) -- list of named MD3Ds
    returns:
        md3d_indices (list(int)) -- list of indices for desired MD3Ds
"""
        m2i = {
            'pmi1': 0, 'pmi2': 1, 'pmi3': 2,  
            'rmd02': 3, 'rmd24': 4, 'rmd46': 5, 'rmd68': 6, 'rmd8p': 7
        }
        return [m2i[_] for _ in md3ds]


    def featurize(self, features, custom_mqns=None, custom_md3ds=None):
        """
DMD.featurize
    description:
        Generates features for the dataset,
        Sets the following instance variables:
            self.name_          (name or annotation)
            self.adduct_        (MS adduct)
            self.met_n_         (metabolite number)
            self.mz_            (m/z)
            self.n_             (number of entries)
            self.X_             (full array of features)
            self.y_             (full array of labels)
            self.n_features_    (number of features)
            self.n_parent_      (number of parent drugs)
            self.n_metab_       (number of metabolites)
    parameters:
        features (str) -- specify the feature set to use (must be one of 'mz', 'mqn', 'md3d', 'combined', 'custom')
        custom_mqns (None or list(str)) -- used with 'custom' feature set to determine specific MQNs to include
        custom_md3ds (None or list(str)) -- used with 'custom' feature set to determine specific MD3Ds to include
"""
        # make sure the feature set 
        if features not in ['mz', 'mqn', 'md3d', 'combined', 'custom']:
            m = 'DMD: featurize: feature set "{features}" invalid'.format(features=features)
            raise ValueError(m)
        qry = {
            'mz': qry_mz,
            'mqn': qry_mqn,
            'md3d': qry_md3d,
            'combined': qry_combined, 
            'custom': qry_combined
        }
        X, y, name, adduct, met_n, mz = qry[features](self.db_path_)
        # further processing required for 'mz' and 'custom' feature sets
        if features == 'mz':
            X = X.reshape(-1, 1)
        elif features == 'custom':
            # make sure both of the sets of indices have been provided
            if custom_mqns is None or custom_md3ds is None:
                m = 'DMD: featurize: both custom_mqns and custom_md3ds must be set when using "custom" feature set'
                raise ValueError(m)
            mqn_indices = self.mqns_to_indices(custom_mqns)
            md3d_indices = self.md3ds_to_indices(custom_md3ds)
            # shift md3d indices by 42 since they are added at the end 
            md3d_indices = [i + 42 for i in md3d_indices]
            all_indices = mqn_indices + md3d_indices
            # keep only the desired indices
            X = array([_[all_indices] for _ in X])
        self.name_ = name
        self.adduct_ = adduct
        self.met_n_ = met_n
        self.mz_ = mz
        self.X_ = X
        self.y_ = y
        self.n_, self.n_features_ = self.X_.shape
        self.n_parent_ = len([_ for _ in self.met_n_ if _ == 0])
        self.n_metab_ = self.n_ - self.n_parent_


    def train_test_split(self, test_frac=0.2):
        """
DMD.train_test_split
    description:
        Shuffles the data then splits it into a training set and a test set, storing each in self.X_train_,
        self.y_train_, self.X_test_, self.y_test_ instance variables. The splitting is done in a stratified manner
        based on either CCS or dataset source. In the former case, the CCS distribution in the complete dataset is
        binned into a rough histogram (8 bins) and the train/test sets are split such that they each contain similar 
        proportions of this roughly binned CCS distribution. In the latter case, the train/test sets are split such
        that they each preserve the rough proportions of all dataset sources present in the complete dataset.
        This method DOES NOT get called on DMD objects in the self.datasets_ instance variable.

        Sets the following instance variables:
            self.X_train_       (training set split of features)
            self.y_train_       (training set split of labels)
            self.n_train_       (training set size) 
            self.X_test_        (test set split of features)
            self.y_test_        (test set split of labels)
            self.n_test_        (test set size)
            self.SSSplit_       (StratifiedShuffleSplit instance)

        ! self.featurize(...) must be called first to generate the features and labels (self.X_, self.y_) !
    parameters:
        [test_frac (float)] -- fraction of the complete dataset to reserve as a test set, defaults to an 80 % / 20 %
                               split for the train / test sets, respectively [optional, default=0.2]
"""
        # make sure self.featurize(...) has been called
        if self.X_ is None:
            msg = 'DMD: train_test_split: self.X_ is not initialized, self.featurize(...) must be called before ' + \
                    'calling self.train_test_split(...)'
            raise RuntimeError(msg)
        y_cat = self.get_categorical_y()
        # initialize StratifiedShuffleSplit
        self.SSSplit_ = StratifiedShuffleSplit(n_splits=1, test_size=test_frac, random_state=self.seed_)
        # split and store the X and y train/test sets as instance variables
        for train_index, test_index in self.SSSplit_.split(self.X_, y_cat):
            self.X_train_, self.X_test_ = self.X_[train_index], self.X_[test_index]
            self.y_train_, self.y_test_ = self.y_[train_index], self.y_[test_index]
        # store the size of the train/test sets in instance variables
        self.n_train_ = self.X_train_.shape[0] 
        self.n_test_ = self.X_test_.shape[0] 


    def get_categorical_y(self):
        """
DMD.get_categorical_y
    description:
        transforms the labels into 'categorical' data (required by StratifiedShuffleSplit) by performing a binning
        operation on the continuous label data. The binning is performed using the rough distribution of label values
        with the following bounds (based on quartiles):
                Q1            Q2             Q3
                |             |              |
          bin1  | bin2 | bin3 | bin4 | bin 5 | bin 6
                       |             |
              (Q2 - Q1 / 2) + Q1     |
                            (Q3 - Q2 / 2) + Q2
        Uses labels stored in the self.y_ instance variable
    returns:
        (np.ndarray(int)) -- categorical (binned) label data
"""
        # get the quartiles from the label distribution
        q1, q2, q3 = percentile(self.y_, [25, 50, 75])
        # get midpoints
        mp12 = q1 + (q2 - q1) / 2.
        mp23 = q2 + (q3 - q2) / 2.
        # bin boundaries
        bounds = [q1, mp12, q2, mp23, q3]
        return digitize(self.y_, bounds)


    def center_and_scale(self):
        """
DMD.center_and_scale
    description:
        Centers and scales the training set features such that each has an average of 0 and variance of 1. Applies
        this transformation to the training and testing features, storing the results in the self.X_train_ss_ and 
        self.X_test_ss_ instance variables, respectively. Also stores a reference to the fitted StandardScaler
        instance for use with future data.
        This method DOES NOT get called on DMD objects in the self.datasets_ instance variable.

        sets the following instance variables:
            self.SScaler_       (StandardScaler instance)
            self.X_train_ss_    (centered/scaled training set features)
            self.X_test_ss_     (centered/scaled test set features)

        ! self.train_test_split(...) must be called first to generate the training features and labels (self.X_train_, 
        self.y_train_) that are used to initialize the StandardScaler !
"""
        if self.X_train_ is None:
            msg = 'DMD: center_and_scale: self.X_train_ is not initialized, self.train_test_split(...) must be ' + \
                  'called before calling self.center_and_scale(...)'
            raise RuntimeError(msg)

        # perform the scaling
        self.SScaler_ = StandardScaler()
        self.X_train_ss_ = self.SScaler_.fit_transform(self.X_train_)
        self.X_test_ss_ = self.SScaler_.transform(self.X_test_)








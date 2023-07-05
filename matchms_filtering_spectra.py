import pickle
import numpy as np

from time import time
from tqdm import tqdm
from functools import wraps
from joblib import Parallel, delayed
from matchms.filtering import add_losses
from matchms.filtering import add_parent_mass
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import reduce_to_number_of_peaks
from matchms.filtering import require_minimum_number_of_peaks
from matchms.filtering import select_by_mz
from matchms.importing import load_from_mgf

# Decorator to time functions
def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        print(f'func:{f.__name__} args:[{args}, {kw}] took: {te-ts:2.4f} sec')
        return result
    return wrap

# Apply filters to clean and enhance each spectrum
def peak_processing(s):
    s = default_filters(s)
    s = add_parent_mass(s)
    s = normalize_intensities(s)
    s = reduce_to_number_of_peaks(s, n_required=10, ratio_desired=0.5, n_max=500)
    #s = select_by_intensity(spectrum, intensity_from=0.01)
    s = select_by_mz(s, mz_from=20, mz_to=2000)
    s = add_losses(s, loss_mz_from=10.0, loss_mz_to=200.0)
    s = require_minimum_number_of_peaks(s, n_required=10)
    return s

@timing
def main():
    processed_spectrums = Parallel(n_jobs=NPROC)(delayed(peak_processing)(spectrum) for spectrum in tqdm(spectrums))
    processed_spectrums = [s for s in processed_spectrums if s is not None]

    number_of_peaks = [len(spec.peaks) for spec in processed_spectrums]
    print('## AFTER PROCESSING ##')
    print("Maximum number of peaks in one spectrum:", np.max(number_of_peaks))
    print("Number of spectra with > 1000 peaks:", np.sum(np.array(number_of_peaks)>1000))
    print("Number of spectra with > 2000 peaks:", np.sum(np.array(number_of_peaks)>2000))
    print("Number of spectra with > 5000 peaks:", np.sum(np.array(number_of_peaks)>5000))
    print("Careful: Number of spectra with < 10 peaks:", np.sum(np.array(number_of_peaks)<10))
    
    # Save in a pickle object
    #pickle.dump(processed_spectrums, open("TOTAL_CFMID_MGF.matchms_processed.pickle", "wb"))
    #pickle.dump(processed_spectrums, open("GNPS_CFMID.min10peaks.mgf.matchms_processed.pickle", "wb"))
    pickle.dump(processed_spectrums, open("David_MGF.min10peaks.mgf.matchms_processed.pickle", "wb"))

if __name__ == "__main__":
    
    # Number of processors 
    NPROC = 32
    
    # Read spectrums from a MGF formatted file, for other formats see https://matchms.readthedocs.io/en/latest/api/matchms.importing.html
    #spectrums = load_from_mgf("/home/dpflieger/Project/Gaquerel/Elser/Databases/David_MGF/TOTAL_COMPOUNDS_DB.energies_merged_name.mgf")
    #spectrums = load_from_mgf("/home/dpflieger/tmp/matchms2/GNPS_CFMID.mgf")
    spectrums = load_from_mgf("/home/dpflieger/Project/Gaquerel/Elser/Databases/David_MGF/alltissues27042021-py_correct.mgf")

    # number_of_peaks = [len(spec.peaks) for spec in spectrums]
    # print('## BEFORE PROCESSING ##')
    # print("Maximum number of peaks in one spectrum:", np.max(number_of_peaks))
    # print("Number of spectra with > 1000 peaks:", np.sum(np.array(number_of_peaks)>1000))
    # print("Number of spectra with > 2000 peaks:", np.sum(np.array(number_of_peaks)>2000))
    # print("Number of spectra with > 5000 peaks:", np.sum(np.array(number_of_peaks)>5000))
    # print("Careful: Number of spectra with < 5 peaks:", np.sum(np.array(number_of_peaks)<5))

    main()

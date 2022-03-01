import pickle
import gensim
import rdkit

import numpy as np
import multiprocessing
from functools import partial

from time import time
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from functools import wraps
from joblib import Parallel, delayed
from scipy import spatial

from matchms import calculate_scores
from matchms.similarity import CosineGreedy
from matchms.similarity import ModifiedCosine
from matchms.similarity.spectrum_similarity_functions import collect_peak_pairs
from matchms.filtering import normalize_intensities
from matchms.filtering import select_by_relative_intensity
from matchms.filtering import select_by_mz

from spec2vec import Spec2Vec
from spec2vec import SpectrumDocument
from spec2vec.model_building import train_new_word2vec_model

import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt

NPROC = 16

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

def chunks_index(length, ncores):
    njobs = ncores-1 
    n = round(length/njobs)
    i, ind = 0, []
    while i < length-1:
        if(i+n < length-1):
            ind.append((i, i+n-1))
        else:
            ind.append((i, length-1))
        i += n
    return(ind)

def plot_spectra_comparison(spectrum1_in, spectrum2_in,
                            model,
                            intensity_weighting_power=0.5,
                            num_decimals=2,
                            min_mz=5,
                            max_mz=500,
                            intensity_threshold=0.01,
                            method='cosine',
                            tolerance=0.005,
                            wordsim_cutoff=0.5,
                            circle_size=5,
                            circle_scaling='wordsim',
                            padding=10,
                            display_molecules="smiles",
                            figsize=(12, 12),
                            filename=None):
    """ In-depth visual comparison of spectral similarity scores,
    calculated based on cosine/mod.cosine and Spev2Vec.
    Parameters
    ----------
    method: str
        'cosine' or 'modcos' (modified cosine score)
    circle_scaling: str
        Scale circles based on 'wordsim' or 'peak_product'
    """

    def apply_filters(s):
        s = normalize_intensities(s)
        s = select_by_mz(s, mz_from=min_mz, mz_to=max_mz)
        s = select_by_relative_intensity(s, intensity_from=intensity_threshold)
        s.losses = None
        return s

    spectrum1 = apply_filters(spectrum1_in)
    spectrum2 = apply_filters(spectrum2_in)

    plt.style.use("seaborn-white")#('ggplot')
    plot_colors = ['darkcyan', 'purple']

    # Definitions for the axes
    left, width = 0.1, 0.6
    bottom, height = 0.1, 0.6
    spacing = 0.01

    rect_wordsim = [left, bottom, width, height]
    rect_specx = [left, bottom + height + spacing, width, 0.2]
    rect_specy = [left + width + spacing, bottom, 0.25, height]

    document1 = SpectrumDocument(spectrum1, n_decimals=num_decimals)
    document2 = SpectrumDocument(spectrum2, n_decimals=num_decimals)

    # Remove words/peaks that are not in dictionary
    select1 = np.asarray([i for i, word in enumerate(document1.words) if word in model.wv.vocab])
    select2 = np.asarray([i for i, word in enumerate(document2.words) if word in model.wv.vocab])
    peaks1 = np.asarray(spectrum1.peaks[:]).T
    peaks2 = np.asarray(spectrum2.peaks[:]).T
    peaks1 = peaks1[select1, :]
    peaks2 = peaks2[select2, :]
    min_peaks1 = np.min(peaks1[:, 0])
    min_peaks2 = np.min(peaks2[:, 0])
    max_peaks1 = np.max(peaks1[:, 0])
    max_peaks2 = np.max(peaks2[:, 0])
    possible_grid_points = np.arange(0, 2000, 50)
    grid_points1 = possible_grid_points[(possible_grid_points > min_peaks1 - padding) \
                                        & (possible_grid_points < max_peaks1 + padding)]
    grid_points2 = possible_grid_points[(possible_grid_points > min_peaks2 - padding) \
                                        & (possible_grid_points < max_peaks2 + padding)]

    word_vectors1 = model.wv[[document1.words[x] for x in select1]]
    word_vectors2 = model.wv[[document2.words[x] for x in select2]]

    csim_words = 1 - spatial.distance.cdist(word_vectors1, word_vectors2, 'cosine')
    csim_words[csim_words < wordsim_cutoff] = 0  # Remove values below cutoff
    print(np.min(csim_words), np.max(csim_words))

    # Plot spectra
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=figsize)
    # Word similariy plot (central)
    ax_wordsim = plt.axes(rect_wordsim)
    ax_wordsim.tick_params(direction='in', top=True, right=True)
    # Spectra plot (top)
    ax_specx = plt.axes(rect_specx)
    ax_specx.tick_params(direction='in', labelbottom=False)
    # Spectra plot 2 (right)
    ax_specy = plt.axes(rect_specy)
    ax_specy.tick_params(direction='in', labelleft=False)

    # Spec2Vec similarity plot:
    # -------------------------------------------------------------------------
    data_x = []
    data_y = []
    data_z = []
    data_peak_product = []
    for i in range(len(select1)):
        for j in range(len(select2)):
            data_x.append(peaks1[i, 0])
            data_y.append(peaks2[j, 0])
            data_z.append(csim_words[i, j])
            data_peak_product.append(peaks1[i, 1] * peaks2[j, 1])

    # Sort by word similarity
    data_x = np.array(data_x)
    data_y = np.array(data_y)
    data_z = np.array(data_z)
    data_peak_product = np.array(data_peak_product)
    idx = np.lexsort((data_x, data_y, data_z))

    cm = plt.cm.get_cmap('RdYlBu_r')  # 'YlOrRd') #'RdBu_r')

    # Plot word similarities
    if circle_scaling == 'peak_product':
        wordsimplot = ax_wordsim.scatter(data_x[idx],
                                         data_y[idx],
                                         s=100 * circle_size *
                                         (0.01 + data_peak_product[idx]**2),
                                         marker="o",
                                         c=data_z[idx],
                                         cmap=cm,
                                         alpha=0.6)
    elif circle_scaling == 'wordsim':
        wordsimplot = ax_wordsim.scatter(data_x[idx],
                                         data_y[idx],
                                         s=100 * circle_size *
                                         (0.01 + data_z[idx]**2),
                                         marker="o",
                                         c=data_z[idx],
                                         cmap=cm,
                                         alpha=0.6)

    # (Modified) Cosine similarity plot:
    # -------------------------------------------------------------------------
    if method == 'cosine':
        score_classical, used_matches = cosine_score(spectrum1, spectrum2, tolerance, modified_cosine=False)
    elif method == 'modcos':
        score_classical, used_matches = cosine_score(spectrum1, spectrum2, tolerance, modified_cosine=True)
    else:
        print("Given method unkown.")

    idx1, idx2, _ = zip(*used_matches)
    cosine_x = []
    cosine_y = []
    for i in range(len(idx1)):
        if idx1[i] in select1 and idx2[i] in select2:
            cosine_x.append(peaks1[int(idx1[i]), 0])
            cosine_y.append(peaks2[int(idx2[i]), 0])

    # Plot (mod.) cosine similarities
    ax_wordsim.scatter(cosine_x, cosine_y, s=100, c='black', marker=(5, 2))
    ax_wordsim.set_xlim(min_peaks1 - padding, max_peaks1 + padding)
    ax_wordsim.set_ylim(min_peaks2 - padding, max_peaks2 + padding)
    ax_wordsim.set_xlabel('spectrum 1 - fragment mz', fontsize=16)
    ax_wordsim.set_ylabel('spectrum 2 - fragment mz', fontsize=16)
    ax_wordsim.tick_params(labelsize=13)
    ax_wordsim.set_xticks(grid_points1)
    ax_wordsim.set_yticks(grid_points2)
    ax_wordsim.grid(True)

    # Plot spectra 1
    ax_specx.vlines(peaks1[:, 0], [0], peaks1[:, 1], color=plot_colors[0])
    ax_specx.plot(peaks1[:, 0], peaks1[:, 1], '.')  # Stem ends
    ax_specx.plot([peaks1[:, 0].max(), peaks1[:, 0].min()], [0, 0],
                  '--')  # Middle bar
    ax_specx.set_xlim(min_peaks1 - padding, max_peaks1 + padding)
    ax_specx.set_yticks([0,0.25,0.5,0.75,1])
    ax_specx.set_xticks(grid_points1)
    ax_specx.set_ylabel('peak intensity (relative)', fontsize=16)
    ax_specx.tick_params(labelsize=13)

    ax_specx.grid(True)

    # Plot spectra 2
    ax_specy.hlines(peaks2[:, 0], [0], peaks2[:, 1], color=plot_colors[1])
    ax_specy.plot(peaks2[:, 1], peaks2[:, 0], '.')  # Stem ends
    ax_specy.plot([0, 0], [peaks2[:, 0].min(), peaks2[:, 0].max()],
                  '--')  # Middle bar
    ax_specy.set_ylim(min_peaks2 - padding, max_peaks2 + padding)
    ax_specy.set_xticks([0,0.25,0.5,0.75,1])
    ax_specy.set_yticks(grid_points2)
    ax_specy.set_xlabel('peak intensity (relative)', fontsize=16)
    ax_specy.tick_params(labelsize=13)

    ax_specy.grid(True)

    fig.colorbar(wordsimplot, ax=ax_specy)
    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()

    # Plot molecules
    # -------------------------------------------------------------------------
    # TODO: Add conversion from InchI or 
    if display_molecules == "smiles":
        smiles = [spectrum1.get("smiles"), spectrum2.get("smiles")]
        molecules = [rdkit.Chem.MolFromSmiles(x) for x in smiles]
        display(rdkit.Draw.MolsToGridImage(molecules, molsPerRow=2, subImgSize=(400, 400)))
    elif display_molecules == "inchikey":
        inchikeys = [spectrum1.get("inchikey"), spectrum2.get("inchikey")]
        molecules = [rdkit.Chem.inchi.MolFromInchi(x) for x in inchikeys]
        display(rdkit.Draw.MolsToGridImage(molecules, molsPerRow=2, subImgSize=(400, 400)))

def cosine_score(spectrum1, spectrum2, tolerance, modified_cosine=False):
    """
    Parameters
    ----------
    spectrum1 : TYPE
        DESCRIPTION.
    spectrum2 : TYPE
        DESCRIPTION.
    tolerance : TYPE
        DESCRIPTION.
    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    def get_peaks_arrays():
        """Get peaks mz and intensities as numpy array."""
        spec1 = np.vstack((spectrum1.peaks.mz, spectrum1.peaks.intensities)).T
        spec2 = np.vstack((spectrum2.peaks.mz, spectrum2.peaks.intensities)).T
        assert max(spec1[:, 1]) <= 1, ("Input spectrum1 is not normalized. ",
                                       "Apply 'normalize_intensities' filter first.")
        assert max(spec2[:, 1]) <= 1, ("Input spectrum2 is not normalized. ",
                                       "Apply 'normalize_intensities' filter first.")
        return spec1, spec2

    def get_matching_pairs():
        """Get pairs of peaks that match within the given tolerance."""
        zero_pairs = collect_peak_pairs(spec1, spec2, tolerance, shift=0.0)
        if modified_cosine:
            message = "Precursor_mz missing. Apply 'add_precursor_mz' filter first."
            assert spectrum1.get("precursor_mz") and spectrum2.get("precursor_mz"), message
            mass_shift = spectrum1.get("precursor_mz") - spectrum2.get("precursor_mz")
            nonzero_pairs = collect_peak_pairs(spec1, spec2, tolerance, shift=mass_shift)
            unsorted_matching_pairs = zero_pairs + nonzero_pairs
        else:
            unsorted_matching_pairs = zero_pairs
        return sorted(unsorted_matching_pairs, key=lambda x: x[2], reverse=True)

    def calc_score():
        """Calculate cosine similarity score."""
        used1 = set()
        used2 = set()
        score = 0.0
        used_matches = []
        for match in matching_pairs:
            if not match[0] in used1 and not match[1] in used2:
                score += match[2]
                used1.add(match[0])  # Every peak can only be paired once
                used2.add(match[1])  # Every peak can only be paired once
                used_matches.append(match)
        # Normalize score:
        score = score/max(np.sum(spec1[:, 1]**2), np.sum(spec2[:, 1]**2))
        return score, used_matches

    spec1, spec2 = get_peaks_arrays()
    matching_pairs = get_matching_pairs()
    return calc_score()

if __name__ == "__main__":
    
    # My model trained from CFMID mgf
    #model = gensim.models.Word2Vec.load("spec2vec_AllPositive_ratio05_filtered_201101_iter_15.model")
    #model = gensim.models.Word2Vec.load("GNPS_CFMID.min5peaks.model")
    model = gensim.models.Word2Vec.load("GNPS_CFMID.min10peaks.model")
     
    #similarity_measure = ModifiedCosine(tolerance=0.005)
    similarity_measure = Spec2Vec(model=model, intensity_weighting_power=0.5, allowed_missing_percentage=100.0)
        
    print("Loading spectra from pickle objects")
    #with open("GNPS_CFMID.min10peaks.mgf.matchms_processed.pickle", "rb") as cfmid, open("David_MGF.min10peaks.mgf.matchms_processed.pickle", "rb") as user:
    with open("TOTAL_CFMID_MGF.matchms_processed.pickle", "rb") as cfmid, open("David_MGF.min10peaks.mgf.matchms_processed.pickle", "rb") as user:
        cfmid_spectrums, user_spectrums = pickle.load(cfmid), pickle.load(user)
    
    # Create spectrum documents
    #print("Creating spectrum document")
    #reference_documents = [SpectrumDocument(s, n_decimals=2) for s in cfmid_spectrums]

    # print("Creating word2vec model")
    # model_file = "GNPS_CFMID.min10peaks.model"
    # model = train_new_word2vec_model(reference_documents, iterations=[10, 20, 30], filename=model_file, workers=24, progress_logger=True)
    # print("Model done!")
    # exit()

    # spectrum = user_spectrums[0]
    # print(spectrum.peaks.mz[0])
    # print(spectrum.peaks.intensities[0])
    # print(spectrum.get('title'))

    # csim = plot_spectra_comparison(user_spectrums[1], user_spectrums[1],
    #                             model,
    #                             intensity_weighting_power=0.5,
    #                             num_decimals=2,
    #                             min_mz=50,
    #                             max_mz=2000,
    #                             intensity_threshold=0.05,
    #                             method="modcos",#"cosine", "modcos", #
    #                             tolerance=0.005,
    #                             wordsim_cutoff=0.05,
    #                             circle_size=5,
    #                             circle_scaling='wordsim',
    #                             padding=30,
    #                             display_molecules=False,
    #                             figsize=(12, 12),
    #                             filename="TEST.pdf")#None)#

   
    
    #print(spectrums[5]) # <matchms.Spectrum.Spectrum object at 0x7fea7408a0d0>

    # for i, spectrum in enumerate(spectrums):
    #     print(i, spectrum)
    #     scores = calculate_scores([spectrum], spectrums, similarity_measure, is_symmetric=False)
    #     best_matches = scores.scores_by_query(spectrum, sort=True)[:5]
    #     print("Best matches are", best_matches)

    #for i, spectrum in enumerate(spectrums):
        #print(i, spectrum)
        #scores = calculate_scores(spectrums[0:10000], spectrums[0:10000], similarity_measure, is_symmetric=False)
    
    t = time()
    number_of_spectra = len(list(cfmid_spectrums))
    print(f"Computing scores with spec2vec for {number_of_spectra} spectra")
    
    parts = chunks_index(length = number_of_spectra, ncores = 64)

    def get_scores(parts):
        start, end = parts
        scores = calculate_scores(cfmid_spectrums[start:end], user_spectrums, similarity_measure, is_symmetric=False)
        return scores

    # Using joblib package
    #scores = Parallel(n_jobs=NPROC)(delayed(calculate_scores)(cfmid_spectrums, user_spectrums[start:end], similarity_measure, is_symmetric=False) for start, end in chunks_index(length = number_of_spectra, ncores = NPROC))

    # Without parallelization 
    scores = calculate_scores(cfmid_spectrums, user_spectrums, similarity_measure, is_symmetric=False)

    # # Using threading package
    # with multiprocessing.Pool(processes = NPROC) as pool:
    #     scores = pool.map(get_scores, parts)

    print("Scoring done. Saving now.")
    # Save scores with pickle 
    #pickle.dump(scores, open("spec2vec_scores.pickle", "wb"))
    pickle.dump(scores, open("spec2vec_scores_CFMID_vs_DAVID_MGF.pickle", "wb"))

    #scores = process_map(get_scores, parts, max_workers = NPROC)

    print("Done in", round(time() - t), "seconds", sep = " ")

    #for i, spectrum in enumerate(user_spectrums[0:10000]):
    # best_scores = scores.scores_by_query(user_spectrums[0], sort=True)[:5]
    # print("Best matches are", best_scores)
    # best_spectrums = [x[0] for x in best_scores]
    # print([x.get('title') for x in best_spectrums])

    # fig = plt.figure(figsize=(6,6), dpi=150)
    # plt.imshow(scores.scores[:50, :50]["score"], cmap="viridis")
    # plt.colorbar(shrink=0.7)
    # plt.title("Modified Cosine spectra similarities")
    # plt.xlabel("Spectrum #ID")
    # plt.ylabel("Spectrum #ID")
    # fig.tight_layout()
    # fig.savefig("matrix.pdf")


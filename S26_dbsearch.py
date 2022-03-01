from matchms import calculate_scores
from matchms.similarity import CosineGreedy
from matchms.similarity import ModifiedCosine
import os
import sys
from matchms.importing import load_from_mgf
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import select_by_intensity
from matchms.filtering import select_by_mz
from matchms.filtering import default_filters
from matchms.filtering import repair_inchi_inchikey_smiles
from matchms.filtering import derive_inchikey_from_inchi
from matchms.filtering import derive_smiles_from_inchi
from matchms.filtering import derive_inchi_from_smiles
from matchms.filtering import harmonize_undefined_inchi
from matchms.filtering import harmonize_undefined_inchikey
from matchms.filtering import harmonize_undefined_smiles

def peak_processing(spectrum):
    spectrum = default_filters(spectrum)
    spectrum = normalize_intensities(spectrum)
    spectrum = select_by_intensity(spectrum, intensity_from=0.01)
    spectrum = select_by_mz(spectrum, mz_from=10, mz_to=2000)
    return spectrum

def metadata_processing(spectrum):
    spectrum = default_filters(spectrum)
    spectrum = repair_inchi_inchikey_smiles(spectrum)
    spectrum = derive_inchi_from_smiles(spectrum)
    spectrum = derive_smiles_from_inchi(spectrum)
    spectrum = derive_inchikey_from_inchi(spectrum)
    spectrum = harmonize_undefined_smiles(spectrum)
    spectrum = harmonize_undefined_inchi(spectrum)
    spectrum = harmonize_undefined_inchikey(spectrum)
    return spectrum



path_data = "C:/Users/Public/database/"  # enter path to library mgf file
file_mgf = os.path.join(path_data,"{}.mgf".format(str(sys.argv[1])))
library = list(load_from_mgf(file_mgf))

#print(len(library))
library = [peak_processing(s) for s in library]
#library = [metadata_processing(s) for s in library]

#print(len(library))


path_data = "C:/Users/Public/input/"  # enter path to Measured spectra mgf file
file_mgf = os.path.join(path_data,"{}".format(str(sys.argv[2])))
spectrums = list(load_from_mgf(file_mgf))


#print(len(spectrums))
spectrums = [peak_processing(s) for s in spectrums]
spectrums = [metadata_processing(s) for s in spectrums]
#print(len(spectrums))

similarity_measure = ModifiedCosine(tolerance=0.005)
scores = calculate_scores(library, spectrums, similarity_measure,
                          is_symmetric=False)
scores.scores[:5, :5]["matches"]

counter=0
with open("C:\\Users\\Public\\results\\p{}-results.tsv".format(str(sys.argv[1])), 'w') as f:
    f.write('Feature_id''\t''m/z')
    f.write('\t''Library Hit 1''\t''m/z 1''\t''cousine-score Hit 1''\t''numoffpeaks 1''\t''Smiles 1''\t')
    f.write('\n')

    for i in scores.queries:
        #print(i.get('feature_id'))
        f.write(i.get('feature_id')) #print feature_id to file
        f.write('\t')
        #print(dir(i))
        x=i.get('pepmass')
        f.write(str(x[0])) #print m\z to file
        f.write('\t')
        x=scores.scores_by_query(spectrums[counter], sort=True)[:1] # get top 1 hit
        #print(x[0][0].get('title'))
        for a in x:                                # get metadata of top 10 hits an print to file
            f.write(a[0].get('title')) #print Name
            f.write('\t')
            y=a[0].get('pepmass') 
            f.write(str(y[0])) # print m/z
            f.write('\t')
            #print(b[1]) #get hit values and matching peaks of top 10
            f.write(str(a[1][0].item())) #cousin score
            f.write('\t')
            f.write(str(a[1][1].item())) #matching peaks
            f.write('\t')
            f.write(a[0].get('inchi')) #print smiles
            f.write('\t')
        counter+=1
        #print(counter)
    
        f.write('\n')
f.close()
print("wrote file to p{}-results.csv".format(str(sys.argv[1])))
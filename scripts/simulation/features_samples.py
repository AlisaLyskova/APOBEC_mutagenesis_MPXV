import pandas as pd

def my_filtering_function(pair):
    key, value = pair
    if "dRNA" in value or "mock" in value:
        return False  # filter pair out of the dictionary
    else:
        return True  # keep pair in the filtered dictionary
        
## motifs
motif_list = ['TCT', 'TCA', 'TCC', 'TCG']

## samples features
samples_features = pd.read_excel("41597_2023_2149_MOESM2_ESM.xlsx", skiprows=2)
samples_features = samples_features.fillna(method='ffill', axis=0)

samples_features_dict = samples_features.set_index('Run accession')['File names'].to_dict()

# filter mock and dRNA
filtered_samples_features = dict(filter(my_filtering_function, samples_features_dict.items()))

# LIST WITH FEATURES
features = []
for key, value in filtered_samples_features.items():
    new_value = value.split('_')[2:]
    new_value[2] = new_value[2][5:]
    for motif in motif_list:
        new_feature = key+'_'+'_'.join(new_value)+'_'+motif
        features.append(new_feature)

print(' '.join(features))

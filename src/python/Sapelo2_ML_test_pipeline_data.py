import os
import sys
from collections import Counter
from pathlib import Path

import joblib  # # for save the model
import numpy as np
import pandas as pd

# %%
import sklearn
from scipy import stats
from six import StringIO
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, train_test_split

### this script is used to use the model we train and give the label to the mutations data
### this test dataset doesn't contains true labels

# %%
model_folder = sys.argv[1]
data_test_file = sys.argv[2]
final_output = sys.argv[3]


model_file = (
    model_folder
    + "/"
    + "WithVACutNewlabel_MT_overlap_comparison_Human_remained_02_09_23.joblib"
)

### because each number represent different meaning, and we need to re-encode to it's orignal meaning
## 02_09_23, we change the encoding dict because lack of HQ WT
reverse_results_encoding_dict_file = (
    model_folder + "/" + "New_label_Model_encoding_results_02_09.txt"
)
reverse_results_encoding_dict = pd.read_csv(
    reverse_results_encoding_dict_file, sep="\t"
)
reverse_results_encoding_dict = dict(
    zip(reverse_results_encoding_dict.Encoding, reverse_results_encoding_dict.Results)
)

needed_column = ["Chrom", "Pos", "Ref", "Alt", "Ref_reads", "Alt_reads", "VAF"]
### load the data and prepare for model testing

somatic_test_info = pd.read_csv(data_test_file, sep="\t")
### only grep somatic results from the our previous pipeline filtering
somatic_test_info = somatic_test_info[somatic_test_info.Status == "Somatic"]
col_pos = np.where(np.array(somatic_test_info.columns) == "Status")[0][0]
col_name = list(somatic_test_info.columns)
col_name[col_pos] = "Pipeline_Results"
col_pos = np.where(np.array(somatic_test_info.columns) == "Start")[0][0]
col_name[col_pos] = "Pos"
somatic_test_info.columns = col_name
somatic_test_data = somatic_test_info[needed_column]
test_data = somatic_test_data
test_data = test_data.drop_duplicates()
test_data_meta = test_data.copy()

# %%
#### Encode the test set for the model, the label encoding number must be the same as the training (ex: results_encoding_dict)
#### 10/18/22, test data only test with the model, and we don't know the ground true of the data (RNA-seq MT, OSA, OM dataset)
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

le = LabelEncoder()
test_data.loc[:, "Chrom"] = le.fit_transform(test_data.loc[:, "Chrom"])
test_data.loc[:, "Pos"] = test_data.loc[:, "Pos"].astype(int)
test_data.loc[:, "Ref"] = le.fit_transform(test_data.loc[:, "Ref"])
test_data.loc[:, "Alt"] = le.fit_transform(test_data.loc[:, "Alt"])

# test_data_meta.loc[:,"Encoding_results"] = le.fit_transform(test_data.iloc[:, 7])
X = test_data.iloc[:, 0:7].to_numpy()
# y = test_data.Results.values
# le.fit(y)
# y= le.transform(y)

# %%
## load the model and see the test
loaded_model = joblib.load(model_file)
result = loaded_model.predict(X)


# %%
#### merge the predict results with the original df
def mapBackResults(each_value):
    results = reverse_results_encoding_dict[each_value]
    return str(results)


### merge back the predicted results with meta_data
test_data_meta["Model_prediction"] = result
test_data_meta.loc[:, "Model_prediction"] = test_data_meta["Model_prediction"].apply(
    mapBackResults
)
overlap_col = list(set(somatic_test_info.columns) & set(test_data_meta.columns))
final_test_data = pd.merge(
    somatic_test_info, test_data_meta, how="left", on=overlap_col
)

### remove unnecessary columns
remove_col = [
    "Chrom_mut_info",
    "Transcript_mut_info",
    "Homo_Hetero_Ratio",
    "Gene_mut_info",
    "Hetero_ratio",
    "Variant_in_all_sample_ratio",
    "Variants_in_number_tumor_type",
    "Each_tumor_type_total_sample",
    "Each_tumor_type_ratio",
    "Sample_number",
    "Homo_ratio",
    "Variants_in_each_tumor_type",
]
final_test_data = final_test_data.drop(remove_col, axis=1)
final_test_data = final_test_data.drop_duplicates()
final_test_data.to_csv(final_output, sep="\t", index=False)

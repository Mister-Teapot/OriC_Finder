import pandas as pd
import numpy as np
from sklearn.svm import SVC
import joblib

all_samples_df = pd.read_csv('Hyperparameter tuning/tuning.csv')
all_refseqs_df = pd.read_csv('DoriC data prep/DoriC_oriC_concat_entries.csv')['RefSeq']

train_refseqs = all_refseqs_df.sample(frac = 0.75, random_state=42)
test_refseqs = all_refseqs_df.drop(train_refseqs.index)
train_samples_refseqs = []
for i, sample in all_samples_df.iterrows():
    if sample['RefSeq_oriC'][:-2] not in test_refseqs:
        train_samples_refseqs.append(sample['RefSeq_oriC'])

# .to_numpy() to get rid of feature name warning
X_train = all_samples_df[all_samples_df['RefSeq_oriC'].isin(train_samples_refseqs)][['Z_occurance', 'D_occurance']].to_numpy()
y_train = all_samples_df[all_samples_df['RefSeq_oriC'].isin(train_samples_refseqs)]['Correct'].to_numpy()

model = SVC(C=500, kernel='rbf', random_state=42).fit(X_train, y_train)

joblib.dump(model, '75_train_model_no_G.pkl')

from sklearn.linear_model import LinearRegression, LogisticRegressionCV
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier

from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split

from scipy.stats import spearmanr

import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from ast import literal_eval
import math
import sys

sys.path.append('../OriC_Finder')
from peak import Peak

CSV_PATH = 'Comparison/v5/in_both_sets_all.csv'
CSV_OUT  = 'Hyperparameter tuning/tuning.csv'

def remake_comparator_df(csv_path):
    df = pd.read_csv(csv_path)

    new_df = {
        'RefSeq_oriC': [],
        'Z_occurance': [],
        'G_occurance': [],
        'D_occurance': [],
        'Correct'    : [],
    }

    for i, sample_original in df.iterrows():
        sample = sample_original.copy()
        sample.dropna(inplace=True)
        seq_len = sample['Sequence_length']

        D_oriC_cols = sorted( [i for i in sample.axes[0] if 'DoriC_oriC' in i], key=lambda x:x[-1] )
        D_oriC_middles = [ Peak.get_middle( int( literal_eval(sample[D_oriC_col])[0] ), int( literal_eval(sample[D_oriC_col])[1] ), seq_len ) for D_oriC_col in D_oriC_cols ]

        num_of_Z_oriCs = max([int(i[-1]) for i in sample.axes[0] if 'oriC_middles' in i]) + 1

        for i in range(num_of_Z_oriCs):
            Z_pos = sample[f'oriC_middles_{i}']
            dist = float('inf')
            for D_pos in D_oriC_middles:
                new_dist = Peak.calc_dist(Z_pos, D_pos, seq_len)
                if new_dist < dist:
                    dist = new_dist
            new_df['RefSeq_oriC'].append(sample['RefSeq'] + f'_{i}')
            new_df['Z_occurance'].append(sample[f'Z_Occurance_oriC_{i}'])
            new_df['G_occurance'].append(sample[f'G_Occurance_oriC_{i}'])
            new_df['D_occurance'].append(sample[f'D_Occurance_oriC_{i}'])
            # if dist <= (seq_len / 100) * 2.5:
            new_df['Correct'].append( dist <= (seq_len / 100) * 2.5 )
            # else:
            #     new_df['Correct'].append('false')
    return pd.DataFrame(new_df)

def polyfit(x, y, degree):
    results = {}
    coeffs = np.polyfit(x, y, degree)
    p = np.poly1d(coeffs)
    yhat = p(x)
    ybar = np.sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((y - ybar)**2)
    results['formula'] = str(p)
    results['r_squared'] = ssreg / sstot

    return results

if __name__ == '__main__':
    # df = remake_comparator_df(CSV_PATH)
    # df.to_csv(CSV_OUT)
    df = pd.read_csv(CSV_OUT)
    Z = df['Z_occurance']
    G = df['G_occurance']
    D = df['D_occurance']

    # cov_mat = np.cov(Z, G, rowvar=False)
    # print(cov_mat)
    # pearson = cov_mat / (Z.std() * G.std())
    # print(pearson)
    # spearman = spearmanr(Z, G)
    # print(spearman)

    # sns.set_theme(style="dark")
    # f, ax = plt.subplots(figsize=(7, 5))
    # sns.scatterplot(x=Z, y=G, s=5, color=".15")
    # sns.histplot(data=df[['Z_occurance', 'G_occurance', 'D_occurance']], bins=50)
    # sns.kdeplot(data=df[['Z_occurance', 'D_occurance']], levels=5, color="b", linewidths=1, fill=True)

    # plt.xlabel('Confidence (%)')
    # plt.xlim(0, 1)
    # plt.show()

    X = df[['Z_occurance', 'G_occurance', 'D_occurance']]
    y = df['Correct']
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.8, random_state=0)

    # X2 = sm.add_constant(X)
    # est = sm.OLS(y, X2)
    # est2 = est.fit()
    # print(est2.summary())
    def sigmoid(x):
        return np.where(x < 0, np.exp(x)/(1 + np.exp(x)), 1/(1 + np.exp(-x)))

    lin_model = LogisticRegressionCV(random_state=0).fit(X_train, y_train)
    one, two, three = 0.4, 0.5, 0.5
    # print(dir(lin_model))
    print('model scoring  :', lin_model.predict_proba( pd.DataFrame({'Z_occurance': [one], 'G_occurance':[two], 'D_occurance': [three]}) )[0][-1])
    print('model decision :', lin_model.predict( pd.DataFrame({'Z_occurance': [one], 'G_occurance':[two], 'D_occurance': [three]}) )[0])
    print('my scoring     :', sigmoid( np.dot([one, two, three], lin_model.coef_.T) + lin_model.intercept_ )[0])
    print('my decision    :', sigmoid( np.dot([one, two, three], lin_model.coef_.T) + lin_model.intercept_ )[0] > 0.5)
    
    # y_pred_lin = lin_model.predict(X_test)

    # print('coef_     :', lin_model.coef_)
    # print('intercept_:', lin_model.intercept_)

    # print('Linear Model Attributes')
    # # print(model.get_params())
    # print('features     :', lin_model.feature_names_in_)
    # print('coefficients :', lin_model.coef_)
    # print('accuracy     :', accuracy_score(y_test, y_pred_lin))
    # print()

    # forest = RandomForestClassifier(random_state=0).fit(X_train, y_train)
    # y_pred_rf = forest.predict(X_test)

    # print('Random Forest Attributes')
    # print('features     :', forest.feature_names_in_)
    # print('coefficients :', forest.feature_importances_)
    # print('accuracy     :', accuracy_score(y_test, y_pred_rf))

'''
    RefSeq      Z_oc    G_oc    D_oc    Correct
    NC_000913   0.8     0.4     0.7     True
    NC_000964   0.5     0.9     0.3     False
    ...
'''
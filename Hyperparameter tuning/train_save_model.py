import pandas as pd
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import joblib

df = pd.read_csv('Hyperparameter tuning/tuning.csv')

# .to_numpy() to get rid of feature name warning
X = df[['Z_occurance', 'G_occurance', 'D_occurance']].to_numpy()
y = df['Correct'].to_numpy()

X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.7, random_state=42)

model = SVC(C=500, kernel='rbf', random_state=42).fit(X_train, y_train)

y_pred_model = model.predict(X_test)
y_pred_dummy = [False] * X_test.shape[0]

print('model acc', accuracy_score(y_test, y_pred_model))
print('dummy_F acc', accuracy_score(y_test, y_pred_dummy))

use_model = SVC(C=500, kernel='rbf', random_state=42).fit(X, y)

joblib.dump(use_model, 'model.pkl')

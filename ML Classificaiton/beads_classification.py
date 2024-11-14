import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, accuracy_score
from imblearn.over_sampling import SMOTE
from sklearn.neighbors import KNeighborsClassifier

# Load both datasets
data1 = pd.read_csv('Path/to/Dataset1')
data2 = pd.read_csv('Path/to/Dataset1')

# Calculate Peak-to-Peak Voltage for each dataset
data1['PeakToPeakVoltage'] = data1['PosPeak'] %- data1['NegPeak']
data2['PeakToPeakVoltage'] = data2['PosPeak'] %- data2['NegPeak']

# Filter for Channel 1 data only (as an example)
data1_filtered = data1[data1['Channel'] == 2][['PeakToPeakVoltage', 'Time_ms']].copy()
data2_filtered = data2[data2['Channel'] == 2][['PeakToPeakVoltage', 'Time_ms']].copy()

# Label each dataset
data1_filtered['Label'] = 0  # Label 0 for 4um dataset
data2_filtered['Label'] = 1  # Label 1 for 7um dataset

# Combine both datasets
combined_data = pd.concat([data1_filtered, data2_filtered])

# Separate features and labels
X = combined_data[['PeakToPeakVoltage', 'Time_ms']]
y = combined_data['Label']

# Standardize the features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

smote = SMOTE(random_state=42)
X_resampled, y_resampled = smote.fit_resample(X_scaled, y)

# Then, split the resampled data
X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.3, random_state=42)

# Initialize and train a Random Forest/SVM/kNN/or other types of classifier
# classifier = RandomForestClassifier(class_weight='balanced', random_state=42)
classifier = KNeighborsClassifier(n_neighbors=5)
classifier.fit(X_train, y_train)

# Predict on the test set
y_pred = classifier.predict(X_test)

# Calculate accuracy
accuracy = accuracy_score(y_test, y_pred)
print(f'Accuracy: {accuracy:.2f}')

# Generate confusion matrix
conf_matrix = confusion_matrix(y_test, y_pred)
ConfusionMatrixDisplay(confusion_matrix=conf_matrix, display_labels=['4um', '7um']).plot()
# plt.title("Confusion Matrix for Classifying 4um and 7um Particles")
plt.show()
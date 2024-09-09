import pandas as pd
from sklearn.linear_model import LinearRegression
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import PolynomialFeatures

k = 1.96
MSE = 0.029188325596117455
corr_fact = 1.0000153508435288
diff_mass_sum = 21495293839.22711
mono_mean = 905.1084025286648
beta0 = -0.06034476939748856
beta1 = 0.99942198


def PredNomMass(monoMass, k, MSE, corr_fact, diff_mass_sum):
    
    nomMass = beta0 + beta1 * monoMass
    #nomMass = beta0 + beta1[0] * monoMass + beta2[0] * monoMass**2

    rounded_nomMass = np.round(nomMass)

    lwb = nomMass - k * (MSE)**0.5 * (corr_fact + (monoMass - mono_mean)**2 / diff_mass_sum)**0.5
    upb = nomMass + k * (MSE)**0.5 * (corr_fact + (monoMass - mono_mean)**2 / diff_mass_sum)**0.5

    return rounded_nomMass, nomMass, lwb, upb




# Calculated using this :
'''
df = pd.read_csv('database_folder/Lipid_Maps_and_SWISS_lipids.csv')
X = df['mass1'].values.reshape(-1, 1)
y = df['nominal_mass'].values

degree = 1
poly = PolynomialFeatures(degree=degree)
X_poly = poly.fit_transform(X)
model = LinearRegression()
model.fit(X_poly, y)

beta0 = model.intercept_
beta1 = model.coef_[1:]  
#beta2 = model.coef_[2:]

y_pred = model.predict(X_poly)

#print(f"Intercept: {beta0}")
#print(f"beta1: {beta1}")
#print(f"beta2: {beta2}")

k = 1.96
n = len(df)
corr_fact = 1 + 1/n
MSE = mean_squared_error(y, y_pred)
#print("MSE:", MSE)
mono_mean = np.mean(X)
diff_mass_sum = np.sum((X - mono_mean)**2)

print(k, MSE, corr_fact, diff_mass_sum, mono_mean, beta0, beta1)'''


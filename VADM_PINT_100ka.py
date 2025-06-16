# **Rolling window - 100.000 years**

This study focuses on analyzing paleointensity data from the PINT database, specifically the *pint10* subset spanning the past 10 million years. The primary objective is to explore the evolution of the Earth's magnetic field intensity through time using statistical methods, including cubic spline interpolation for estimating missing values, and Monte Carlo simulations and bootstrap approaches for uncertainty analysis.

import pandas as pd
import matplotlib.pyplot as plt
import tabulate
import datetime
import numpy as np
import scipy as sci
import scipy.interpolate as si
import seaborn as sns
import altair as alt
import gc
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import CubicSpline
from datetime import datetime
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.base import TransformerMixin
from sklearn.pipeline import make_pipeline
import os
from patsy import dmatrix
import statsmodels.api as sm

import warnings
warnings.filterwarnings('ignore')


## Data Preparation


# Función para crear directorios
def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directorio creado: {directory_path}")

# Montar Google Drive
from google.colab import drive
drive.mount('/content/drive/', force_remount=True)

# Directorio base para resultados
results_base_dir = '/content/drive/MyDrive/Results'

# Crear directorio único por run con timestamp
current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
run_dir = os.path.join(results_base_dir, f"run_{current_time}")
create_directory(run_dir)

# Subdirectorios
plots_dir = os.path.join(run_dir, "plots")
dataframes_dir = os.path.join(run_dir, "dataframes")
results_dir = os.path.join(run_dir, "results")

# Crear subdirectorios
create_directory(plots_dir)
create_directory(dataframes_dir)
create_directory(results_dir)

# 1. Guardar DataFrames
def save_dataframe(df, filename, directory=dataframes_dir):
    csv_path = os.path.join(directory, f"{filename}.csv")

    df.to_csv(csv_path, index=False)
    print(f"DataFrame guardado en: {csv_path}")

# 2. Guardar gráficos
def save_plot(fig, filename, directory=plots_dir):
  plot_path = os.path.join(directory, f"{filename}.png")
  plt.savefig(plot_path, dpi=300)
  print(f"Gráfico guardado en: {plot_path}")
  plt.close(fig)

# 3. Guardar resultados intermedios
def save_results(data, filename, directory=results_dir):
    results_path = os.path.join(directory, f"{filename}.csv")
    np.savetxt(results_path, data, delimiter=",")
    print(f"Resultados guardados en: {results_path}")

Selected the PINT10 subset from the PINT database, including data entries up to 10 million years old.

from google.colab import drive
drive.mount('/content/drive/',force_remount=True),

pint = pd.read_excel('/content/drive/MyDrive/Colab Notebooks/PINTv811.xlsx')

# ----------------------- FILTER BY SAMPLE'S TYPE ----------------------
igneous_groups = ['Volcanic', 'volcanic', 'Plutonic', 'Intrusive']
pint_10my = pint_10my[pint_10my['GROUP'].isin(igneous_groups)]
# ----------------------- FILTER BY SAMPLE'S AGE ----------------------
filtro_edad = 5
pint_10my = pint[pint['AGE'] <= filtro_edad]
pint_10my

group_counts = pint_10my.groupby('GROUP').size()
sorted_groups = group_counts.sort_values(ascending=False)

top_10_groups = sorted_groups.head(10)

# Imprime los resultados
print(top_10_groups.to_markdown(numalign="left", stralign="left"))

columns_to_check = ['VDM', 'VADM', 'VDM/VADM']
df_all_null_selected_cols = pint_10my[pint_10my[columns_to_check].isnull().all(axis=1)]

print(df_all_null_selected_cols[columns_to_check].isnull().all(axis=1).sum())

print(df_all_null_selected_cols.head(20).to_markdown(index=False, numalign="left", stralign="left"))


Focused on relevant columns, including the Virtual Dipole Moment (VDM), Virtual Axial Dipole Moment (VADM), and age of the data point.
Sorted the dataset according to age in ascending order.

columns_vadm = ['AGE', 'VADM', 'VDM', 'VDM/VADM']
data_vadm = pint_10my[columns_vadm]
data_vadm.sort_values(by='AGE', ascending=True, inplace=True)
data_vadm

MCADAM dataset

### Load New Data ###
new_data_path = '/content/drive/MyDrive/Colab Notebooks/MCADAMv1a.xlsx'  # Replace with the actual path
try:
    new_data = pd.read_excel(new_data_path, index_col =0)
    # Asegúrate de que el nuevo dataframe tenga AGE y VADM
    new_data = new_data[['age', 'mean']]  # Ajusta los nombres de las columnas si es necesario
    new_data = new_data.rename(columns={'age': 'AGE', 'mean': 'VADM'}) #Rename for standarization
    # Filtra los datos para incluir solo edades menores o iguales a 10 My
    new_data = new_data[new_data['AGE'] <= 10]

except FileNotFoundError:
    print(f"Error: File not found at {new_data_path}")
    new_data = pd.DataFrame({'AGE': [], 'VADM': []})#Dataframe  vacio


new_data.head()

### Chrons definition for Polarity Reversals


VADM_actual = 8.22

chrones = {
    "C1 (Brunhes)": {"base": 0.000, "techo": 0.781},
    "C2 (Matuyama)": {"base": 0.781, "techo": 2.581},
    "C3 (Gauss)": {"base": 2.581, "techo": 3.596},
    "C4 (Gilbert)": {"base": 3.596, "techo": 7.537},
    "C5 ": {"base": 7.537, "techo": 9.786},
}
# Crear la figura y el eje
fig, ax = plt.subplots(figsize=(16, 2))

# Colores en tonos de negro y blanco
colores = ["black", "white", "black", "white", "black"]

# Graficar cada cron como una barra horizontal
for i, (cron, datos) in enumerate(chrones.items()):
    ax.barh(0, width=datos["techo"] - datos["base"], left=datos["base"], height=0.6, color=colores[i], edgecolor="black", label=cron)

# Ocultar el eje Y
ax.yaxis.set_visible(False)

# Ocultar el marco del gráfico
for spine in ax.spines.values():
    spine.set_visible(False)

# Añadir etiquetas de los crones
for cron, datos in chrones.items():
    ax.text((datos["base"] + datos["techo"]) / 2, 0, cron, ha="center", va="center", color="white" if colores[list(chrones.keys()).index(cron)] == "black" else "black", fontsize=8)

# Añadir título
ax.set_title("Chrons Temporal Scale (GPTS, Ogg 2020)", pad=20)

# Ajustar límites del eje X
ax.set_xlim(0, 10)

# Mostrar la gráfica
plt.tight_layout()
plt.show()

## **Chrones completos**

# Subchrones de polaridad magnética
subchrones = {
    "C1n (Brunhes)": {"techo": 0.000, "base": 0.773},
    "C1r.1r (Matuyama)": {"techo": 0.773, "base": 0.990},
    "C1r.1n (Jaramillo)": {"techo": 0.990, "base": 1.070},
    "C1r.2r": {"techo": 1.070, "base": 1.180},
    "C1r.2n (Cobb Mountain)": {"techo": 1.180, "base": 1.215},
    "C1r.3r": {"techo": 1.215, "base": 1.775},

    "C2n (Olduvai)": {"techo": 1.775, "base": 1.934},

    "C2r.1r": {"techo": 1.934, "base": 2.116},
    "C2r.1n (Feni)": {"techo": 2.116, "base": 2.140},
    "C2r.2r (Matuyama)": {"techo": 2.140, "base": 2.595},
    "C2An.1n (Gauss)": {"techo": 2.595, "base": 3.032},
    "C2An.1r (Keana)": {"techo": 3.032, "base": 3.116},
    "C2An.2n": {"techo": 3.116, "base": 3.207},
    "C2An.2r (Mammoth)": {"techo": 3.207, "base": 3.330},
    "C2An.3n (Gauss)": {"techo": 3.330, "base": 3.596},
    "C2Ar (Gilbert)": {"techo": 3.596, "base": 4.187},

    "C3n.1n (Cochiti)": {"techo": 4.187, "base": 4.300},
    "C3n.1r": {"techo": 4.300, "base": 4.631},
    "C3n.2n (Nunivak)": {"techo": 4.631, "base": 4.799},
    "C3n.2r": {"techo": 4.799, "base": 4.896},
    "C3n.3n (Sidufjall)": {"techo": 4.896, "base": 4.997}
}

# Crear la figura y el eje
fig, ax = plt.subplots(figsize=(16, 2.6))

# Colores en tonos de negro y blanco
colores = ["black", "white"] * (len(subchrones) // 2 + 1)

# Graficar cada subcron como una barra horizontal
for i, (subcron, datos) in enumerate(subchrones.items()):
    ax.barh(0, width=datos["techo"] - datos["base"], left=datos["base"], height=0.6, color=colores[i], edgecolor="black")

# Ocultar el eje Y
ax.yaxis.set_visible(False)

# Ocultar el marco del gráfico
for spine in ax.spines.values():
    spine.set_visible(False)

# Añadir etiquetas de los subchrones (opcional, puede ser demasiado denso)
for subcron, datos in subchrones.items():
    ax.text((datos["base"] + datos["techo"]) / 2, 0, subcron, ha="center", va="center", color="white" if colores[list(subchrones.keys()).index(subcron)] == "black" else "black", fontsize=7, rotation=90)

# Añadir título
ax.set_title("Subchrons Temporal Scale (GPTS, Ogg 2020)", pad=20)

# Ajustar límites del eje X
ax.set_xlim(0, 5)

# Mostrar la gráfica
plt.tight_layout()
plt.show()

### Treatment for null data at VADM, VDM and VDM/VADM


data_vadm.loc[data_vadm['VADM'].isnull() & data_vadm['VDM'].notnull(), 'VADM'] = data_vadm['VDM']
data_vadm.loc[data_vadm['VADM'].isnull() & data_vadm['VDM'].isnull() & data_vadm['VDM/VADM'].notnull(), 'VADM'] = data_vadm['VDM/VADM']
columns= ['AGE', 'VADM']
data_vadm = data_vadm[columns]
data_vadm.to_csv('age_vs_vadm', index=False)
# Save to result folder
save_dataframe(data_vadm, "age_vs_vadm")
np.random.seed(42)

#data_vadm.to_csv('test/age_vs_vadm.csv', index=False)
#data_vadm.to_excel('test/age_vs_vadm.xlsx', index=False)

data_vadm.to_csv('age_vs_vadm', index=False)


### Rolling Window

It was applied a rolling window to segment the dataset in 0.005 intervals between 0 and 10 million years. Also
it was calculated the mean and standard deviation of VADM for each interval.

# 1. Create the 'window' column
data = data_vadm.copy()
sequence = np.arange(0, 5.01, 0.05)
new_column = pd.Series(np.nan, index=data.index)

assignment_limit = min(len(sequence), len(new_column))
new_column[:assignment_limit] = sequence[:assignment_limit]

data['window'] = new_column

data

# 2. Calculate mean VADM for each window

def calculate_mean_vadm(ventana, df, age_column):
    filtered_df = df[(df[age_column] >= 0) & (df[age_column] < ventana)]
    if filtered_df.empty:
        return np.nan
    else:
        return data_vadm['VADM'].mean()

data['mean_VADM'] = data['window'].apply(lambda x: calculate_mean_vadm(x, data, 'AGE'))

# Refine mean_VADM calculation within each window
unique_ventana_values = sorted(data['window'].unique())
for i, ventana_value in enumerate(unique_ventana_values):
    if ventana_value > 5.01:
        data.loc[data['window'] == ventana_value, 'mean_VADM'] = np.nan
    else:
        lower_limit = 0 if i == 0 else unique_ventana_values[i - 1]
        data.loc[data['window'] == ventana_value, 'mean_VADM'] = data[
            (data['AGE'] >= lower_limit) &
            (data['AGE'] <= ventana_value)
        ]['VADM'].mean()


data.head(100)

# 3. Calculate standard deviation of VADM for each window

def calculate_std_vadm(ventana, df, age_column):
    data = df[(df[age_column] >= 0) & (df[age_column] < ventana)]
    if data.empty or len(data) == 1:
        return np.nan
    else:
        return data['VADM'].std()

unique_ventana_values = sorted(data['window'].unique())

for i, ventana_value in enumerate(unique_ventana_values):
    if ventana_value > 5.01:
        data_vadm.loc[data['window'] == ventana_value, 'std_VADM'] = np.nan
    else:
        lower_limit = 0 if i == 0 else unique_ventana_values[i - 1]
        data.loc[data['window'] == ventana_value, 'std_VADM'] = data[
            (data['AGE'] >= lower_limit) &
            (data['AGE'] <= ventana_value)
        ]['VADM'].std(ddof=0)

# Save to result folder
save_dataframe(data, "RW_vs_vadm")

df_plot = data.dropna(subset=['mean_VADM'])

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 7), gridspec_kw={'height_ratios': [4, 1]}, sharex=True)
ax1.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8,linewidth=0.4)
ax1.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.6)

ax1.errorbar(df_plot['window'], df_plot['mean_VADM'], yerr=df_plot['std_VADM'], fmt='none', ecolor='gray', elinewidth=0.3, capsize=3)
ax1.errorbar(df_plot['window'], df_plot['mean_VADM'], xerr=0.05, fmt='none', ecolor='gray', elinewidth=0.3, capsize=3)
ax1.plot(df_plot['window'], df_plot['mean_VADM'], color='blue', linestyle='-', linewidth=1.2, label='Mean VADM 100 ka RW')
ax1.axhline(y=VADM_actual, color='blueviolet', linestyle='--', linewidth=0.5, label='VADM today (8.22 × 10²² Am²)')

ax1.set_title("Mean VADM with Rolling Window (100 ka) for last 5 My ", fontsize=15, pad=20)
ax1.set_ylabel('VADM (x10^22 Am²)', fontsize=12)
ax1.set_xticks(np.arange(0, 5, 1))
ax1.grid(True, linestyle='--', linewidth=0.2)
ax1.legend(loc='upper right', fontsize=12)

ax1.set_xlim(0, 5)
ax1.set_ylim(0, 20)

ax1.fill_between(df_plot['window'], df_plot['mean_VADM'] - df_plot['std_VADM'], df_plot['mean_VADM'] + df_plot['std_VADM'], color='lightblue', alpha=0.3, label='RW ± 1σ')

colores = ["black", "white"] * (len(subchrones) // 2 + 1)

for i, (subcron, datos) in enumerate(subchrones.items()):
    ax2.barh(0, width=datos["techo"] - datos["base"], left=datos["base"], height=0.6, color=colores[i], edgecolor="black")

ax2.yaxis.set_visible(False)
ax2.set_xlabel('Age (My)', fontsize=12)

for spine in ax2.spines.values():
    spine.set_visible(False)

for subcron, datos in subchrones.items():
    ax2.text((datos["base"] + datos["techo"]) / 2, 0, subcron, ha="center", va="center", color="white" if colores[list(subchrones.keys()).index(subcron)] == "black" else "black", fontsize=7, rotation=90)

ax2.set_xlim(0, 5)

plt.subplots_adjust(hspace=0.08)
plt.tight_layout()
save_plot(fig, filename="RW_100ka_5myr_meanVADM")
plt.close()

## Data Preprocessing

Addressed potential issues caused by missing data in the original dataset using cubic spline interpolation. This method leverages the known age and VADM values to estimate missing points in the timeseries, establishing a smooth curve that connects the data points.

### Cubic spline

df_filtered = data.dropna(subset=['mean_VADM', 'window'])

x = df_filtered['window'].values
y = df_filtered['mean_VADM'].values

cs = CubicSpline(x, y)

x_new = np.linspace(x.min(), x.max(), 100)
y_new = cs(x_new)

# Creación del DataFrame cubicspline
cubic_spline_df = pd.DataFrame({'window': x_new, 'mean_VADM': y_new})
save_dataframe(cubic_spline_df, filename="cubic_spline_RW")

plt.figure(figsize=(20, 7))
plt.plot(x_new, y_new, '-', label='Cubic spline interpolation', color='darkorange')
plt.scatter(x,y,label='PINT', facecolors='none', edgecolors='black', s=8,linewidth=0.4)
plt.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.6)
plt.axhline(y=VADM_actual, color='blueviolet', linestyle='--', linewidth=0.4, label='VADM today (8.22 × 10²² Am²)')

plt.title(' Cubic interpolation (mean_VADM vs window)')
plt.xlabel('window')
plt.ylabel('mean_VADM')
ax = plt.gca()
ax.set_xlim(0, 5)
ax.set_ylim(0, 20)


plt.grid(True, linestyle='--',linewidth=0.2)
plt.legend(loc='upper left', fontsize=10)

plt.tight_layout()
fig = plt.gcf()  # Obtener la figura actual para guardarla

#save_plot(fig, filename = "cubicspline")
#plt.savefig('/content/drive/MyDrive/Results/cubicspline.png', dpi=300)
plt.show()

# --- Load Data (Ensure paths are correct) ---
data_vadm = pd.read_csv('age_vs_vadm')
#df_plot = pd.read_csv('RW_vs_vadm.csv').dropna(subset=['mean_VADM'])

# --- Parameters ---
smoothing_factors = [1]  # Range of smoothing factors to explore

# --- Prepare Data ---
x = df_plot['window'].values
y = df_plot['mean_VADM'].values
x_eval = np.linspace(0, 5, 200)  # Evaluate over the full range

# --- Metrics calculation - Evaluate agains a data known -----
x_true = data_vadm['AGE'].values
y_true = data_vadm['VADM'].values

# --- remove nan values for x_true
nan_indices = np.isnan(x_true) | np.isnan(y_true)
x_true = x_true[~nan_indices]
y_true = y_true[~nan_indices]

# --- Store Results ---
spline_results = {}  # Store results for each smoothing factor

# --- Iterate over Smoothing Factors ---
for smoothing_factor in smoothing_factors:
    print(f"Evaluating smoothing factor: {smoothing_factor}")

    try:
        # --- Fit Smoothed Spline ---
        us = UnivariateSpline(x, y, s=smoothing_factor)  # Use UnivariateSpline with smoothing

        # --- Calculate MSE ---
        y_pred = us(x_eval)

        #-- Calculate MSE
        indices = np.where((x_true >= x_eval.min()) & (x_true <= x_eval.max()))

        mse = mean_squared_error(y_true[indices], us(x_true[indices]))
        # Ensure same len

        # --- Store Results ---
        spline_results[smoothing_factor] = {
            'spline': us,
            'mse': mse
        }
    except Exception as e:
        print(f"Error with smoothing factor {smoothing_factor}: {e}")

# --- Plot the Results ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 7), sharex=True, gridspec_kw={'height_ratios': [4, 1]})

# --- Plot Original Data (Subplot 1) ---
ax1.scatter(x_true, y_true, label='PINT Data',facecolors='none', edgecolors='black', s=8,linewidth=0.4)
# --- Plot MCADAM Data ---
ax1.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
ax1.axhline(y=VADM_actual, color='blueviolet', linestyle='--', linewidth=0.6, label='VADM today (8.22 × 10²² Am²)')

# --- Find the best model
best_smoothing = min(spline_results, key=lambda k: spline_results[k]['mse'])
print(f"\n  the best model selected. With low MSE: {best_smoothing}")

# --- Plot with each model
for smoothing_factor, spline_data in spline_results.items():
    y_pred = spline_data['spline'](x_eval)

    if best_smoothing == smoothing_factor:

      # Plot the derivative
      ax1.plot(x_eval, y_pred, label=f"Spline (s={smoothing_factor}, MSE={spline_data['mse']:.2f}) <-- SELECTED " , linewidth = 2, color = 'darkorange')
    else:
        # Plot the derivative
      ax1.plot(x_eval, y_pred, label=f"Spline (s={smoothing_factor}, MSE={spline_data['mse']:.2f})")


ax1.set_title("Cubic Spline with Different Smoothing Factors")
ax1.set_ylabel("VADM")
ax1.legend()
ax1.grid(True, linestyle='--', linewidth=0.2)
ax1.set_xlim(0, 5)
ax1.set_ylim(0, 20)


# --- Plot Subchrones (Subplot 2) ---
colores = ["black", "white"] * (len(subchrones) // 2 + 1)
for i, (subcron, datos) in enumerate(subchrones.items()):
    ax2.barh(0, width=datos["techo"] - datos["base"], left=datos["base"], height=0.6, color=colores[i % len(colores)], edgecolor="black")

ax2.yaxis.set_visible(False)
ax2.set_xlabel('Age (Ma)', fontsize=12)
for spine in ax2.spines.values():
    spine.set_visible(False)
for subcron, datos in subchrones.items():
    ax2.text((datos["base"] + datos["techo"]) / 2, 0, subcron, ha="center", va="center", color="white" if colores[list(subchrones.keys()).index(subcron) % len(colores)] == "black" else "black", fontsize=7, rotation=90)
ax2.set_xlim(0, 5)


# Ajustar el espaciado entre subplots
plt.subplots_adjust(hspace=0.08)
plt.tight_layout()
save_plot(fig, filename = "cubicspline_with_subchrons")
plt.close()

# --- Print Metrics Summary ---
print("\n--- Metrics Summary ---")
for smoothing_factor, spline_data in spline_results.items():
    print(f"Smoothing Factor: {smoothing_factor}, MSE: {spline_data['mse']:.4f}")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 7),
                              gridspec_kw={'height_ratios': [4, 1]},
                              sharex=True)


ax1.plot(x_new, y_new, '-', label='Cubic spline interpolation', color='darkorange', linewidth=1.5)
ax1.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8,linewidth=0.4)
ax1.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
ax1.axhline(y=VADM_actual, color='blueviolet', linestyle='--', linewidth=0.4, label='VADM today (8.22 × 10²² Am²)')

ax1.set_title("Cubic Spline Interpolation of VADM for last 5 My)",
              fontsize=15, pad=20)
ax1.set_ylabel('VADM (x10^22 Am²)', fontsize=12)
ax1.grid(True, linestyle='--', linewidth=0.2)

#ax1.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
ax1.legend(loc='upper right', fontsize=10)
ax1.set_xlim(0, 5)
ax1.set_ylim(0, 20)

colores = ["black", "white"] * (len(subchrones) // 2  )

for i, (subcron, datos) in enumerate(subchrones.items()):
    ax2.barh(0, width=datos["techo"] - datos["base"],
             left=datos["base"], height=0.6,
             color=colores[i % len(colores)],
             edgecolor="black")

ax2.yaxis.set_visible(False)
ax2.set_xlabel('Age (Myr)', fontsize=12)
for spine in ax2.spines.values():
    spine.set_visible(False)

for subcron, datos in subchrones.items():
    ax2.text((datos["base"] + datos["techo"]) / 2, 0, subcron,
             ha="center", va="center",
             color="white" if colores[list(subchrones.keys()).index(subcron) % len(colores)] == "black" else "black",
             fontsize=7, rotation=90)

ax2.set_xlim(0, 5)
plt.tight_layout()
plt.subplots_adjust(hspace=0.08)
save_plot(fig, filename = "cubicspline_with_subchrones")
plt.show()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 7), gridspec_kw={'height_ratios': [4, 1]}, sharex=True)

# --- Plot del VADM, Rolling Window y Cubic Spline (en el primer subplot) ---
# 1. Puntos de datos originales
ax1.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8,linewidth=0.4)
ax1.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)

# 2. Curva de Rolling Window ( VADM promedio)
ax1.plot(df_plot['window'], df_plot['mean_VADM'], color='blue', linestyle='-', linewidth=1.1, label='Window 100 ka')
ax1.axhline(y=VADM_actual, color='blueviolet', linestyle='--', linewidth=0.6, label='VADM today (8.22 × 10²² Am²)')

# Barras de error del rolling window
#ax1.errorbar(df_plot['window'], df_plot['mean_VADM'], yerr=df_plot['std_VADM'], fmt='none', ecolor='lightcoral', elinewidth=0.5, capsize=2, label='_nolegend_')

# 3. Curva de Cubic Spline
ax1.plot(x_new, y_new, color='orange', linestyle='-', linewidth=1, label='Cubic Spline ')

# --- Personalización del primer subplot ---
ax1.set_ylabel('VADM (x10^22 Am²)', fontsize=12)
ax1.set_title('VADM vs. AGE fitted with original data, (100 ka) Rolling Window and Cubic Spline imputation', fontsize=14)
ax1.legend(fontsize=10, loc='upper right')
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.set_xlim(0, 5)
ax1.set_ylim(0, 20)


# --- Plot de la escala de subchrones (en el segundo subplot) ---
colores = ["black", "white"] * (len(subchrones) // 2 + 1)

for i, (subcron, datos) in enumerate(subchrones.items()):
    ax2.barh(0, width=datos["techo"] - datos["base"], left=datos["base"], height=0.6, color=colores[i % len(colores)], edgecolor="black")

# --- Personalización del segundo subplot ---
ax2.yaxis.set_visible(False)
ax2.set_xlabel('AGE (Ma)', fontsize=12)
for spine in ax2.spines.values():
    spine.set_visible(False)

for subcron, datos in subchrones.items():
    ax2.text((datos["base"] + datos["techo"]) / 2, 0, subcron, ha="center", va="center", color="white" if colores[list(subchrones.keys()).index(subcron) % len(colores)] == "black" else "black", fontsize=7, rotation=90)

ax2.set_xlim(0, 5)
plt.subplots_adjust(hspace=0.08)
plt.tight_layout()
# --- Ajustar el espaciado entre subplots ---
save_plot(fig, filename = "cubicspline_with_rw")
plt.close()


# --- Crear la gráfica ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 7), gridspec_kw={'height_ratios': [4, 1]}, sharex=True)

# --- Plot del VADM, Rolling Window y Cubic Spline (en el primer subplot) ---
# 1. Puntos de datos originales
ax1.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8,linewidth=0.4)
ax1.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)

# 2. Curva de Rolling Window (VADM promedio)
ax1.plot(df_plot['window'], df_plot['mean_VADM'], color='blue', linestyle='-', linewidth=0.8, label='Rolling Window 100 ka')

# 3. Área sombreada del Rolling Window (mean +/- std)
ax1.fill_between(df_plot['window'], df_plot['mean_VADM'] - df_plot['std_VADM'], df_plot['mean_VADM'] + df_plot['std_VADM'], color='lightblue', alpha=0.3, label='mean_VADM ± σ')
ax1.axhline(y=VADM_actual, color='blueviolet', linestyle='--', linewidth=0.6, label='VADM today (8.22 × 10²² Am²)')

# 4. Curva de Cubic Spline
ax1.plot(x_new, y_new, color='darkorange', linestyle='--', linewidth=1.4, label='Cubic Spline')

# --- Personalización del primer subplot ---
ax1.set_ylabel('VADM (x10^22 Am²)', fontsize=12)
ax1.set_title('Mean VADM vs. AGE with original data, Rolling Window (100 ka) and Cubic Spline imputation for last 5 My', fontsize=14)
ax1.legend(fontsize=10, loc='upper right')
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.set_xlim(0, 5)
ax1.set_ylim(0, 20)


# --- Plot de la escala de subchrones (en el segundo subplot) ---
colores = ["black", "white"] * (len(subchrones) // 2 + 1)

for i, (subcron, datos) in enumerate(subchrones.items()):
    ax2.barh(0, width=datos["techo"] - datos["base"], left=datos["base"], height=0.6, color=colores[i % len(colores)], edgecolor="black")

# --- Personalización del segundo subplot ---
ax2.yaxis.set_visible(False)
ax2.set_xlabel('AGE (My)', fontsize=12)
for spine in ax2.spines.values():
    spine.set_visible(False)

for subcron, datos in subchrones.items():
    ax2.text((datos["base"] + datos["techo"]) / 2, 0, subcron, ha="center", va="center", color="white" if colores[list(subchrones.keys()).index(subcron) % len(colores)] == "black" else "black", fontsize=7, rotation=90)

ax2.set_xlim(0, 5)

# --- Ajustar el espaciado entre subplots ---
plt.tight_layout()
plt.subplots_adjust(hspace=0.1)

# --- Guardar y mostrar la gráfica ---
#plt.savefig('/content/drive/MyDrive/Results/rolling_window_cubicspline_shaded.png', dpi=300, bbox_inches='tight')

save_plot(fig, filename = "cubicspline_with_subchroness")
plt.show()

# Bootstrap


# --- Parameter Ranges to Explore ---
spline_smoothings = [0]  # Example smoothing values for CubicSpline
n_bootstrap_iterations = [1998, 1999, 2000, 2001, 2002]  # Example number of bootstrap iterations

# --- Store Results ---
spline_results = {}
bootstrap_results = {}
mse_results = []

# --- Define Metric for Evaluation (Example: MSE) ---
def evaluate_spline(x_true, y_true, x_eval, spline):
    """Calculates MSE for a given spline."""
    y_pred = spline(x_eval)

    #HACK para que cuando los datos no se puedan evaluar, no se caiga la simulación
    #y le asigno un valor de perdida grande
    if(len(y_pred) != len(y_true)):
        return 999999999

    mse = mean_squared_error(y_true, y_pred)
    return mse

def resample_with_replacement(x, y, seed=None):
    """Function to resample data with replacement."""
    if seed is not None:
        np.random.seed(seed)

    indices_resample = np.random.choice(len(x), len(x), replace=True)
    x_sample = x[indices_resample]
    y_sample = y[indices_resample]
    sort_indices = np.argsort(x_sample)
    x_sample = x_sample[sort_indices]
    y_sample = y_sample[sort_indices]

    x_unique, indices = np.unique(x_sample, return_index=True)
    x_sample = x_sample[indices]
    y_sample = y_sample[indices]

    return x_sample, y_sample

# --- Spline Parameter Iteration ---
for smoothing in spline_smoothings:
    x = df_plot['window'].values
    y = df_plot['mean_VADM'].values
    try:
        cs = CubicSpline(x, y, extrapolate=True)  # Adjust smoothing if necessary
        x_eval = np.linspace(x.min(), x.max(), 200)
        y_true = df_plot['mean_VADM'].values
        mse = evaluate_spline(y_true, cs(x_eval)[:len(y_true)], x_eval, cs)

        spline_results[smoothing] = {
            'spline': cs,
            'mse': mse,
            'x_eval': x_eval,
        }
    except Exception as e:
        print(f"CubicSpline failed with smoothing={smoothing}: {e}")

# --- Bootstrap Iteration ---
for n_bootstraps in n_bootstrap_iterations:
    print(f"Evaluating Bootstrap with n_bootstraps={n_bootstraps}")
    spline_curves = []
    SIMULATION_OK = True

    for i in range(n_bootstraps):
        max_attempts = 3
        for attempt in range(max_attempts):
            try:
                x_sample, y_sample = resample_with_replacement(df_plot['window'].values, df_plot['mean_VADM'].values)
                cs = CubicSpline(x_sample, y_sample, extrapolate=True)
                spline_curves.append(cs)
                break
            except Exception as e:
                print(f"Attempt {attempt+1} failed on iteration {i}: {e}")
                SIMULATION_OK = attempt < max_attempts - 1
                if not SIMULATION_OK:
                    break

    if SIMULATION_OK and spline_curves:
        x_new = np.linspace(x.min(), x.max(), 200)
        y_values = np.array([cs(x_new) for cs in spline_curves])
        y_lower = np.percentile(y_values, 2.5, axis=0)
        y_upper = np.percentile(y_values, 97.5, axis=0)
        y_mean_bootstrap = np.mean(y_values, axis=0)

        y_true = df_plot['mean_VADM'].values
        mse = mean_squared_error(y_true, y_mean_bootstrap[:len(y_true)])

        bootstrap_results[n_bootstraps] = {
            'x_new': x_new,
            'y_mean_bootstrap': y_mean_bootstrap,
            'y_lower': y_lower,
            'y_upper': y_upper,
            'mse': mse
        }

        bootstrap_df = pd.DataFrame({
            'x_new': x_new,
            'y_mean_bootstrap': y_mean_bootstrap,
            'y_lower': y_lower,
            'y_upper': y_upper
        })
        save_dataframe(bootstrap_df, f'bootstrap_{n_bootstraps}', results_dir)

        mse_results.append({'n_bootstraps': n_bootstraps, 'mse': mse})

# --- Save Results ---
mse_df = pd.DataFrame(mse_results)
save_dataframe(mse_df, 'bootstrap_mse_summary', results_dir)


# --- Visualization ---
fig, axes = plt.subplots(len(spline_smoothings) + len(n_bootstrap_iterations), 1, figsize=(20, 5 * (len(spline_smoothings) + len(n_bootstrap_iterations))), sharex=True)


# Plot Spline Results
for i, (smoothing, result) in enumerate(spline_results.items()):
    ax = axes[i]
    ax.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
    ax.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)

    ax.plot(result['x_eval'], result['spline'](result['x_eval']), label=f'CubicSpline (smoothing={smoothing}, MSE={result["mse"]:.2f})', color='darkorange')
    ax.legend()

    ax.set_xlim(0, 5)
    ax.set_ylim(0, 20)
    ax.set_title(f'CubicSpline with smoothing={smoothing}')
    ax.set_ylabel('VADM')

# Plot Bootstrap Results
for i, (n_bootstraps, result) in enumerate(bootstrap_results.items()):
    ax = axes[len(spline_smoothings) + i]  # Offset index
    ax.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
    ax.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
    ax.axhline(y=VADM_actual, color='blueviolet', linestyle='--', linewidth=0.6, label='VADM today (8.22 × 10²² Am²)')

    ax.plot(result['x_new'], result['y_mean_bootstrap'], label=f'Bootstrap (n={n_bootstraps}, MSE={result["mse"]:.2f})', color='forestgreen')
    ax.fill_between(result['x_new'], result['y_lower'], result['y_upper'], color='lightgreen', alpha=0.4, label='95% CI')
    ax.legend()

    ax.set_xlim(0, 5)
    ax.set_ylim(0, 20)
    ax.set_title(f'Bootstrap with {n_bootstraps} iterations')
    ax.set_ylabel('VADM')
    ax.grid(True, linestyle='--', alpha=0.6)



plt.xlabel("Age (Ma)")
plt.tight_layout()
save_plot(fig, filename = "bootstrap_fit")
plt.show()

# --- Print Summary of Metrics ---

print("\n--- Summary of Bootstrap Results ---")
for n_bootstraps, result in bootstrap_results.items():
    print(f"n_bootstraps={n_bootstraps}: MSE={result['mse']:.4f}")

#### Fitting best n for Bootstrap

# Función para obtener input de usuario sobre el rango
def get_bootstrap_range():
    try:
        min_n_bootstrap = int(input("Ingrese el número mínimo de iteraciones de bootstrap: "))
        max_n_bootstrap = int(input("Ingrese el número máximo de iteraciones de bootstrap: "))
        step_n_bootstrap = int(input("Ingrese el incremento del rango de bootstrap: "))

        if min_n_bootstrap <= 0 or max_n_bootstrap <= 0 or step_n_bootstrap <= 0:
            raise ValueError("Todos los valores deben ser positivos y mayores que cero.")
        return range(min_n_bootstrap, max_n_bootstrap + 1, step_n_bootstrap)
    except ValueError as e:
        print(f"Entrada inválida: {e}")
        exit()

n_bootstrap_iterations = get_bootstrap_range()

# --- Store Results ---
bootstrap_results = {}

# --- Bootstrap Iteration ---
for n_bootstraps in n_bootstrap_iterations:
    print(f"Evaluating Bootstrap with n_bootstraps={n_bootstraps}")
    y_values = []

    for i in range(n_bootstraps):
        try:
            x_sample, y_sample = resample_with_replacement(df_plot['window'].values, df_plot['mean_VADM'].values)
            cs = CubicSpline(x_sample, y_sample, extrapolate=True)
            x_new = np.linspace(x.min(), x.max(), 200)
            y_values.append(cs(x_new))
        except Exception as e:
            print(f"Error during bootstrap iteration {i}: {e}")
            continue

    if y_values:
        y_values = np.array(y_values)
        y_mean_bootstrap = np.mean(y_values, axis=0)
        mse = mean_squared_error(df_plot['mean_VADM'].values, y_mean_bootstrap[:len(df_plot['mean_VADM'].values)])

        # Store Bootstrap Results
        bootstrap_results[n_bootstraps] = {
            'y_values': y_values,
            'x_new': x_new
        }
        mse_results.append({'n_bootstraps': n_bootstraps, 'mse': mse})

        # Save each bootstrap result
        bootstrap_df = pd.DataFrame({
            'x_new': x_new,
            'y_mean_bootstrap': y_mean_bootstrap
        })
        save_dataframe(bootstrap_df, f'bootstrap_{n_bootstraps}', results_dir)

    else:
        print(f"No valid splines could be generated for n_bootstraps={n_bootstraps}.")

# --- Save MSE Results ---
mse_df = pd.DataFrame(mse_results)
save_dataframe(mse_df, 'bootstrap_mse_summary', results_dir)

# --- Visualization ---
fig, ax = plt.subplots(figsize=(20, 7))

# Plot Bootstrap Curves
for n_bootstraps, result in bootstrap_results.items():
    for y_curve in result['y_values']:
        ax.plot(result['x_new'], y_curve, alpha=0.4, linewidth=0.8)

# Plot data points
ax.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
ax.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
ax.axhline(y=VADM_actual, color='blueviolet', linestyle='-.', linewidth=0.4, label='VADM today (8.22 × 10²² Am²)')

ax.set_xlim(0, 5)
ax.set_ylim(0, 20)
ax.set_title('Bootstrap Curves')
ax.set_ylabel('VADM')
ax.set_xlabel('Age (Ma)')
ax.legend()

plt.tight_layout()
#plt.savefig('/content/drive/MyDrive/Results/bootstrap_curves.png', dpi=300)
save_plot(fig, 'bootstrap_curves', plots_dir)

plt.show()

# --- Print Summary of Metrics ---
#print("\n--- Summary of Bootstrap Results ---")
#for n_bootstraps, result in bootstrap_results.items():
#    print(f"n_bootstraps={n_bootstraps}: MSE={result['mse']:.4f}")


# --- Store Results ---
bootstrap_results = {}
mse_results = []

# --- Bootstrap Iteration ---
for n_bootstraps in n_bootstrap_iterations:
    print(f"Evaluating Bootstrap with n_bootstraps={n_bootstraps}")
    y_values = []

    for i in range(n_bootstraps):
        try:
            x_sample, y_sample = resample_with_replacement(df_plot['window'].values, df_plot['mean_VADM'].values)
            cs = CubicSpline(x_sample, y_sample, extrapolate=True)
            x_new = np.linspace(x.min(), x.max(), 200)
            y_values.append(cs(x_new))
        except Exception as e:
            print(f"Error during bootstrap iteration {i}: {e}")
            continue

    if y_values:
        y_values = np.array(y_values)
        y_mean_bootstrap = np.mean(y_values, axis=0)
        y_lower = np.percentile(y_values, 2.5, axis=0)  # 2.5% percentile
        y_upper = np.percentile(y_values, 97.5, axis=0)  # 97.5% percentile

        mse = mean_squared_error(df_plot['mean_VADM'].values, y_mean_bootstrap[:len(df_plot['mean_VADM'].values)])

        # Store Bootstrap Results
        bootstrap_results[n_bootstraps] = {
            'y_values': y_values,
            'y_mean': y_mean_bootstrap,
            'y_lower': y_lower,
            'y_upper': y_upper,
            'x_new': x_new,
            'mse': mse
        }
        mse_results.append({'n_bootstraps': n_bootstraps, 'mse': mse})

        # Save each bootstrap result
        bootstrap_df = pd.DataFrame({
            'x_new': x_new,
            'y_mean_bootstrap': y_mean_bootstrap,
            'y_lower': y_lower,
            'y_upper': y_upper
        })
        save_dataframe(bootstrap_df, f'bootstrap_{n_bootstraps}', results_dir)

    else:
        print(f"No valid splines could be generated for n_bootstraps={n_bootstraps}.")

# --- Save MSE Results ---
mse_df = pd.DataFrame(mse_results)
save_dataframe(mse_df, 'bootstrap_mse_summary', results_dir)

# --- Find Best Bootstrap Model based on MSE ---
if mse_results:
    best_model = min(mse_results, key=lambda x: x['mse'])
    best_n_bootstraps = best_model['n_bootstraps']
    best_mse = best_model['mse']
    print(f"\nBest Bootstrap Model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")

    # --- Save Best Model Information ---
    with open(os.path.join(results_dir, 'best_model.txt'), 'w') as f:
        f.write(f"Best Bootstrap Model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")
else:
    print("No valid bootstrap models were found.")
    best_n_bootstraps = None

# --- Visualization of All Bootstrap Curves ---
fig_all, ax_all = plt.subplots(figsize=(20, 7))

# Plot all bootstrap curves (semitransparent)
for n_bootstraps, result in bootstrap_results.items():
    for y_curve in result['y_values']:
        ax_all.plot(result['x_new'], y_curve, alpha=0.1, linewidth=0.5, color='gray')

# Plot data points
ax_all.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
ax_all.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
ax_all.axhline(y=VADM_actual, color='blueviolet', linestyle='-.', linewidth=0.4, label='VADM today (8.22 × 10²² Am²)')

ax_all.set_xlim(0, 5)
ax_all.set_ylim(0, 20)
ax_all.set_title('All Bootstrap Curves')
ax_all.set_ylabel('VADM')
ax_all.set_xlabel('Age (Ma)')
ax_all.legend()
plt.tight_layout()
save_plot(fig_all, 'all_bootstrap_curves', plots_dir)
plt.close()




# --- Visualización de la mejor curva de bootstrap ---
if best_n_bootstraps:
    best_result = bootstrap_results[best_n_bootstraps]
    fig_best, ax_best = plt.subplots(figsize=(20, 7))

    # Plot data points
    ax_best.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
    ax_best.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
    ax_best.axhline(y=VADM_actual, color='blueviolet', linestyle='-.', linewidth=0.4, label='VADM today (8.22 × 10²² Am²)')
    ax_best.plot(df_plot['window'], df_plot['mean_VADM'], color='cornflowerblue', linestyle='-', linewidth=1.2, label='Mean VADM 100 ka RW')

    # Plot mean curve and confidence interval
    ax_best.plot(best_result['x_new'], best_result['y_mean'], color='cyan', linewidth=2,
                label=f'Best Bootstrap (n={best_n_bootstraps}, MSE={best_mse:.4f})')
    ax_best.fill_between(best_result['x_new'], best_result['y_lower'], best_result['y_upper'],
                        color='lightblue', alpha=0.4, label='95% CI')

    ax_best.set_xlim(0, 5)
    ax_best.set_ylim(0, 20)
    ax_best.set_title(f'Best Bootstrap Model: n_bootstraps={best_n_bootstraps}')
    ax_best.set_ylabel('VADM')
    ax_best.set_xlabel('Age (Ma)')
    ax_best.grid(True, linestyle='--', alpha=0.6)
    ax_best.legend()

    plt.tight_layout()
    save_plot(fig_best, 'best_bootstrap_model', plots_dir)
    plt.show()

    # --- Print Summary of Metrics ---
    print("\n--- Summary of Bootstrap Results ---")
    print(f"Best model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")

    # Optional: plot other metrics or summaries for all models
    mse_fig, mse_ax = plt.subplots(figsize=(10, 6))
    mse_ax.plot([item['n_bootstraps'] for item in mse_results],
               [item['mse'] for item in mse_results], 'o-')
    mse_ax.set_xlabel('Number of Bootstrap Iterations')
    mse_ax.set_ylabel('MSE')
    mse_ax.set_title('MSE vs. Number of Bootstrap Iterations')
    mse_ax.grid(True)
    save_plot(mse_fig, 'mse_vs_iterations', plots_dir)
    plt.show()


# --- Store Results ---
bootstrap_results = {}
mse_results = []

# --- Bootstrap Iteration ---
for n_bootstraps in n_bootstrap_iterations:
    print(f"Evaluating Bootstrap with n_bootstraps={n_bootstraps}")
    y_values = []

    for i in range(n_bootstraps):
        try:
            x_sample, y_sample = resample_with_replacement(df_plot['window'].values, df_plot['mean_VADM'].values)
            cs = CubicSpline(x_sample, y_sample, extrapolate=True)
            x_new = np.linspace(x.min(), x.max(), 200)
            y_values.append(cs(x_new))
        except Exception as e:
            print(f"Error during bootstrap iteration {i}: {e}")
            continue

    if y_values:
        y_values = np.array(y_values)
        y_mean_bootstrap = np.mean(y_values, axis=0)
        y_lower = np.percentile(y_values, 2.5, axis=0)  # 2.5% percentile
        y_upper = np.percentile(y_values, 97.5, axis=0)  # 97.5% percentile

        mse = mean_squared_error(df_plot['mean_VADM'].values, y_mean_bootstrap[:len(df_plot['mean_VADM'].values)])

        # Store Bootstrap Results
        bootstrap_results[n_bootstraps] = {
            'y_values': y_values,
            'y_mean': y_mean_bootstrap,
            'y_lower': y_lower,
            'y_upper': y_upper,
            'x_new': x_new,
            'mse': mse
        }
        mse_results.append({'n_bootstraps': n_bootstraps, 'mse': mse})

        # Save each bootstrap result
        bootstrap_df = pd.DataFrame({
            'x_new': x_new,
            'y_mean_bootstrap': y_mean_bootstrap,
            'y_lower': y_lower,
            'y_upper': y_upper
        })
        save_dataframe(bootstrap_df, f'bootstrap_{n_bootstraps}', results_dir)

    else:
        print(f"No valid splines could be generated for n_bootstraps={n_bootstraps}.")

# --- Save MSE Results ---
mse_df = pd.DataFrame(mse_results)
save_dataframe(mse_df, 'bootstrap_mse_summary', results_dir)

# --- Find Best Bootstrap Model based on MSE ---
if mse_results:
    best_model = min(mse_results, key=lambda x: x['mse'])
    best_n_bootstraps = best_model['n_bootstraps']
    best_mse = best_model['mse']
    print(f"\nBest Bootstrap Model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")

    # --- Save Best Model Information ---
    with open(os.path.join(results_dir, 'best_model.txt'), 'w') as f:
        f.write(f"Best Bootstrap Model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")
else:
    print("No valid bootstrap models were found.")
    best_n_bootstraps = None

# --- Visualization of All Bootstrap Curves ---
fig_all, ax_all = plt.subplots(figsize=(20, 7))

# Plot all bootstrap curves (semitransparent)
for n_bootstraps, result in bootstrap_results.items():
    for y_curve in result['y_values']:
        ax_all.plot(result['x_new'], y_curve, alpha=0.1, linewidth=0.5, color='gray')

# Plot data points
ax_all.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
ax_all.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
ax_all.axhline(y=VADM_actual, color='blueviolet', linestyle='-.', linewidth=0.4, label='VADM today (8.22 × 10²² Am²)')

ax_all.set_xlim(0, 5)
ax_all.set_ylim(0, 20)
ax_all.set_title('All Bootstrap Curves')
ax_all.set_ylabel('VADM')
ax_all.set_xlabel('Age (Ma)')
ax_all.legend()
plt.tight_layout()
save_plot(fig_all, 'all_bootstrap_curves', plots_dir)
plt.close()





# --- Visualización de la mejor curva de bootstrap ---
if best_n_bootstraps:


    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 7),
                                  gridspec_kw={'height_ratios': [4, 1]},
                                  sharex=True)

    # --- PRIMER SUBPLOT: TU GRAFICO PRINCIPAL ---
    ax1.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
    ax1.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
    ax1.axhline(y=VADM_actual, color='blueviolet', linestyle='-.', linewidth=0.4, label='VADM today (8.22 × 10²² Am²)')
    ax1.plot(df_plot['window'], df_plot['mean_VADM'], color='cornflowerblue', linestyle='-', linewidth=1.2, label='Mean VADM 100 ka RW')

    ax1.plot(best_result['x_new'], best_result['y_mean'], color='cyan', linewidth=2,
        label=f'Best Bootstrap (n={best_n_bootstraps}, MSE={best_mse:.4f})')
    ax1.fill_between(best_result['x_new'], best_result['y_lower'], best_result['y_upper'],
        color='lightblue', alpha=0.4, label='95% CI')

    ax1.set_xlim(0, 5)
    ax1.set_ylim(0, 20)
    ax1.set_title(f'Best Bootstrap Model: n_bootstraps={best_n_bootstraps}')
    ax1.set_ylabel('VADM')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(fontsize=10, loc='upper right')

    # --- SEGUNDO SUBPLOT: BARRAS DE SUBCHRONES (TIMELINE) ---
    colores = ["black", "white"] * (len(subchrones) // 2 + 1)
    for i, (subcron, datos) in enumerate(subchrones.items()):
        ax2.barh(0, width=datos["techo"] - datos["base"], left=datos["base"], height=0.6,
                color=colores[i % len(colores)], edgecolor="black")
        # Texto centrado
        ax2.text((datos["base"] + datos["techo"]) / 2, 0, subcron, ha="center", va="center",
                color="white" if colores[i % len(colores)] == "black" else "black",
                fontsize=7, rotation=90)

    ax2.yaxis.set_visible(False)
    ax2.set_xlabel('AGE (Ma)', fontsize=12)
    for spine in ax2.spines.values():
        spine.set_visible(False)
    ax2.set_xlim(0, 5)

    # --- Empalme visual atractivo ---
    plt.subplots_adjust(hspace=0.08)
    plt.tight_layout()

    # --- Guardar y mostrar ---
    save_plot(fig, "best_bootstrap_model_with_subchrones", plots_dir)
    plt.show()


    # --- Print Summary of Metrics ---
    print("\n--- Summary of Bootstrap Results ---")
    print(f"Best model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")

    # Optional: plot other metrics or summaries for all models
    mse_fig, mse_ax = plt.subplots(figsize=(10, 6))
    mse_ax.plot([item['n_bootstraps'] for item in mse_results],
               [item['mse'] for item in mse_results], 'o-')
    mse_ax.set_xlabel('Number of Bootstrap Iterations')
    mse_ax.set_ylabel('MSE')
    mse_ax.set_title('MSE vs. Number of Bootstrap Iterations')
    mse_ax.grid(True)
    save_plot(mse_fig, 'mse_vs_iterations', plots_dir)
    plt.show()














# --- Store Results ---
bootstrap_results = {}

# --- Bootstrap Iteration ---
for n_bootstraps in n_bootstrap_iterations:
    print(f"Evaluating Bootstrap with n_bootstraps={n_bootstraps}")
    y_values = []

    for i in range(n_bootstraps):
        try:
            x_sample, y_sample = resample_with_replacement(df_plot['window'].values, df_plot['mean_VADM'].values)
            cs = CubicSpline(x_sample, y_sample, extrapolate=True)
            x_new = np.linspace(x.min(), x.max(), 200)
            y_values.append(cs(x_new))
        except Exception as e:
            print(f"Error during bootstrap iteration {i}: {e}")
            continue

    if y_values:
        y_values = np.array(y_values)
        y_mean_bootstrap = np.mean(y_values, axis=0)
        mse = mean_squared_error(df_plot['mean_VADM'].values, y_mean_bootstrap[:len(df_plot['mean_VADM'].values)])

        # Store Bootstrap Results
        bootstrap_results[n_bootstraps] = {
            'y_values': y_values,
            'x_new': x_new
        }
        mse_results.append({'n_bootstraps': n_bootstraps, 'mse': mse})

        # Save each bootstrap result
        bootstrap_df = pd.DataFrame({
            'x_new': x_new,
            'y_mean_bootstrap': y_mean_bootstrap
        })
        save_dataframe(bootstrap_df, f'bootstrap_{n_bootstraps}', results_dir)

    else:
        print(f"No valid splines could be generated for n_bootstraps={n_bootstraps}.")

# --- Save MSE Results ---
mse_df = pd.DataFrame(mse_results)
save_dataframe(mse_df, 'bootstrap_mse_summary', results_dir)

# --- Visualization ---
fig, ax = plt.subplots(figsize=(20, 7))

# Plot Bootstrap Curves
for n_bootstraps, result in bootstrap_results.items():
    for y_curve in result['y_values']:
        ax.plot(result['x_new'], y_curve, alpha=0.4, linewidth=0.8)

# Plot data points
ax.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
ax.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
ax.axhline(y=VADM_actual, color='blueviolet', linestyle='-.', linewidth=0.4, label='VADM today (8.22 × 10²² Am²)')

ax.set_xlim(0, 5)
ax.set_ylim(0, 20)
ax.set_title('Bootstrap Curves')
ax.set_ylabel('VADM')
ax.set_xlabel('Age (Ma)')
ax.legend()

plt.tight_layout()
#plt.savefig('/content/drive/MyDrive/Results/bootstrap_curves.png', dpi=300)
save_plot(fig, 'bootstrap_curves', plots_dir)

plt.show()

# --- Print Summary of Metrics ---
#print("\n--- Summary of Bootstrap Results ---")
#for n_bootstraps, result in bootstrap_results.items():
#    print(f"n_bootstraps={n_bootstraps}: MSE={result['mse']:.4f}")


# --- Store Results ---
bootstrap_results = {}
# --- Bootstrap Iteration ---
for n_bootstraps in n_bootstrap_iterations:
    print(f"Evaluating Bootstrap with n_bootstraps={n_bootstraps}")
    y_values = []

    for i in range(n_bootstraps):
        try:
            x_sample, y_sample = resample_with_replacement(df_plot['window'].values, df_plot['mean_VADM'].values)
            cs = CubicSpline(x_sample, y_sample, extrapolate=True)
            x_new = np.linspace(x.min(), x.max(), 200)
            y_values.append(cs(x_new))
        except Exception as e:
            print(f"Error during bootstrap iteration {i}: {e}")
            continue

    if y_values:
        y_values = np.array(y_values)
        y_mean_bootstrap = np.mean(y_values, axis=0)
        mse = mean_squared_error(df_plot['mean_VADM'].values, y_mean_bootstrap[:len(df_plot['mean_VADM'].values)])

        # Store Bootstrap Results
        bootstrap_results[n_bootstraps] = {
            'y_values': y_values,
            'x_new': x_new
        }
        mse_results.append({'n_bootstraps': n_bootstraps, 'mse': mse})

        # Save each bootstrap result
        bootstrap_df = pd.DataFrame({
            'x_new': x_new,
            'y_mean_bootstrap': y_mean_bootstrap
        })
        save_dataframe(bootstrap_df, f'bootstrap_{n_bootstraps}', results_dir)

    else:
        print(f"No valid splines could be generated for n_bootstraps={n_bootstraps}.")

# --- Save MSE Results ---
mse_df = pd.DataFrame(mse_results)
save_dataframe(mse_df, 'bootstrap_mse_summary', results_dir)

# --- Visualization ---
fig, ax = plt.subplots(figsize=(20, 7))

# Plot Bootstrap Curves
for n_bootstraps, result in bootstrap_results.items():
    for y_curve in result['y_values']:
        ax.plot(result['x_new'], y_curve, alpha=0.4, linewidth=0.8)

# Plot data points
ax.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
ax.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
ax.axhline(y=VADM_actual, color='blueviolet', linestyle='--', linewidth=0.6, label='VADM today (8.22 × 10²² Am²)')

ax.set_xlim(0, 5)
ax.set_ylim(0, 20)
ax.set_title('Bootstrap Curves')
ax.set_ylabel('VADM')
ax.set_xlabel('Age (Ma)')
ax.legend()

plt.tight_layout()
#plt.savefig('/content/drive/MyDrive/Results/bootstrap_curves.png', dpi=300)
save_plot(fig, 'bootstrap_curves', plots_dir)

plt.show()

# --- Print Summary of Metrics ---
#print("\n--- Summary of Bootstrap Results ---")
#for n_bootstraps, result in bootstrap_results.items():
#    print(f"n_bootstraps={n_bootstraps}: MSE={result['mse']:.4f}")



# --- Store Results ---
bootstrap_results = {}
mse_results = []

# --- Bootstrap Iteration ---
for n_bootstraps in n_bootstrap_iterations:
    print(f"Evaluating Bootstrap with n_bootstraps={n_bootstraps}")
    y_values = []

    for i in range(n_bootstraps):
        try:
            x_sample, y_sample = resample_with_replacement(df_plot['window'].values, df_plot['mean_VADM'].values)
            cs = CubicSpline(x_sample, y_sample, extrapolate=True)
            x_new = np.linspace(x.min(), x.max(), 200)
            y_values.append(cs(x_new))
        except Exception as e:
            print(f"Error during bootstrap iteration {i}: {e}")
            continue

    if y_values:
        y_values = np.array(y_values)
        y_mean_bootstrap = np.mean(y_values, axis=0)
        y_lower = np.percentile(y_values, 2.5, axis=0)  # 2.5% percentile
        y_upper = np.percentile(y_values, 97.5, axis=0)  # 97.5% percentile

        mse = mean_squared_error(df_plot['mean_VADM'].values, y_mean_bootstrap[:len(df_plot['mean_VADM'].values)])

        # Store Bootstrap Results
        bootstrap_results[n_bootstraps] = {
            'y_values': y_values,
            'y_mean': y_mean_bootstrap,
            'y_lower': y_lower,
            'y_upper': y_upper,
            'x_new': x_new,
            'mse': mse
        }
        mse_results.append({'n_bootstraps': n_bootstraps, 'mse': mse})

        # Save each bootstrap result
        bootstrap_df = pd.DataFrame({
            'x_new': x_new,
            'y_mean_bootstrap': y_mean_bootstrap,
            'y_lower': y_lower,
            'y_upper': y_upper
        })
        save_dataframe(bootstrap_df, f'bootstrap_{n_bootstraps}', results_dir)

    else:
        print(f"No valid splines could be generated for n_bootstraps={n_bootstraps}.")

# --- Save MSE Results ---
mse_df = pd.DataFrame(mse_results)
save_dataframe(mse_df, 'bootstrap_mse_summary', results_dir)

# --- Find Best Bootstrap Model based on MSE ---
if mse_results:
    best_model = min(mse_results, key=lambda x: x['mse'])
    best_n_bootstraps = best_model['n_bootstraps']
    best_mse = best_model['mse']
    print(f"\nBest Bootstrap Model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")

    # --- Save Best Model Information ---
    with open(os.path.join(results_dir, 'best_model.txt'), 'w') as f:
        f.write(f"Best Bootstrap Model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")
else:
    print("No valid bootstrap models were found.")
    best_n_bootstraps = None

# --- Visualization of All Bootstrap Curves ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 7),
                                  gridspec_kw={'height_ratios': [4, 1]},
                                  sharex=True)
# Plot all bootstrap curves (semitransparent)
for n_bootstraps, result in bootstrap_results.items():
     # --- PRIMER SUBPLOT: GRAFICO PRINCIPAL --
    for y_curve in result['y_values']:
        ax1.plot(result['x_new'], y_curve, alpha=0.1, linewidth=0.4, color='gray')

ax1.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
ax1.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
ax1.axhline(y=VADM_actual, color='blueviolet', linestyle='--', linewidth=0.6, label='VADM today (8.22 × 10²² Am²)')

ax1.plot(df_plot['window'], df_plot['mean_VADM'], color='blue', linestyle='-', linewidth=1.2, label='Mean VADM 100 ka RW')
ax1.set_xlim(0, 5)
ax1.set_ylim(0, 20)
ax1.set_title('All Bootstrap Curves')
ax1.set_ylabel('VADM')
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.legend(fontsize=10, loc='upper right')

# --- SEGUNDO SUBPLOT: BARRAS DE SUBCHRONES (TIMELINE) ---
colores = ["black", "white"] * (len(subchrones) // 2 + 1)
for i, (subcron, datos) in enumerate(subchrones.items()):
    ax2.barh(0, width=datos["techo"] - datos["base"], left=datos["base"], height=0.6,
            color=colores[i % len(colores)], edgecolor="black")
    # Texto centrado
    ax2.text((datos["base"] + datos["techo"]) / 2, 0, subcron, ha="center", va="center",
            color="white" if colores[i % len(colores)] == "black" else "black",
            fontsize=7, rotation=90)

ax2.yaxis.set_visible(False)
ax2.set_xlabel('AGE (Ma)', fontsize=12)
for spine in ax2.spines.values():
    spine.set_visible(False)
ax2.set_xlim(0, 5)

# --- Empalme visual atractivo ---
plt.subplots_adjust(hspace=0.08)
plt.tight_layout()
save_plot(fig, 'all_bootstrap_curves', plots_dir)
plt.close()



# --- Visualización de la mejor curva de bootstrap ---
if best_n_bootstraps:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 7),
                                  gridspec_kw={'height_ratios': [4, 1]},
                                  sharex=True)
    best_result = bootstrap_results[best_n_bootstraps]
    fig_best, ax_best = plt.subplots(figsize=(20, 7))

    # --- PRIMER SUBPLOT: TU GRAFICO PRINCIPAL ---
    ax1.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
    ax1.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
    ax1.axhline(y=VADM_actual, color='blueviolet', linestyle='-.', linewidth=0.4, label='VADM today (8.22 × 10²² Am²)')
    ax1.plot(df_plot['window'], df_plot['mean_VADM'], color='blue', linestyle='-', linewidth=1.2, label='Mean VADM 100 ka RW')

    ax1.plot(best_result['x_new'], best_result['y_mean'], color='cyan', linewidth=2,
        label=f'Best Bootstrap (n={best_n_bootstraps}, MSE={best_mse:.4f})')
    ax1.fill_between(best_result['x_new'], best_result['y_lower'], best_result['y_upper'],
        color='lightblue', alpha=0.4, label='95% CI')

    ax1.set_xlim(0, 5)
    ax1.set_ylim(0, 20)
    ax1.set_title(f'Best Bootstrap Model: n_bootstraps={best_n_bootstraps}')
    ax1.set_ylabel('VADM')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(fontsize=10, loc='upper right')

    # --- SEGUNDO SUBPLOT: BARRAS DE SUBCHRONES (TIMELINE) ---
    colores = ["black", "white"] * (len(subchrones) // 2 + 1)
    for i, (subcron, datos) in enumerate(subchrones.items()):
        ax2.barh(0, width=datos["techo"] - datos["base"], left=datos["base"], height=0.6,
                color=colores[i % len(colores)], edgecolor="black")
        # Texto centrado
        ax2.text((datos["base"] + datos["techo"]) / 2, 0, subcron, ha="center", va="center",
                color="white" if colores[i % len(colores)] == "black" else "black",
                fontsize=7, rotation=90)

    ax2.yaxis.set_visible(False)
    ax2.set_xlabel('AGE (Ma)', fontsize=12)
    for spine in ax2.spines.values():
        spine.set_visible(False)
    ax2.set_xlim(0, 5)

    # --- Empalme visual atractivo ---
    plt.subplots_adjust(hspace=0.08)
    #plt.tight_layout()

     # --- Guardar y mostrar ---
    save_plot(fig, "best_bootstrap_model_with_subchrones", plots_dir)
    plt.show()

    # --- Print Summary of Metrics ---
    print("\n--- Summary of Bootstrap Results ---")
    print(f"Best model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")

    # Optional: plot other metrics or summaries for all models
    mse_fig, mse_ax = plt.subplots(figsize=(10, 6))
    mse_ax.plot([item['n_bootstraps'] for item in mse_results],
               [item['mse'] for item in mse_results], 'o-')
    mse_ax.set_xlabel('Number of Bootstrap Iterations')
    mse_ax.set_ylabel('MSE')
    mse_ax.set_title('MSE vs. Number of Bootstrap Iterations')
    mse_ax.grid(True)
    save_plot(mse_fig, 'mse_vs_iterations', plots_dir)
    plt.show()









# --- Store Results ---
bootstrap_results = {}
mse_results = []

# --- Bootstrap Iteration ---
for n_bootstraps in n_bootstrap_iterations:
    print(f"Evaluating Bootstrap with n_bootstraps={n_bootstraps}")
    y_values = []

    for i in range(n_bootstraps):
        try:
            x_sample, y_sample = resample_with_replacement(df_plot['window'].values, df_plot['mean_VADM'].values)
            cs = CubicSpline(x_sample, y_sample, extrapolate=True)
            x_new = np.linspace(x.min(), x.max(), 200)
            y_values.append(cs(x_new))
        except Exception as e:
            print(f"Error during bootstrap iteration {i}: {e}")
            continue

    if y_values:
        y_values = np.array(y_values)
        y_mean_bootstrap = np.mean(y_values, axis=0)
        y_lower = np.percentile(y_values, 2.5, axis=0)  # 2.5% percentile
        y_upper = np.percentile(y_values, 97.5, axis=0)  # 97.5% percentile

        mse = mean_squared_error(df_plot['mean_VADM'].values, y_mean_bootstrap[:len(df_plot['mean_VADM'].values)])

        # Store Bootstrap Results
        bootstrap_results[n_bootstraps] = {
            'y_values': y_values,
            'y_mean': y_mean_bootstrap,
            'y_lower': y_lower,
            'y_upper': y_upper,
            'x_new': x_new,
            'mse': mse
        }
        mse_results.append({'n_bootstraps': n_bootstraps, 'mse': mse})

        # Save each bootstrap result
        bootstrap_df = pd.DataFrame({
            'x_new': x_new,
            'y_mean_bootstrap': y_mean_bootstrap,
            'y_lower': y_lower,
            'y_upper': y_upper
        })
        save_dataframe(bootstrap_df, f'bootstrap_{n_bootstraps}', results_dir)

    else:
        print(f"No valid splines could be generated for n_bootstraps={n_bootstraps}.")

# --- Save MSE Results ---
mse_df = pd.DataFrame(mse_results)
save_dataframe(mse_df, 'bootstrap_mse_summary', results_dir)

# --- Find Best Bootstrap Model based on MSE ---
if mse_results:
    best_model = min(mse_results, key=lambda x: x['mse'])
    best_n_bootstraps = best_model['n_bootstraps']
    best_mse = best_model['mse']
    print(f"\nBest Bootstrap Model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")

    # --- Save Best Model Information ---
    with open(os.path.join(results_dir, 'best_model.txt'), 'w') as f:
        f.write(f"Best Bootstrap Model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")
else:
    print("No valid bootstrap models were found.")
    best_n_bootstraps = None



# --- Visualization of All Bootstrap Curves ---
fig_all, ax_all = plt.subplots(figsize=(20, 7))

# Plot all bootstrap curves (semitransparent)
for n_bootstraps, result in bootstrap_results.items():
    for y_curve in result['y_values']:
        ax_all.plot(result['x_new'], y_curve, alpha=0.1, linewidth=0.5, color='gray')

# Plot data points
ax_all.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
ax_all.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
ax_all.axhline(y=VADM_actual, color='blueviolet', linestyle='-.', linewidth=0.4, label='VADM today (8.22 × 10²² Am²)')

ax_all.set_xlim(0, 5)
ax_all.set_ylim(0, 20)
ax_all.set_title('All Bootstrap Curves')
ax_all.set_ylabel('VADM')
ax_all.set_xlabel('Age (Ma)')
ax_all.legend()
plt.tight_layout()
save_plot(fig_all, 'all_bootstrap_curves', plots_dir)
plt.close()
# --- Visualization of All Bootstrap Curves ---
fig_all, ax_all = plt.subplots(figsize=(20, 7))

# Plot all bootstrap curves (semitransparent)
for n_bootstraps, result in bootstrap_results.items():
    for y_curve in result['y_values']:
        ax_all.plot(result['x_new'], y_curve, alpha=0.1, linewidth=0.5, color='gray')

# Plot data points
ax_all.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
ax_all.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
ax_all.axhline(y=VADM_actual, color='blueviolet', linestyle='--', linewidth=0.6, label='VADM today (8.22 × 10²² Am²)')

ax_all.set_xlim(0, 5)
ax_all.set_ylim(0, 20)
ax_all.set_title('All Bootstrap Curves')
ax_all.set_ylabel('VADM')
ax_all.set_xlabel('Age (Ma)')
ax_all.legend()
plt.tight_layout()
save_plot(fig_all, 'all_bootstrap_curves', plots_dir)
plt.close()
# --- Visualización de la mejor curva de bootstrap ---
if best_n_bootstraps:


    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 7),
                                  gridspec_kw={'height_ratios': [4, 1]},
                                  sharex=True)

    # --- PRIMER SUBPLOT: TU GRAFICO PRINCIPAL ---
    ax1.scatter(data_vadm['AGE'], data_vadm['VADM'], label='PINT', facecolors='none', edgecolors='black', s=8, linewidth=0.4)
    ax1.scatter(new_data['AGE'], new_data['VADM'], label='MCADAM', facecolors='red', edgecolors='red', s=8, linewidth=0.4)
    ax1.axhline(y=VADM_actual, color='blueviolet', linestyle='--', linewidth=0.4, label='VADM today (8.22 × 10²² Am²)')
    ax1.plot(df_plot['window'], df_plot['mean_VADM'], color='blue', linestyle='-', linewidth=1.2, label='Mean VADM 100 ka RW')

    ax1.plot(best_result['x_new'], best_result['y_mean'], color='cyan', linewidth=2,
        label=f'Best Bootstrap (n={best_n_bootstraps}, MSE={best_mse:.4f})')
    ax1.fill_between(best_result['x_new'], best_result['y_lower'], best_result['y_upper'],
        color='lightblue', alpha=0.4, label='95% CI')

    ax1.set_xlim(0, 5)
    ax1.set_ylim(0, 20)
    ax1.set_title(f'Best Bootstrap Model: n_bootstraps={best_n_bootstraps}')
    ax1.set_ylabel('VADM')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(fontsize=10, loc='upper right')

    # --- SEGUNDO SUBPLOT: BARRAS DE SUBCHRONES (TIMELINE) ---
    colores = ["black", "white"] * (len(subchrones) // 2 + 1)
    for i, (subcron, datos) in enumerate(subchrones.items()):
        ax2.barh(0, width=datos["techo"] - datos["base"], left=datos["base"], height=0.6,
                color=colores[i % len(colores)], edgecolor="black")
        # Texto centrado
        ax2.text((datos["base"] + datos["techo"]) / 2, 0, subcron, ha="center", va="center",
                color="white" if colores[i % len(colores)] == "black" else "black",
                fontsize=7, rotation=90)

    ax2.yaxis.set_visible(False)
    ax2.set_xlabel('AGE (Ma)', fontsize=12)
    for spine in ax2.spines.values():
        spine.set_visible(False)
    ax2.set_xlim(0, 5)

    # --- Empalme visual atractivo ---
    plt.subplots_adjust(hspace=0.08)
    plt.tight_layout()

    # --- Guardar y mostrar ---
    save_plot(fig, "best_bootstrap_model_with_subchrones", plots_dir)
    plt.show()


    # --- Print Summary of Metrics ---
    print("\n--- Summary of Bootstrap Results ---")
    print(f"Best model: n_bootstraps={best_n_bootstraps}, MSE={best_mse:.4f}")

    # Optional: plot other metrics or summaries for all models
    mse_fig, mse_ax = plt.subplots(figsize=(10, 6))
    mse_ax.plot([item['n_bootstraps'] for item in mse_results],
               [item['mse'] for item in mse_results], 'o-')
    mse_ax.set_xlabel('Number of Bootstrap Iterations')
    mse_ax.set_ylabel('MSE')
    mse_ax.set_title('MSE vs. Number of Bootstrap Iterations')
    mse_ax.grid(True)
    save_plot(mse_fig, 'mse_vs_iterations', plots_dir)
    plt.show()

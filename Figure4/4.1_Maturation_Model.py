# ------------------------------------------------------------------------------
# Title: Maturation Model for CM Data
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script trains a machine learning model to predict CM (cardiomyocyte) maturation 
# using multiple datasets. It includes data preprocessing, feature engineering, 
# model training, evaluation, and SHAP analysis.
#
# Key Functions:
# - Loads and preprocesses CM datasets
# - Handles missing values through imputation
# - Builds predictive models using Linear Regression and Random Forest
# - Evaluates model performance using RMSE and R-squared
# - Performs SHAP analysis to interpret feature importance
#
# Dependencies:
# - pandas, scikit-learn, shap
#
# Output:
# - Trained models
# - Model performance metrics
# - SHAP feature importance visualization
#
# Usage:
# - Ensure CM data files are present in the `./data` directory.
# - Run the script to train and evaluate the models.
# ------------------------------------------------------------------------------

import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer

# ---------------------------------------------
# 1. Load & Preprocess the Data
# ---------------------------------------------

data_dir = "./data"  
        
cm_dataframes = {}  # Store only cm_data
for file_name in os.listdir(data_dir):
    if file_name.endswith(".csv") and "_cm_" in file_name:
        file_path = os.path.join(data_dir, file_name)
        df = pd.read_csv(file_path, index_col=0)  # Load each file and set the first column as index
        key_name = file_name.replace("_data.csv", "").replace(".csv", "")
        cm_dataframes[key_name] = df

# Display loaded CM datasets
print("Loaded CM datasets:")
for key, df in cm_dataframes.items():
    print(f"{key} -> {df.shape[0]} rows, {df.shape[1]} columns")
    
# ------------------------------------------------------------------------------
# Step 2: Handle Missing Values
# ------------------------------------------------------------------------------

def check_and_fix_missing_values(cm_dataframes):
    """Check and impute missing values in all CM datasets."""
    
    for file_name, df in cm_dataframes.items():
        # Identify missing values
        numeric_features = df.select_dtypes(include=['float64', 'int64']).columns
        categorical_features = df.select_dtypes(include=['object']).columns

        # Impute missing numeric values with median
        df[numeric_features] = df[numeric_features].fillna(df[numeric_features].median())

        # Impute missing categorical values with most frequent value
        df[categorical_features] = df[categorical_features].fillna(df[categorical_features].mode().iloc[0])

        print(f"Processed missing values for {file_name}")

    return cm_dataframes

# Apply missing value handling
cm_dataframes = check_and_fix_missing_values(cm_dataframes)

# ------------------------------------------------------------------------------
# Step 3: Feature Engineering & Model Training
# ------------------------------------------------------------------------------

# ligand_receptor_only
 
BASE_SAVE_DIR = "saved_models_ligand_receptor_only"
TRAIN_TEST_DIR = os.path.join(BASE_SAVE_DIR, "train_test_data")
os.makedirs(TRAIN_TEST_DIR, exist_ok=True)

target = "Maturation_Features1"

def train_rf_with_shap(df, dataset_name):
    """Train a RandomForest model using only 'commot.' ligand-receptor pairs and run SHAP analysis."""
    print(f"Processing dataset: {dataset_name}")

    # **Filter for commot ligand-receptor only features**
    lr_pattern = re.compile(r"^commot\.user_database\.[^.]+\.[^.]+$")  # Matches 'commot.user_database.Angptl2.Itga5_Itgb1'

    commot_lr_features = [col for col in df.columns if lr_pattern.match(col)]

    if target in df.columns:
        commot_lr_features.append(target)

    if not commot_lr_features:
        print(f"⚠️ Skipping {dataset_name} - No commot ligand-receptor features found.")
        return None

    df_filtered = df[commot_lr_features]

    if target not in df_filtered.columns:
        print(f"⚠️ Skipping {dataset_name} - Target '{target}' not found.")
        return None

    # Fill missing values in commot ligand-receptor features with 0
    df_filtered = df_filtered.fillna(0)

    # Paths for train/test split
    X_train_path = os.path.join(TRAIN_TEST_DIR, f"X_train_{dataset_name}.joblib")
    X_test_path = os.path.join(TRAIN_TEST_DIR, f"X_test_{dataset_name}.joblib")
    y_train_path = os.path.join(TRAIN_TEST_DIR, f"y_train_{dataset_name}.joblib")
    y_test_path = os.path.join(TRAIN_TEST_DIR, f"y_test_{dataset_name}.joblib")

    # Train-test split
    X = df_filtered.drop(columns=[target])
    y = df_filtered[target]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Save train/test data
    joblib.dump(X_train, X_train_path)
    joblib.dump(X_test, X_test_path)
    joblib.dump(y_train, y_train_path)
    joblib.dump(y_test, y_test_path)

    # Define preprocessing (Standardizing numerical features)
    preprocessor = ColumnTransformer(
        transformers=[
            ('num', StandardScaler(), X_train.columns)
        ]
    )

    # Train model
    model = Pipeline(steps=[
        ('preprocessor', preprocessor),
        ('regressor', RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1))
    ])
    
    model.fit(X_train, y_train)

    # Evaluate model
    y_pred = model.predict(X_test)
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)
    print(f"{dataset_name}: RMSE={rmse:.4f}, R²={r2:.4f}")
    model_path = os.path.join(BASE_SAVE_DIR, f"random_forest_{dataset_name}.joblib")
    joblib.dump(model, model_path)
    print(f"💾 Model saved: {model_path}")

    # -----------------------------
    # SHAP Analysis (Fast Mode)
    # -----------------------------
    print(f"\n🔍 Running SHAP analysis for {dataset_name} (Fast Mode)")

    # Define PDF report for this dataset
    pdf_path = os.path.join(BASE_SAVE_DIR, f"shap_analysis_{dataset_name}.pdf")

    with PdfPages(pdf_path) as pdf:
        X_test_numeric = model.named_steps['preprocessor'].transform(X_test)
        feature_names = X_test.columns

        explainer = shap.TreeExplainer(
            model.named_steps['regressor'], 
            feature_perturbation="tree_path_dependent",  # Fastest for RF
            approximate=True
        )

        X_test_sampled = shap.sample(X_test_numeric, 300)  # Sample 300 rows for efficiency
        shap_values = explainer.shap_values(X_test_sampled)

        importance_df = pd.DataFrame({
            "Feature": feature_names,
            "Mean |SHAP|": np.abs(shap_values).mean(axis=0)
        }).sort_values(by="Mean |SHAP|", ascending=False)

        # **1️ SHAP Importance Plot (Top 30 Features)**
        plt.figure(figsize=(16, 7))  # **Increased Size**
        plt.barh(importance_df["Feature"].head(30), importance_df["Mean |SHAP|"].head(30), color="blue")
        plt.xlabel("Mean |SHAP| Value")
        plt.ylabel("Features")
        plt.title(f"Top 30 SHAP Feature Importances: {dataset_name}")
        plt.gca().invert_yaxis()
        pdf.savefig()
        plt.close()

        # **2️ SHAP Summary Plot (Top 30 Features)**
        plt.figure(figsize=(12, 6))
        shap.summary_plot(shap_values, X_test_sampled, feature_names=feature_names, show=False, max_display=30)
        pdf.savefig()
        plt.close()

        # **3️ SHAP Dependence Plots (Top 30 Features)**
        for feature in importance_df["Feature"].head(30):
            plt.figure(figsize=(8, 5))
            shap.dependence_plot(feature, shap_values, X_test_sampled, feature_names=feature_names, show=False)
            pdf.savefig()
            plt.close()

    print(f" SHAP analysis completed for {dataset_name}, saved to {pdf_path}")

for dataset_name, df in cm_dataframes.items():
    train_rf_with_shap(df, dataset_name)



# ------------------------------------------------------------------------------
# Step 4: Post training Evaluation and SHAP Analysis for Feature Importance
# ------------------------------------------------------------------------------

def test_saved_model(dataset_name):
    """Loads saved model and test set, evaluates performance."""
    
    model_path = os.path.join(BASE_SAVE_DIR, f"random_forest_{dataset_name}.joblib")
    X_test_path = os.path.join(TRAIN_TEST_DIR, f"X_test_{dataset_name}.joblib")
    y_test_path = os.path.join(TRAIN_TEST_DIR, f"y_test_{dataset_name}.joblib")

    if not os.path.exists(model_path):
        print(f"⚠️ Model not found: {model_path}")
        return

    if not os.path.exists(X_test_path) or not os.path.exists(y_test_path):
        print(f"⚠️ Test data missing for {dataset_name}")
        return

    X_test = joblib.load(X_test_path)
    y_test = joblib.load(y_test_path)

    model = joblib.load(model_path)

    y_pred = model.predict(X_test)

    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)

    print(f"{dataset_name}: RMSE={rmse:.4f}, R²={r2:.4f}")


for dataset_name in os.listdir(TRAIN_TEST_DIR):
    if dataset_name.startswith("X_test_") and dataset_name.endswith(".joblib"):
        dataset_name = dataset_name.replace("X_test_", "").replace(".joblib", "")
        test_saved_model(dataset_name)


performance_results = []
for filename in os.listdir(TRAIN_TEST_DIR):
    if filename.startswith("X_test_") and filename.endswith(".joblib"):
        dataset_name = filename.replace("X_test_", "").replace(".joblib", "")

        X_test_path = os.path.join(TRAIN_TEST_DIR, f"X_test_{dataset_name}.joblib")
        y_test_path = os.path.join(TRAIN_TEST_DIR, f"y_test_{dataset_name}.joblib")
        model_path = os.path.join(BASE_SAVE_DIR, f"random_forest_{dataset_name}.joblib")

        if not os.path.exists(model_path) or not os.path.exists(X_test_path) or not os.path.exists(y_test_path):
            continue  # Skip if any required file is missing

        X_test = joblib.load(X_test_path)
        y_test = joblib.load(y_test_path)
        model = joblib.load(model_path)

        y_pred = model.predict(X_test)

        rmse = mean_squared_error(y_test, y_pred, squared=False)
        r2 = r2_score(y_test, y_pred)

        performance_results.append({"Dataset": dataset_name, "RMSE": rmse, "R^2": r2})

performance_df = pd.DataFrame(performance_results)
desired_order = ["p0_cm", "p7_cm", "p14_cm", "p21_cm"]
performance_df = performance_df.set_index("Dataset").loc[desired_order].reset_index()
plot_save_path = os.path.join(BASE_SAVE_DIR, "model_performance_R2.pdf")
fig, ax = plt.subplots(figsize=(10, 5))
x_labels = np.arange(len(performance_df["Dataset"]))
bars_r2 = ax.bar(x_labels, performance_df["R^2"], color='red', alpha=0.6, label="R²")

for bar in bars_r2:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, height, f"{height:.2f}", ha='center', va='bottom', fontsize=10)

ax.set_ylabel("R² Score")
ax.set_xticks(x_labels)
ax.set_xticklabels(performance_df["Dataset"], ha="right")
ax.set_title("Model R² Performance Across Datasets")
ax.legend()
plt.tight_layout()
plt.savefig(plot_save_path, dpi=300, format="pdf", bbox_inches='tight')
plt.close()

print(f"Saved R² performance plot at {plot_save_path}.")


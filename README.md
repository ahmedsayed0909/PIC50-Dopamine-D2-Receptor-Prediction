# 🧬 Dopamine D2 Receptor Bioactivity Predictor

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![Streamlit](https://img.shields.io/badge/Streamlit-App-FF4B4B.svg)
![Machine Learning](https://img.shields.io/badge/Machine%20Learning-scikit--learn-orange.svg)
![Cheminformatics](https://img.shields.io/badge/Cheminformatics-RDKit-green.svg)

## 📌 Project Overview
This project is an end-to-end computational drug discovery and machine learning pipeline. It is designed to predict the bioactivity (specifically the **pIC50** value) of chemical compounds targeting the **Dopamine D2 receptor**, a critical target for treating neurological and psychiatric disorders like schizophrenia and Parkinson's disease.

The pipeline covers everything from automated data retrieval from the ChEMBL database to feature engineering (Morgan fingerprints, Lipinski's Rule of 5), exploratory data analysis, machine learning model tuning, and final deployment as a web application.

## ✨ Features
* **Automated Data Mining:** Retrieves raw bioactivity data directly from the ChEMBL database using the `chembl_webresource_client`.
* **Cheminformatics Feature Engineering:** Uses `RDKit` to calculate molecular descriptors (Molecular Weight, LogP, NumHDonors, NumHAcceptors) and 1024-bit Morgan Fingerprints from SMILES strings.
* **Exploratory Data Analysis (EDA):** Analyzes chemical space and evaluates the drug-likeness of compounds based on Lipinski's Rule of 5.
* **Machine Learning Pipeline:** Employs `LazyPredict` for initial algorithm screening, followed by rigorous hyperparameter tuning (`GridSearchCV`) of a `GradientBoostingRegressor`. 
* **Interactive Web Application:** A user-friendly Streamlit app that allows researchers to input SMILES strings and instantly predict the pIC50 bioactivity of novel compounds.

## 📂 Project Structure

The project is divided into 8 sequential Jupyter Notebooks documenting the research process, plus the final deployment script:

* **`ML_part1_drug_discovery_Dopamine_D2_receptor.ipynb`**: Data collection and downloading from the ChEMBL database.
* **`ML_part2_drug_discovery_Dopamine_D2_receptor.ipynb`**: Data cleaning, handling missing values, and data validation.
* **`ML_part3_drug_discovery_Dopamine_D2_receptor.ipynb`**: Data preprocessing, calculating Lipinski descriptors, and transforming IC50 to pIC50.
* **`ML_part4_drug_discovery_Dopamine_D2_receptor.ipynb`**: Exploratory Data Analysis (EDA) and visualizing chemical space.
* **`ML_part5_drug_discovery_Dopamine_D2_receptor.ipynb`**: Feature engineering: Generating Morgan Fingerprints for modeling.
* **`ML_part6_drug_discovery_Dopamine_D2_receptor.ipynb`**: Initial model building, applying VarianceThreshold, and screening algorithms.
* **`ML_part7_drug_discovery_Dopamine_D2_receptor.ipynb`**: Hyperparameter tuning and selection of the best-performing model.
* **`ML_part8_drug_discovery_Dopamine_D2_receptor.ipynb`**: Final model training, saving the pipeline (model + feature selector), and residual analysis.
* **`app.py`**: The Streamlit web application for real-time predictions.

## 🛠️ Tech Stack & Libraries
* **Language:** Python
* **Data Processing & EDA:** Pandas, NumPy, Matplotlib, Seaborn
* **Cheminformatics:** RDKit, ChEMBL Webresource Client
* **Machine Learning:** scikit-learn, LazyPredict, Joblib
* **Deployment:** Streamlit

## 🚀 Getting Started

```bash
git clone [https://github.com/your-username/dopamine-d2-bioactivity-predictor.git](https://github.com/your-username/dopamine-d2-bioactivity-predictor.git)
cd dopamine-d2-bioactivity-predictor
pip install pandas numpy scikit-learn rdkit chembl_webresource_client matplotlib seaborn streamlit joblib
streamlit run app.py

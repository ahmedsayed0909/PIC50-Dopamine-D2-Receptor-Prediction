import os
import streamlit as st
import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# --- 1. Page Configuration ---
st.set_page_config(page_title="Dopamine D2 Bioactivity Predictor",  layout="centered")

st.title("Dopamine D2 Receptor Bioactivity Predictor")
st.write("""
This application predicts the pIC50 bioactivity of chemical compounds against the Dopamine D2 receptor. 
Simply enter the SMILES strings of your molecules below to run them through our tuned Gradient Boosting model.
""")

# --- 2. Load Models ---
@st.cache_resource
def load_components():
    # model = joblib.load("../Final Model/final_evaluated_GBR_model.pkl")
    # selector = joblib.load("../Final Model/final_variance_selector.pkl")
    current_dir = os.path.dirname(os.path.abspath(__file__))
    model_path = os.path.join(current_dir, "..", "Final Model", "final_evaluated_GBR_model.pkl")
    selector_path = os.path.join(current_dir, "..", "Final Model", "final_variance_selector.pkl")
    model = joblib.load(model_path)
    selector = joblib.load(selector_path)

    return model, selector

model, selector = load_components()

# --- 3. Feature Engineering Function ---
def get_fingerprints(smiles_list):
    """Converts SMILES into 1024 bits + MW + LogP + NumHDonors + NumHAcceptors"""
    features = []
    valid_smiles = []
    
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is not None:
            # 1. Generate 1024-bit Morgan Fingerprint
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
            fp_array = list(fp)
            
            # 2. Calculate the 4 Lipinski Descriptors
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            num_hdonors = Descriptors.NumHDonors(mol)
            num_hacceptors = Descriptors.NumHAcceptors(mol)
            
            # 3. Combine them all into a single row
            # We add a dummy 'Unnamed: 0' (just an index number `i`) at the start 
            # because that column was present during training!
            row = [i] + fp_array + [mw, logp, num_hdonors, num_hacceptors]
            
            features.append(row)
            valid_smiles.append(smiles)
        else:
            st.warning(f"Invalid SMILES ignored: {smiles}")
            
    if not features:
        return None, None
        
    # 4. Create exact column names to match the trained Variance Selector
    col_names = ['Unnamed: 0'] + [f'Bit_{j}' for j in range(1024)] + ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
    
    df_features = pd.DataFrame(features, columns=col_names)
    return df_features, valid_smiles

# --- 4. User Interface ---
st.subheader("Input Molecules")
user_input = st.text_area(
    "Enter SMILES strings (one per line):", 
    "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c3ccc(F)cc3\nCN1CCc2cccc3c2C1Cc4ccc(O)cc34"
)

if st.button("Predict pIC50", type="primary"):
    if user_input:
        # Split the text area input by newlines and remove empty lines
        smiles_list = [s for s in user_input.split('\n') if s.strip()]
        
        with st.spinner("Calculating chemical descriptors and making predictions..."):
            # 1. Calculate all 1029 features
            features_df, valid_smiles = get_fingerprints(smiles_list)
            
            if features_df is not None:
                try:
                    # 2. Filter features using the saved VarianceThreshold
                    filtered_features = selector.transform(features_df)
                    
                    # 3. Make the Predictions
                    predictions = model.predict(filtered_features)
                    
                    # 4. Display Results
                    st.subheader("Prediction Results")
                    results_df = pd.DataFrame({
                        "Molecule (SMILES)": valid_smiles,
                        "Predicted pIC50": np.round(predictions, 4)
                    })
                    
                    st.dataframe(results_df, use_container_width=True)
                    
                    # Download button
                    csv = results_df.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        label="Download Results as CSV",
                        data=csv,
                        file_name="dopamine_d2_predictions.csv",
                        mime="text/csv",
                    )
                except ValueError as e:
                    st.error(f"Error during prediction: {e}")
    else:
        st.error("Please enter at least one SMILES string.")
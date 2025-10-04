import streamlit as st
import pandas as pd
import subprocess
import os
import time
import datetime

# --- Helper Function for Logging ---
def log_message(message):
    """Appends a timestamped message to the debug.log file."""
    with open("debug.log", "a") as log_file:
        log_file.write(f"{datetime.datetime.now()}: {message}\n")

# --- Page Configuration ---
st.set_page_config(page_title="GAPO", layout="wide", initial_sidebar_state="expanded")

st.title("🧬 GAPO - Genetic Algorithm for Protein Optimization")
st.markdown("A web interface for sequence-based protein engineering using ESM and Genetic Algorithms.")

# --- Sidebar for User Input ---
st.sidebar.header("Configuration")

initial_seq = st.sidebar.text_input("Initial Amino Acid Sequence", "MTSGTGSWIGSSGT")
apt_function = st.sidebar.selectbox(
    "Aptitude Function", 
    ["esm", "esm_penalty", "esm_shannon_penalty"], 
    index=0,
    help="Select the function to evaluate the fitness of each sequence."
)
direction = st.sidebar.selectbox(
    "Optimization Direction", 
    ["up", "down"], 
    index=0,
    help="'up' to maximize the score, 'down' to minimize it."
)

st.sidebar.subheader("GA Parameters")
residues_to_mute = st.sidebar.text_area(
    "Residues to Mutate (comma or space separated)", 
    "3, 4, 5, 6, 7, 8, 9",
    help="Provide the indices of the amino acids in the sequence that are allowed to mutate."
)
pop_size = st.sidebar.number_input("Population Size", min_value=5, max_value=500, value=50, step=10)
cycles = st.sidebar.number_input("Number of Cycles (Generations)", min_value=1, max_value=1000, value=20, step=5)
mutation_rate = st.sidebar.slider("Mutation Rate", min_value=0.0, max_value=1.0, value=0.9, step=0.05)
cpus = st.sidebar.number_input("Number of CPUs", min_value=1, max_value=16, value=4, step=1)

run_button = st.sidebar.button("🚀 Run GAPO Optimization")

# --- Main Area for Results and Debugging ---
st.header("Results & Debug")

if run_button:
    # Clean up old log file for the new run
    if os.path.exists("debug.log"):
        os.remove("debug.log")

    # --- 1. Build the Command ---
    log_message("Building command...")
    
    # --- Input Validation for Residues ---
    residue_list = residues_to_mute.replace(',', ' ').split()
    validated_residues = []
    try:
        for res in residue_list:
            if res.strip():
                validated_residues.append(str(int(res.strip())))
        
        if not validated_residues:
            raise ValueError("Residue input cannot be empty.")
            
    except ValueError:
        st.error("Invalid input for 'Residues to Mutate'. Please provide a comma or space-separated list of valid numbers only.")
        st.stop()
    
    # Define a unique output filename for this run
    output_filename = f"results/run_{int(time.time())}.csv"
    os.makedirs("results", exist_ok=True)
    
    command = [
        "python", "GA_main.py",
        "sequence",
        "--seq", initial_seq,
        "--residues_to_mute"
    ]
    command.extend(validated_residues)
    command.extend([
        "--apt_function", apt_function,
        "--direction", direction,
        "--pop_size", str(pop_size),
        "--cycles", str(cycles),
        "--mutation_rate", str(mutation_rate),
        "--cpus", str(cpus),
        "--output_file", output_filename
    ])
    
    # --- CORRECTED PART FOR DISPLAYING THE COMMAND ---
    # Explicitly convert all items in the list to strings before joining
    command_for_display = [str(item) for item in command]
    command_str = " ".join(command_for_display)
    
    st.info(f"Running GAPO with command:\n`{command_str}`")
    log_message(f"Executing command: {command_str}")
    
    # --- 2. Execute the Subprocess and Wait for It to Complete ---
    with st.spinner("GAPO is running... This may take a while. Please wait."):
        # We still pass the original 'command' list to subprocess.run
        result = subprocess.run(command, capture_output=True, text=True)

    # --- 3. Log and Display the Debugging Information ---
    log_message("Process finished.")
    log_message(f"Return Code: {result.returncode}")
    log_message(f"--- STDOUT ---\n{result.stdout}")
    log_message(f"--- STDERR ---\n{result.stderr}")
    
    st.subheader("Debug Output")
    
    if result.returncode == 0:
        st.success("✅ Process finished successfully!")
    else:
        st.error(f"❌ Process failed with return code {result.returncode}. Check the error logs below.")

    with st.expander("Show Standard Output (stdout)"):
        st.text(result.stdout if result.stdout else "No standard output.")
    
    with st.expander("Show Standard Error (stderr) - LOOK FOR ERRORS HERE!"):
        st.code(result.stderr if result.stderr else "No standard error.", language=None)
        
    # --- 4. Display the Final Results if Successful ---
    if os.path.exists(output_filename):
        st.subheader("Results Data")
        results_df = pd.read_csv(output_filename, header=None, names=['sequence', 'score', 'population'])
        st.dataframe(results_df)

        st.subheader("Best Score Convergence")
        try:
            if direction == 'down':
                best_scores = results_df.groupby('population')['score'].min()
            else:
                best_scores = results_df.groupby('population')['score'].max()
            
            st.line_chart(best_scores)
        except Exception as e:
            st.warning(f"Could not generate convergence plot. Error: {e}")
            
    else:
        st.warning(f"Output file '{output_filename}' was not created. Please check the error logs above for clues.")


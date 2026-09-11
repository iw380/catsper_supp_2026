import pandas as pd
import re
import os
import matplotlib.pyplot as plt
import seaborn as sns

# --- PLOT CONFIGURATION - TWEAK ANY PARAMETER HERE ---
PLOT_CONFIGURATION = {
    'figure': {'figsize': (14, 5), 'dpi': 300},
    # Configuration for the filled area plot
    'fill_plot': {'color': 'darkblue', 'alpha': 0.8, 'edgecolor': 'none'},
    'font_sizes': {'title': 20, 'axis_labels': 20, 'axis_ticks': 20},
    'cutoff_line': {'y_value': 0.6, 'color': 'red', 'style': '--', 'width': 1.5}
}
# -----------------------------------------------------------

def create_profile_plots(input_folder, output_folder):
    """
    Reads profile data and generates clean, continuous profile plots for the
    pDP and Sbind scores of each subunit, with a horizontal cutoff line.
    """
    subunits = ['catsper1', 'catsper2', 'catsper3', 'catsper4']
    
    plt.style.use('seaborn-v0_8-white') 
    
    print(f"Reading data from '{input_folder}/' and generating plots...")
    
    for subunit in subunits:
        file_path = os.path.join(input_folder, f'{subunit}_fuzdrop.tsv')
        
        if not os.path.exists(file_path):
            print(f"  - Warning: Data file for {subunit} not found. Skipping.")
            continue
            
        try:
            df = pd.read_csv(file_path, sep='\t', comment='#')
            print(f"  - Successfully loaded {subunit} data.")

            # --- Generate the pDP Profile Plot ---
            fig_pdp, ax_pdp = plt.subplots(figsize=PLOT_CONFIGURATION['figure']['figsize'])
            
            # --- REPLACEMENT FOR ax.bar() ---
            # Plot a line connecting the tops of the bars (linewidth=0 makes it invisible)
            ax_pdp.plot(df['position'], df['pDP'], lw=0) 
            # Fill the area between the line and the x-axis (y=0)
            ax_pdp.fill_between(df['position'], 0, df['pDP'], 
                                color=PLOT_CONFIGURATION['fill_plot']['color'],
                                alpha=PLOT_CONFIGURATION['fill_plot']['alpha'])

            # Add cutoff line
            ax_pdp.axhline(y=PLOT_CONFIGURATION['cutoff_line']['y_value'],
                         color=PLOT_CONFIGURATION['cutoff_line']['color'],
                         linestyle=PLOT_CONFIGURATION['cutoff_line']['style'],
                         linewidth=PLOT_CONFIGURATION['cutoff_line']['width'])
            
            # Styling
            ax_pdp.set_title(f'{subunit.upper()} pDP Profile', fontsize=PLOT_CONFIGURATION['font_sizes']['title'], pad=20)
            ax_pdp.set_xlabel('Residue Position', fontsize=PLOT_CONFIGURATION['font_sizes']['axis_labels'])
            ax_pdp.set_ylabel('pDP Score', fontsize=PLOT_CONFIGURATION['font_sizes']['axis_labels'])
            ax_pdp.tick_params(axis='both', which='major', labelsize=PLOT_CONFIGURATION['font_sizes']['axis_ticks'])
            ax_pdp.set_xlim(0, len(df) + 1)
            
            pdp_output_path = os.path.join(output_folder, f'{subunit}_pdp_profile.png')
            fig_pdp.savefig(pdp_output_path, dpi=PLOT_CONFIGURATION['figure']['dpi'], bbox_inches='tight')
            print(f"    - Saved pDP plot to '{pdp_output_path}'")
            plt.close(fig_pdp)

            # --- Generate the Sbind Profile Plot ---
            fig_sbind, ax_sbind = plt.subplots(figsize=PLOT_CONFIGURATION['figure']['figsize'])

            # --- REPLACEMENT FOR ax.bar() ---
            ax_sbind.plot(df['position'], df['Sbind'], lw=0)
            ax_sbind.fill_between(df['position'], 0, df['Sbind'],
                                  color=PLOT_CONFIGURATION['fill_plot']['color'],
                                  alpha=PLOT_CONFIGURATION['fill_plot']['alpha'])

            # Add cutoff line
            ax_sbind.axhline(y=PLOT_CONFIGURATION['cutoff_line']['y_value'],
                           color=PLOT_CONFIGURATION['cutoff_line']['color'],
                           linestyle=PLOT_CONFIGURATION['cutoff_line']['style'],
                           linewidth=PLOT_CONFIGURATION['cutoff_line']['width'])

            # Styling
            ax_sbind.set_title(f'{subunit.upper()} Sbind Profile', fontsize=PLOT_CONFIGURATION['font_sizes']['title'], pad=20)
            ax_sbind.set_xlabel('Residue Position', fontsize=PLOT_CONFIGURATION['font_sizes']['axis_labels'])
            ax_sbind.set_ylabel('Sbind Score', fontsize=PLOT_CONFIGURATION['font_sizes']['axis_labels'])
            ax_sbind.tick_params(axis='both', which='major', labelsize=PLOT_CONFIGURATION['font_sizes']['axis_ticks'])
            ax_sbind.set_xlim(0, len(df) + 1)
            
            sbind_output_path = os.path.join(output_folder, f'{subunit}_sbind_profile.png')
            fig_sbind.savefig(sbind_output_path, dpi=PLOT_CONFIGURATION['figure']['dpi'], bbox_inches='tight')
            print(f"    - Saved Sbind plot to '{sbind_output_path}'")
            plt.close(fig_sbind)

        except Exception as e:
            print(f"  - ERROR: Could not process file '{file_path}'. Error: {e}")
            continue

    print("\n--- All profile plots generated successfully! ---")

if __name__ == '__main__':
    INPUT_DATA_FOLDER = 'profile_data'
    OUTPUT_FOLDER = 'finalplots'
    
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    try:
        import pandas, matplotlib, seaborn
    except ImportError as e:
        print(f"Error: Required library '{e.name}' is not installed.")
        print("Please run: pip install pandas matplotlib seaborn")
    else:
        create_profile_plots(INPUT_DATA_FOLDER, OUTPUT_FOLDER)
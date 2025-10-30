import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import pandas as pd
import sys
import numpy as np

colors = {
    "All"      : "blue",
    "Positron" : "orange",
    "Electron" : "gray",
    "Photon"   : "yellow",
    "Muon"     : "cyan",
    "Proton"   : "green"
}

from matplotlib.ticker import ScalarFormatter

def plot_particles(df, initial_energy, output_file, show_zero_markers=True):
    plt.figure(figsize=(15,6))
    
    # Add two empty spaces before and after the x-axis
    x = list(range(-2, len(df['Tree']) + 2))
    lines = []
    
    ecal_front_row = df[df['Tree'] == 'ECALFront']
    if not ecal_front_row.empty:
        total_particles_at_ecal_front = ecal_front_row['All'].values[0] if 'All' in ecal_front_row else 1
    else:
        print("Warning: ECALFront not found in Tree column. Defaulting to 1 for denominator in percentage calculation.")
        total_particles_at_ecal_front = 1

    for particle_type, color in list(colors.items()):
        y_values = list(df[particle_type])
        y = [0, 0] + y_values + [0, 0]  # Adjusting y range to include the empty spaces
        if not show_zero_markers:
            y = [value if value > 0 else None for value in y]  # Do not show zero markers after adjusting y range
        particle_count_at_ecal_front = ecal_front_row[particle_type].values[0] if not ecal_front_row.empty else 0
        particle_percentage = 100 * particle_count_at_ecal_front / float(total_particles_at_ecal_front)
        label = "{} ({:.2f}%)".format(particle_type, particle_percentage)
        line, = plt.plot(x, y, color=color, marker='o', label=label)
        lines.append(line)
    
    x_ticks_labels = [' ', ' '] + list(df['Tree']) + [' ', ' ']
    plt.xticks(x, x_ticks_labels, rotation=45, ha="right", rotation_mode="anchor")

    # Set y-axis to scientific notation
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    plt.gca().yaxis.set_major_formatter(formatter)

    plt.ylabel('Number Of Particles')
    plt.title('Particles With PTotal Less Than %3.0 Of Nominal Energy [' + initial_energy + ' - e+]')
    plt.legend(handles=lines)

    plt.tight_layout()
    plt.savefig(output_file, format='pdf')
    print("Plot saved to", output_file)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plot_particle_production.py <csv_file> <initial_energy>")
        sys.exit(1)

    csv_file = sys.argv[1]
    initial_energy = sys.argv[2]
    output_file = "results_%s.pdf" % (initial_energy)

    df = pd.read_csv(csv_file)

    if 'Tree' not in df.columns:
        print("Error: 'Tree' column not found in", csv_file)
        print("Available columns are:")
        print(df.columns.tolist())
        sys.exit(1)

    print("DataFrame structure:")
    print(df.head())

    # Plot for all detectors
    plot_particles(df, initial_energy, output_file, show_zero_markers=False)

    # Plot for the last 7 detectors
    if len(df) >= 7:
        last_seven_df = df.tail(7)
        last_seven_output_file = "last_seven_" + output_file
        plot_particles(last_seven_df, initial_energy, last_seven_output_file, show_zero_markers=False)
    else:
        print("Warning: Less than 7 detectors found. Skipping the zoomed-in plot.")



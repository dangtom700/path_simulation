import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import subprocess
import os

def get_time_str():
    import datetime
    now = datetime.datetime.now()
    return now.strftime("%Y%m%d%H%M%S")

image_root = "bar_image/"
os.makedirs(image_root, exist_ok=True)
iteration = 1000

subprocess.run(["g++", "-o", "bar", "bar.cpp", "-lm"])

for i in range(iteration):
    subprocess.run(["./bar"])
    df = pd.read_csv("bar_data/joint_positions.csv")

    series_name = get_time_str()
    
    for bar in range(int(len(df.columns) / 2)-2):
        bar = bar + 1
        # Time and signals
        t = df["time"].to_numpy()
        x = df[f"x{bar}"].to_numpy()
        y = df[f"y{bar}"].to_numpy()

        num_samples = len(df)
        num_phases = 9  # You can change this to any number of phases you want
        samples_per_phase = num_samples // num_phases

        # Create phase labels
        phases = np.repeat(np.arange(1, num_phases + 1), samples_per_phase)

        # Handle leftover samples (if any) by appending the last phase label
        if len(phases) < num_samples:
            remaining = num_samples - len(phases)
            phases = np.append(phases, [num_phases] * remaining)

        # Add to DataFrame
        df["phase"] = phases
        # Define figure and GridSpec with 3 rows and 6 columns (3 left, 3 right)
        fig = plt.figure(figsize=(24, 12))
        gs = gridspec.GridSpec(3, 6, figure=fig, wspace=0.4, hspace=0.6)

        # Phase plots: occupy 3x3 slots on the left (columns 0–2)
        phase_axes = [fig.add_subplot(gs[i, j]) for i in range(3) for j in range(3)]
        phases = df["phase"].unique()

        # Plot each phase into its subplot
        for ax, phase in zip(phase_axes, phases):
            mask = df["phase"] == phase
            ax.plot(x[mask], y[mask], label=f"Phase {phase}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_title(f"Phase {phase}")
            ax.legend(loc = 'upper right')

        # Full trajectory plot: occupy all 3 rows of columns 3–5
        full_ax = fig.add_subplot(gs[:, 3:6])  # rows all, columns 3 to 5 (inclusive)
        for phase in phases:
            mask = df["phase"] == phase
            full_ax.plot(x[mask], y[mask], label=f"Phase {phase}")
        full_ax.set_xlabel("x")
        full_ax.set_ylabel("y")
        full_ax.set_title("Full Trajectory")
        full_ax.legend(loc = 'upper right')
        plt.title(f"Trajectory of Bar {bar}")

        plt.savefig(image_root + "trajectory_" + series_name + "_" + str(bar) + ".png")
        plt.close(fig)
    
    print(f"Finished iteration {i}")
import os

with open("report.md", "w") as f:
    f.write("# Path Simulation Using Computational Methods\n\n")
    for file in os.listdir("bar_image/"):
        if file.endswith(".png"):
            f.write(f"![{file}](bar_image/{file})\n")
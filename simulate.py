import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

days = 100
time = np.linspace(0, days, days)
tumor = np.zeros(days)
immune = np.zeros(days)
cytokine = np.zeros(days)

tumor[0] = 200
immune[0] = 50
cytokine[0] = 30

alpha = 0.05
beta = 0.02
gamma = 0.01
delta = 0.005

for t in range(1, days):
    tumor[t] = tumor[t-1] + alpha * tumor[t-1] - beta * tumor[t-1] * immune[t-1] * 0.0001
    immune[t] = immune[t-1] + gamma * cytokine[t-1] - delta * tumor[t-1] * 0.001
    cytokine[t] = cytokine[t-1] + 0.1 * (immune[t-1] - cytokine[t-1])

df = pd.DataFrame({
    "Day": time,
    "Tumor Cells": tumor,
    "Immune Cells": immune,
    "Cytokines": cytokine
})

df.to_csv("example_output.csv", index=False)

plt.plot(time, tumor, label="Tumor Cells")
plt.plot(time, immune, label="Immune Cells")
plt.plot(time, cytokine, label="Cytokines")
plt.xlabel("Day")
plt.ylabel("Concentration")
plt.title("Cancer-Immune Dynamics Simulation")
plt.legend()
plt.savefig("output_plot.png")
plt.show()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
data = pd.read_csv("stage2_average_score.csv")

sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))

sns.lineplot(data=data, x="Step", y="Value", marker="o", color="b")


plt.title("Average Total Score", fontsize=16)
plt.xlabel("Step", fontsize=14)
plt.ylabel("Score", fontsize=14)

plt.tight_layout()
plt.show()
plt.savefig("stage2_results_plot.png", dpi=300)




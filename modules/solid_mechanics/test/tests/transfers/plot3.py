
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("gold/coarse3_out_fine0.csv")
plt.plot(df["time"], df["end_disp_x"], "o-", label = 'Multigrid approach')
df = pd.read_csv("gold/fine3_out.csv")
plt.plot(df["time"], df["end_disp_x"], "o-", label = 'Mass-scaled fine-grid only')
plt.xlabel("Time")
plt.ylabel("End displacement (x)")
plt.grid(True)
plt.title("Harmonic oscillator")
plt.legend()
plt.savefig("result3.png", bbox_inches = 'tight')
plt.show()
                                                                                                                

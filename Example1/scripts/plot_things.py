import math
import matplotlib.pyplot as plt

in_file = snakemake.input["isoforms"]
out_file = snakemake.output["plot"]

# Pretend like we did something with `in_file` here...

# Make a plot!

x = range(100)
y = [math.sin(i/6) for i in x]

plt.figure()
plt.plot(x, y)
plt.savefig(out_file)


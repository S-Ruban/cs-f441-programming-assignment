import matplotlib.pyplot as plt

I = [6.25, 6.9, 8, 9, 10, 11, 12, 13, 14, 15]
f = [
    50.4413,
    57.4053,
    62.0347,
    65.1466,
    67.7966,
    70.0935,
    72.314,
    74.2312,
    76.087,
    77.8643,
]
plt.plot(I, f, "--")
plt.xticks(range(0, 17))
plt.yticks(range(0, 120, 20))
plt.title("f-I curve")
plt.xlabel("Current Density (μA/cm²)")
plt.ylabel("Frequency (Hz)")
plt.show()

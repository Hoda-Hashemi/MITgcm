import matplotlib.pyplot as plt
import numpy as np

ranks = np.array([8, 8, 12])
times = np.array([13:05:20, 9:37:16, 0.0])   # replace with your measured times

Tref = times[0]
speedup = Tref / times
efficiency = speedup / (ranks / ranks[0])

plt.figure()
plt.plot(ranks, times, marker='o')
plt.xlabel('MPI ranks')
plt.ylabel('Runtime (s)')
plt.title('MITgcm runtime vs MPI ranks')
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(ranks, speedup, marker='o', label='Measured')
plt.plot(ranks, ranks / ranks[0], linestyle='--', label='Ideal')
plt.xlabel('MPI ranks')
plt.ylabel('Speedup')
plt.title('Strong scaling speedup')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(ranks, efficiency, marker='o')
plt.xlabel('MPI ranks')
plt.ylabel('Parallel efficiency')
plt.title('Strong scaling efficiency')
plt.grid(True)
plt.tight_layout()
plt.show()

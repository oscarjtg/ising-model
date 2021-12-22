""" Module with plotting routines. """

import numpy as np
import matplotlib.pyplot as plt

def magnetisation_temperature(results, debug):
	""" Plots |magnetisation| / au against Temperature / K. """
	marker = ['s', 'o', 'v', 'd', '^']
	colour = ['b', 'm', 'g', 'r', 'y']
	fig, ax = plt.subplots()

	# keys index each Simulation class object, using ascending integers from 1 to n.
	keys = [int(key) for key in results.keys()]
	keys.sort()

	for key in keys:
		# Extract the appropriate Simulation object from the results dict,
		# then get the quantities for plotting.
		result = results[key]
		T = result.temperatures
		magnetisation = result.magnetisation
		error = result.magnetisation_error
		J = result.J
		mark = marker[key - 1] # key starts at 1, index starts at 0, hence '-1'.
		col = colour[key - 1]
		legend_label = f'J = {J:.2E} meV'

		# Plot
		ax.errorbar(x=T, y=np.abs(magnetisation), yerr=error, marker=mark, color=col, ms=4, label=legend_label)

	plt.xlabel('Temperature / K')
	plt.ylabel('|Magnetisation| / au')
	plt.title('Ising Model')
	ax.legend(loc='lower left')

	plt.savefig('outputs/Ising_model.png')
	plt.savefig('outputs/Ising_model.pdf')

	plt.show()
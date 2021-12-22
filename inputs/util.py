""" Module containing utility functions for reading inputs from inputs.txt. """

import numpy as np

debug = True
debug2 = False
debug3 = False

def read_inputs(filename, debug=False):
	"""Reads the input file and creates a dict with appropriate entries."""
	input_parameters = {}
	with open(filename, 'r') as f:
		file_text = f.read()
	file_lines = file_text.split('\n')
	for line in file_lines:
		if debug2: print(line)
		if line.startswith('#') or line == '':
			if debug2: print(f"skipped line {line}")
			continue

		pair = line.split()
		if debug: print(pair)
		# Expect two elements in pair (or possibly more for J).
		if len(pair) == 2:
			key = pair[0]
			value = pair[1]
		elif len(pair) < 2:
			if pair == []:
				if debug2: print("expect empty list: {pair}")
				continue
			raise ValueError(f"read_input: too few values in line - {pair}")
		elif len(pair) > 2:
			if pair[0] == "J": # User can specify multiple coupling constants J for separate calculations.
				key = pair[0]
				value = list(pair[1:])
			else:
				raise ValueError(f"read_inputs: too many values in line - {pair}")

		try:
			input_parameters[key] = float(value)
		except ValueError:
			if value == "True":
				input_parameters[key] = True
			elif value == "False":
				input_parameters[key] = False
			else:
				input_parameters[key] = value
		except TypeError:
			input_parameters[key] = value

	# The following code block modifies the input parameters as appropriate.
	for key in input_parameters.keys():
		value = input_parameters[key]
		if isinstance(value, str):
			if value.endswith('{e}'):
				input_parameters[key] = multiply_by_e(value, input_parameters)
			elif ':' in value:
				start, end, interval = [float(x) for x in value.split(':')]
				if debug2: print(start, end, interval)
				input_parameters[key] = np.arange(start, end, interval)
		if isinstance(value, list):
			for elem in value:
				if elem.endswith('{e}'):
					index = value.index(elem)
					value[index] = multiply_by_e(elem, input_parameters)
			input_parameters[key] = value

	if debug:
		print('input_parameters', input_parameters)
	return input_parameters



def multiply_by_e(value, input_parameters):
	if debug2: print(value)
	try:
		electronic_charge = input_parameters['e']
	except KeyError:
		raise KeyError(f"Have not specified electronic charge in {filename}")
	factor = float(value[:-3])
	if debug2:
		print(f'electronic charge = {electronic_charge}')
		print(f'factor = {factor}')

	return factor * electronic_charge
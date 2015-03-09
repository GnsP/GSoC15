#!/usr/bin/python
from numpy import array
from scipy.linalg import lu, solve
import random


def generate_input():
	m = random.randint(2,20)
	x = random.choice([True, False])
	n = m
	mat = []
	v = []
	for i in range(m):
		row = []
		for j in range(n):
			row.append(random.uniform(0,20))
		mat.append(row)
		v.append(random.uniform(0,20))

	return [m, n, array(mat), array(v)]

def verify_input(inp):
	try:
		pl, u = lu(inp[2], permute_l=True)
		inp.append(u)
		x = solve(inp[2], inp[3])
		inp.append(x)
		return inp
	except:
		return False


def print_testcase(inp):
	if inp != False:
		for i in range(3):
			print inp[0], inp[1]
			for i in range(inp[0]):
				print ' '.join([str(x) for x in inp[2][i]])
			print  ' '.join([str(x) for x in inp[3]])
			for i in range(inp[0]):
				print ' '.join([str(x) for x in inp[4][i]])
			print  ' '.join([str(x) for x in inp[5]])



for i in range(200):
	print_testcase(verify_input(generate_input()))


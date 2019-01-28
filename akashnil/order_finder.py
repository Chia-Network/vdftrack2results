'''
Copyright 2018 Chia Network Inc
Modifications copyright (C) 2019 Akashnil Dutta

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

from itertools import accumulate, repeat
from inkfish.classgroup import ClassGroup
from inkfish.create_discriminant import create_discriminant
from binascii import b2a_hex
import math
from inkfish.primes import odd_primes_below_n
from inkfish.primes import miller_rabin_test
import random
from sys import argv
import time
from multiprocessing import Pool
import subprocess
import os
import threading

# Store the output of giant steps in files in this directory temporarily.
hash_folder = 'hashes/'

def remove_all_hash_files():
    if not os.path.isdir(hash_folder):
        os.makedirs(hash_folder)
        return
    for the_file in os.listdir(hash_folder):
        file_path = os.path.join(hash_folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception as e:
            print(e)

def count_hash_files():
    if not os.path.isdir(hash_folder):
        os.makedirs(hash_folder)
    return len(os.listdir(hash_folder))

# The old baby_step_giant step algorithm only used for verification for testing
def compute_order(g: ClassGroup) -> int:
    d = g.discriminant()
    quadRoot = int((-d) ** 0.25)
    size = quadRoot
    order = 0
    while order == 0:
        babylist = list(accumulate(repeat(g, quadRoot - 1), ClassGroup.multiply))
        babyset = set(babylist)
        gSqrt = (g ** quadRoot).normalized()
        bigStep = gSqrt
        result = next(filter(lambda giant: giant[1] in babyset,
                      accumulate(repeat((1, bigStep), size),
                      lambda old, new: (old[0] + new[0], old[1] * new[1]))), None)
        if result is not None:
            order = (result[0] * quadRoot - babylist.index(result[1]) - 1)
            return order
        size *= 2

def find_group_order(discriminant):
	generator = ClassGroup.from_ab_discriminant(2, 1, d)
	order = compute_order(generator)
	return generator, order

# These are some custom parameters selected to optimize performance on target hardware
approx_time = 0.
prime_factor = 9.
# This is only limited by memory / num_approxers
max_primes_per_proc = 1<<27

baby_time = 0.
# If number of CPU cores is large, baby_factor must be increased
baby_factor = .07
# depends on available memory. Use test_memory function in calculator.cpp to find out.
baby_max = 1<<30

giant_time = 0.

def find_order_parameterized(abs_disc, num_baby_steps, middle, num_hashes_per_giant, num_giants):
	global baby_time, giant_time
	baby_time -= time.time()
	args = ['./a.out', 'baby', str(abs_disc), str(num_baby_steps)]
	baby_process = subprocess.Popen(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
	giant_step = num_baby_steps * 4
	giant_total = giant_step * num_hashes_per_giant
	indices = [x//2 if x%2 == 0 else -((x+1)//2) for x in range(1000000)]
	od = 0
	remove_all_hash_files()
	def run_giant(idx):
		x = 0
		while od == 0:
			y = idx + x * num_giants
			start = indices[y] * giant_total + middle
			if start < 0:
				return
			fname = hash_folder + 'hashes_' + str(y)
			f = open(fname, 'w+')
			args = ['./a.out', 'giant', str(abs_disc), str(start), str(giant_step), str(num_hashes_per_giant)]
			# print (' '.join(args))
			subprocess.call(args, stdout = f)
			baby_process.stdin.write((fname + '\n').encode('utf-8'))
			baby_process.stdin.flush()
			x += 1
			while count_hash_files() > num_giants * 3:
				time.sleep(0.1)
	threads = [threading.Thread(target=run_giant, args=[x]) for x in range(num_giants)]

	line = baby_process.stdout.readline().decode('utf-8')
	assert(line.strip() == 'populated')
	baby_time += time.time()

	giant_time -= time.time()

	for t in threads:
		t.start()

	line = baby_process.stdout.readline().decode('utf-8')
	od = int(line.strip())

	for t in threads:
		t.join()

	giant_time += time.time()

	baby_process.terminate()
	return od

# This is an universal constant. Significant upto 1.0471975 See writeup.
correction = 1.0471975223663745
# We update the coonstant a little bit as we calculate order of large primes.
learn_factor = 0.3

def find_approx_order_from_cpp(abs_disc, num_primes = 100):
	args = ['./a.out', 'approx', str(abs_disc), str(0), str(num_primes)]
	return int(correction * float(subprocess.check_output(args).strip()) * (abs_disc ** 0.5) / 2.)

def find_order_from_cpp(abs_disc):
	global correction, approx_time
	g = ClassGroup.from_ab_discriminant(2, 1, -abs_disc)
	num_baby_steps = int(abs_disc**.2 * baby_factor)
	num_baby_steps = min (baby_max, num_baby_steps)
	while True:
		fa = (g ** (num_baby_steps * 4))[0]
		if miller_rabin_test(fa) and fa >= num_baby_steps:
			break
		num_baby_steps += 1
	approx_time -= time.time()
	num_approxers = 16
	max_prime = int(abs_disc**.2 * prime_factor)
	primes_per_process = 1 + max_prime // num_approxers
	primes_per_process = min(primes_per_process, max_primes_per_proc)
	approx_partial = [0.] * num_approxers
	def approx(idx):
		val = 1.
		start = idx * primes_per_process
		while start < max_prime:
			args = ['./a.out', 'approx', str(abs_disc), str(start), str(start + primes_per_process)]
			# print (' '.join(args))
			val *= float(subprocess.check_output(args).strip())
			start += primes_per_process * num_approxers
		approx_partial[idx] = val

	threads = [threading.Thread(target=approx, args=[x]) for x in range(num_approxers)]
	for t in threads:
		t.start()
	for t in threads:
		t.join()

	# print (approx_partial)
	middle = (abs_disc ** 0.5) / 2.
	for v in approx_partial:
		middle *= v
	middle *= correction
	middle = int(middle)
	if middle % 2 == 1:
		middle += 1
	num_hashes_per_giant = min(int(abs_disc**.20), 1000000)
	num_giants = 16

	approx_time += time.time()
	od = find_order_parameterized(abs_disc, num_baby_steps, middle, num_hashes_per_giant, num_giants)
	correction_required = (correction * od) / middle
	correction = correction_required * learn_factor + correction * (1. - learn_factor)
	# print the true order, our approximate guess, and the updated value of correction.
	# Currently we manually set the initial value of this from the output back into the code.
	# Then for future order calculation the initial estimate is better. Can be automated.
	print (od, middle, correction)
	return od

def is_quadratic_residue(discriminant, prime):
	if prime == 2:
		return True
	def moduler_exponentiate(a, n, d):
		if n == 0:
			return 1
		if n % 2 == 0:
			x = moduler_exponentiate(a, n/2, d)
			return (x * x) % d
		else:
			return (moduler_exponentiate(a, n-1, d) * a) % d
	return moduler_exponentiate((discriminant % prime) + prime, (prime - 1) / 2, prime) == 1

def odd_primes_in_range(start, end):
	sieve = [True] * (end - start)
	if start == 0:
		sieve[0] = False
		sieve[1] = False
	if start == 1:
		sieve[0] = False
	for i in range(2, int(end ** 0.5) + 1):
		s = start - (start % i)
		if s < start:
			s += i
		s = max(s, i * i)
		sieve[s-start::i] = [False] * ((end - s + i - 1) // i)
	return [start + i for i in range(end - start) if sieve[i]]

small_odd = odd_primes_in_range(3, 100000)

def find_approx_order(abs_disc, num_primes = 5):
	ret = correction * (abs_disc ** 0.5) / 2.
	for i in range(len(small_odd)):
		p = small_odd[i]
		if is_quadratic_residue(-abs_disc, p):
			ret += ret / p
		else:
			ret -= ret / p
	return ret

length = int(argv[1])

if len(argv) > 2:
	num_disc = int(argv[2])
else:
	num_disc = 3

if len(argv) > 3:
	seed = int(argv[3])
else:
	seed = 1328946

subprocess.call(["g++ -std=c++11 -O3 calculator.cpp -lgmpxx -pthread -lgmp"], shell=True)

filename = 'entry.txt'
f = open(filename, "w+")

def disc_from_seed(i):
	sd = int.to_bytes(seed+i, 4, 'big')
	return create_discriminant(sd, length)

# choose some parameters to decide how much time should be spent looking for good 
# discriminant candidates (the ones which have smallest order)
num_candidates = 3
max_prime_in_search = 500
num_candidates += int((1 << length)  ** 0.2 / 500000)

exploration_time = time.time()

good_candidates = [x for x in range(num_candidates)]

# Round one: use only a few primes to narrow down the search for candidates
all_candidates = []
for i in range(num_candidates):
	d = disc_from_seed(good_candidates[i])
	all_candidates.append((find_approx_order_from_cpp(-d, max_prime_in_search // 10), good_candidates[i]))
all_candidates.sort()
good_candidates = [all_candidates[x][1] for x in range(len(all_candidates))]

# Round two: For top 10% promising candidate, use more primes to get a more accurate 
# approximation of the order, then we'll select top 3.
all_candidates = []
for i in range(max(num_candidates // 10, 3)):
	d = disc_from_seed(good_candidates[i])
	all_candidates.append((find_approx_order_from_cpp(-d, max_prime_in_search), good_candidates[i]))
all_candidates.sort()
good_candidates = [all_candidates[x][1] for x in range(len(all_candidates))]

# print time taken in phase one (refer to the pdf writeup). Choose parameters in a way such that
# this is small compared to the other phases in following sections.
print ('exploration_time: ', time.time() - exploration_time)

for i in range(num_disc):
	sd = int.to_bytes(seed+good_candidates[i], 4, 'big')
	d = create_discriminant(sd, length)
	od = find_order_from_cpp(-d)
	f.write('%s %i %i %i %i %i\n' % (b2a_hex(sd).decode('latin-1'), length, 2, 1, (-d+1) // 8, od))

# These are time taken in phase 2, 3, 4 (refer to the pdf writeup).
# Choose parameters in such a way that these times are approximately equal.
print ('approx time: ', approx_time)
print ('baby time: ', baby_time)
print ('giant time: ', giant_time)


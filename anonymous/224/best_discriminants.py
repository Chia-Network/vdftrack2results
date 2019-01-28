#!/usr/bin/env python3

from sys import argv

from hashlib import sha256
from binascii import b2a_hex

# --- full, but slow, create_discriminant
from inkfish.create_discriminant import create_discriminant

#### ----- modified code from inkfish.create_discriminant

m = 8 * 3 * 5 * 7 * 11 * 13
residues = [x for x in range(7, m, 8) if all([x % y != 0 for y in (3, 5, 7, 11, 13)])]

def entropy_from_seed(seed, byte_count):
    blob = bytearray()
    extra = 0
    while len(blob) < byte_count:
        extra_bits = extra.to_bytes(2, "big")
        more_entropy = sha256(seed + extra_bits).digest()
        blob.extend(more_entropy)
        extra += 1
    return bytes(blob[:byte_count])

def _create_discriminant(seed, length=2048):
    """
    Return a discriminant of the given length using the given seed.
    Generate a probable prime p where p % 8 == 7.
    Return -p.
    """
    extra = length % 8
    entropy = entropy_from_seed(seed, (length >> 3) + (2 if extra == 0 else 3))
    n = (int.from_bytes(entropy[:-2], 'big') >> (0 if extra == 0 else 8 - extra)) | (1 << (length - 1))
    n -= n % m
    n += residues[int.from_bytes(entropy[-2:], 'big') % len(residues)]

    # --- modification: cut off lengthy calculation which finds next prime
    return -n

#### ----------------------------------------------------------------------


best = [(0,0)]*16
length = int(argv[1])
#limit = 4295 # ~ 1 millionth (with original func, takes ~ 53 seconds)
#limit = 4294967 # ~ 1 thousandth (with redef above, takes ~ 15 seconds)
limit = 1<<32  # 4294967296

for i in range(limit):
    seed = int.to_bytes(i, 4, 'big')
    score = -_create_discriminant(seed, length)
    if score > best[0][0]:
        best[0] = (score,seed)
        best.sort()

for score,seed in best:
    d = create_discriminant(seed, length)
    print('%s %i %i' % (b2a_hex(seed).decode('latin-1'), length, d))


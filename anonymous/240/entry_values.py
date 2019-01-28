from itertools import accumulate, repeat
from inkfish.classgroup import ClassGroup
from inkfish.create_discriminant import create_discriminant
from binascii import b2a_hex
from hashlib import sha256

# results from initial run
best_seeds = [0xb1c21759, 0x1e7f8827, 0x1b03272b]

def find_seeds():
    results = []
    for i in range(3):
        seed = int.to_bytes(i, 4, 'big')
        r = bytearray(sha256(seed + b'\0\0').digest())
        r[0] |= 0x80
        results.append((r,seed))
    results.sort(reverse=True)
    _minR = results[-1][0]
    for i in range(3,0x100000000):
        seed = int.to_bytes(i, 4, 'big')
        r = bytearray(sha256(seed + b'\0\0').digest())
        r[0] |= 0x80
        if r > _minR:
            results[-1] = (r,seed)
            results.sort(reverse=True)
            _minR = results[-1][0]
    for r,seed in results:
        print('%s %s' % (b2a_hex(seed).decode('latin-1'),
                         b2a_hex(r).decode('latin-1')))

def print_entry(length):
    for i in range(3):
        seed = int.to_bytes(best_seeds[i], 4, 'big')
        d = create_discriminant(seed, length)
        print('%i' % (d))

if __name__ == '__main__':
    from sys import argv
    if len(argv)==2:
        print_entry(int(argv[1]))
    else:
        find_seeds()


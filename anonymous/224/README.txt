An off the shelf subexponential algorithm is available in pari gp, so it will be hard to beat a simple gp script. According to the documentation, pari uses Buchmann-McCurley's sub-exponential algorithm.

The discriminant and seed info was then given to a pari gp script which computes a generator and class size and outputs a line for inclusion in entry.txt

I overshot and expected 256 bit runs to finish in time. While these ran I searched all the seeds for the best scoring discriminants for smaller bitsizes. Getting nervous near the end, I stopped the larger computations to get some medium sized runs completed as backup. The 224 bit calculations were the longest that finished, and took two and a half days on three cores.

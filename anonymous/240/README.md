
# Additional information as requested

This is how I created the best entry.txt I could with the little compute power
I had access to.

On the wikipedia page for binary quadratic forms it mentions that the problem
we are trying to solve is equivalent to what is called the "class number" of
a quadratic field:  
"Binary quadratic forms are closely related to ideals in quadratic fields, this
allows the class number of a quadratic field to be calculated by counting the
number of reduced binary quadratic forms of a given discriminant."

Searching for "efficiently computing class number" resulted in some interesting
leads and articles that I didn't find when searching for binary quadratic forms.
For instance, apparently this would be weak to a quantum attack?

Jean-Francois Biasse, Fang Song. "Efficient quantum algorithms for computing
class groups and solving the principal ideal problem in arbitrary degree number
fields."  
https://pdfs.semanticscholar.org/959e/ebbc6ebde1cdddb9662b3bfc402f3e2f7a3c.pdf

Not having access to a quantum computer, I kept looking and stumbled upon a
comment here:  
https://math.stackexchange.com/questions/2971/algorithms-to-compute-the-class-number  
Just like factoring (also easy for quantum computers), it turns out this
problem also has a runtime described in terms of
L(N) = exp( (ln(N) ln(ln(N)))^(1/2) ).
The comment also mentions Henri Cohen's work and recommended Pari/GP.

Cohen, et al. discuss the subexponential algorithm here:  
H. Cohen, F. Diaz y Diaz, and M. Olivier. "Subexponential algorithms for class
group and unit computations" J. Symbolic Computation (1997) 24, 433-441  
https://core.ac.uk/download/pdf/82118906.pdf

This is 'quadclassunit' in Pari/GP and runs blazing fast compared to
baby-step giant-step.

However the runtimes didn't scale all that clearly with bitlength, so it was
risky to predict what size to attack for creating my entry.txt for the contest.

I thought (now I see likely incorrectly) that the way Chia was randomly
generating primes from the seeds was also contributing to the prediction noise.
So I decided to search for discriminants that started with as many 1's as
possible, because if I'm going to try to peg to one side of the range, I might
as well go for the top. In doing so, I noticed that the seeds which gave such
primes for one length, also gave it for another length. Looking at the code it
turned out that it took the result from the hash for the upper bits of the
prime. So one search would be universal, so I just searched all 32bits.

I then just entered stuff by hand and copy pasted. Once I got the class size
from quadclassunit, I made the entry.txt by hand with the seed, length, 
(2,1,(1-d)/8) generator, and class size.

# setup

```
sudo apt-get install pari-gp
```

# logs

Here are some snippets of saved console logs in case you are curious about the
raw numbers or run times.

```
$ time python3 entry_values.py | tee seeds.txt
b1c21759 ffffffffcead393239458be80a5aa5984506fa523fc1ddc24c281445627658e0
1e7f8827 fffffffed47688ecfd84420a1e34644ff756429d57ded0c914d476c1402cc448
1b03272b fffffffed394f5b687504cbb44fbe69f6318c12771dae5dea6320b233d61defd

real    107m54.002s
user    107m53.870s
sys     0m0.152s
```

```
$ python3 entry_values.py 240
-1766847064699124768416979591863805874605135527580828533601773796637238703
-1766847064297046971719834031638284082174997642667684668871886630236658287
-1766847064295631014183462482968335916105059566153134446723813027070016271
```

Note: It is necessary to give gp a large stack or it will crash with the larger
discriminants.
```
$ gp -s 5000000000
```

An automatic timer can be turned on with the pound command.
```
? #
   timer = 1 (on)
```

Here are the three runs with 240 bit primes.
```
? quadclassunit(-176684706429563101418346248296833591610505956615313444672381302
7070016271)
time = 231h, 45min, 53,543 ms.
%1 = [1513642557362553694458029240106055123, [1513642557362553694458029240106055
123], [Qfb(40093, 28085, 1101717921018401600144328488120330179897400772050691172
2269554709418)], 1]

? quadclassunit(-176684706429704697171983403163828408217499764266768466887188663
0236658287)
time = 261h, 50min, 12,910 ms.
%1 = [1024017792817320687110806446850261083, [1024017792817320687110806446850261
083], [Qfb(348936850638830590395137149298573937, -191282596524212052029805790309
000103, 1292093463851966445285967109578542552)], 1]

? quadclassunit(-176684706469912476841697959186380587460513552758082853360177379
6637238703)
time = 290h, 19min, 4,788 ms.
%1 = [1014689056833176952292251321171468791, [1014689056833176952292251321171468
791], [Qfb(53171, 23897, 8307381207326948752219158901768848971267869362908486456
911539171768)], 1]
```

As you can see, the run time still varies significantly despite the primes
differring by less than 1 part per billion.


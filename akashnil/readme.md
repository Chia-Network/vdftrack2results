## Chia network VDF competition track 2 submission

Author: Akashnil Dutta

The algorithm overview can be found in the pdf writeup DiscriminantBreak.pdf.

Put the content of this repository in 'vdf-competition/Entry2'. It needs access to inkfish.

Need to install g++ gmp etc. To run the code:

```
python3 order_finder.py [length]
```

It will generate an entry.txt file. There are lot of parameters in the python file specific to the hardware. To test for large lengths, those parameters must be set optimally for best results.

The target hardware used was m4.4xlarge instance type on AWS. 16 cCPUs, 64 GiB memory general purpose computing instance.
#!/bin/bash
# seed bitlength discriminant
gp -q -s 4G <<HEREDOC
{
seed="$1";
size="$2";
d=$3;

result = quadclassunit(d);
print(result);
period = result[2][1];
generator = Vec(result[3][1]);
print(seed," ",size," ",generator[1]," ",generator[2]," ",generator[3]," ",period);

quit;
}
HEREDOC


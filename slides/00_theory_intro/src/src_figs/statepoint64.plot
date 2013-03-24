#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output 'statepoint64.pdf'
set key bottom center
set xlabel 'Slab Position [cm]' 
set ylabel 'Fission Source'
set yrange [0.0:1.6]
set grid
unset key
set title  'Batch number 64'
set style line 1 lt 1 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "blue" lw 3
set style line 3 lt 1 lc rgb "orange" lw 3
plot '-' with lines ls 1

0.0606197268128
0.11334418301
0.16193793148
0.211334905788
0.265928432203
0.312489913256
0.366028206328
0.42052512936
0.472343406017
0.519182413704
0.570118898132
0.605652860608
0.649206165822
0.700246288634
0.748355276432
0.784210040576
0.828628318565
0.862233249102
0.898302718398
0.93878794085
0.968512927491
1.00981082937
1.04336084342
1.07857768523
1.10998221547
1.13629929937
1.16932008251
1.2103390445
1.2170359378
1.23129934069
1.27650040366
1.28525421836
1.31230988664
1.33403106916
1.34505440143
1.35961090063
1.37397092104
1.39408478004
1.40162099643
1.42002817798
1.42970203403
1.45052221743
1.44395523416
1.45996169947
1.4625644859
1.46320700255
1.4688384678
1.46295541023
1.45905998348
1.45916124226
1.47099828197
1.46534151348
1.47276258099
1.46248984024
1.45549884456
1.45868023974
1.44483636211
1.42737334523
1.4185923476
1.41765904665
1.39649620873
1.38842275288
1.38801444631
1.36430712634
1.33731517919
1.33782879828
1.31618313026
1.3053377466
1.27574619342
1.23668325027
1.22849012924
1.19578463266
1.1845958196
1.1386578449
1.14098629903
1.10833747161
1.06230934016
1.02584823696
0.998547452796
0.969741805196
0.931652181009
0.890339518721
0.870654435565
0.820616233672
0.779593224167
0.742170116276
0.698031714821
0.643252848336
0.606418499224
0.56679201219
0.511502089268
0.463452702264
0.415920454606
0.373216284794
0.313343911774
0.260974165097
0.212832119806
0.164598084519
0.114365393326
e

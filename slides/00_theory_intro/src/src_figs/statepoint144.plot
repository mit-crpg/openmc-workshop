#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output 'statepoint144.pdf'
set key bottom center
set xlabel 'Slab Position [cm]' 
set ylabel 'Fission Source'
set yrange [0.0:1.6]
set grid
unset key
set title  'Batch number 144'
set style line 1 lt 1 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "blue" lw 3
set style line 3 lt 1 lc rgb "orange" lw 3
plot '-' with lines ls 1

0.0535832692792
0.107072728861
0.150922117115
0.19787242538
0.255829077365
0.292867172868
0.350727415172
0.392502390726
0.43653485751
0.475361749791
0.51885166184
0.568389724584
0.612897932611
0.652449889466
0.687837957343
0.748238143103
0.780381065965
0.818524426213
0.84224208525
0.899299499661
0.940719921344
0.973364356359
1.01270266557
1.04341604762
1.08260857919
1.1015001258
1.13843882005
1.17089174968
1.21085795027
1.23488363165
1.26379523432
1.28611225083
1.31941319717
1.35085540975
1.35758463041
1.38907823676
1.40711349919
1.40965635828
1.43056900547
1.44632482902
1.46585130175
1.49360764024
1.48465091508
1.51337567715
1.50814140418
1.51334430979
1.52947762838
1.53611943405
1.53284648216
1.51467640672
1.54107542356
1.52726949426
1.53914301877
1.5270255477
1.51448640027
1.49790778772
1.49521047065
1.48180729607
1.47400464077
1.46905436755
1.43894447674
1.42237203006
1.4194795416
1.40352703976
1.3774116362
1.35888105839
1.32315084753
1.31342321915
1.27474369066
1.27247890634
1.22886733482
1.18756376633
1.17801299921
1.1337600774
1.12029583687
1.08477052229
1.06312710766
1.02648959426
0.981547577801
0.939434083436
0.902724717162
0.856943858063
0.821997306597
0.78494420182
0.745760929301
0.701025328957
0.664224692424
0.609751832441
0.578008521387
0.529535803485
0.484608714653
0.436933419323
0.39337168592
0.351621459851
0.299981873352
0.24675768591
0.204714274581
0.160441246784
0.107023365806
e

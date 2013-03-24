#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output 'statepoint128.pdf'
set key bottom center
set xlabel 'Slab Position [cm]' 
set ylabel 'Fission Source'
set yrange [0.0:1.6]
set grid
unset key
set title  'Batch number 128'
set style line 1 lt 1 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "blue" lw 3
set style line 3 lt 1 lc rgb "orange" lw 3
plot '-' with lines ls 1

0.0520335131851
0.107349673007
0.153004483984
0.200473574662
0.248340848067
0.298368962808
0.344580735587
0.386905714154
0.427697477887
0.480733995066
0.535435892586
0.566988974652
0.616168096266
0.653130852611
0.701695226225
0.746518404563
0.784848905311
0.831314225858
0.856639254172
0.895000257579
0.932799844784
0.977587721629
1.01209145149
1.05262100071
1.08353622804
1.10681868534
1.1350369779
1.17862341173
1.19921393197
1.23972779072
1.25962742032
1.28383519447
1.31749240134
1.32748112945
1.34579761403
1.36674494216
1.40972938172
1.42525159608
1.42590461059
1.45156507261
1.44731574072
1.485845939
1.48726773785
1.49073071683
1.50316548517
1.49166011142
1.52731975316
1.53415819722
1.53485110106
1.52486154577
1.52845450957
1.53390946521
1.51797248136
1.52241594614
1.5079222905
1.49261338274
1.48558192901
1.49010641588
1.4727423817
1.46180053073
1.44139144673
1.42264810986
1.39879367698
1.38845366186
1.36891854319
1.36410453488
1.34504574763
1.31607281617
1.27783131876
1.27134211665
1.24093341851
1.20864586601
1.19130461599
1.15180568199
1.12227503996
1.09249628382
1.06016644343
1.01636622617
0.980545549124
0.95677633102
0.912410768398
0.87853115916
0.829586455835
0.785054650721
0.746002407211
0.696442075482
0.668971482425
0.618280345817
0.577848423035
0.531695547418
0.490378301021
0.433526996472
0.389190240896
0.350483400174
0.29782539881
0.248715036627
0.206943870839
0.156971450616
0.105813423976
e

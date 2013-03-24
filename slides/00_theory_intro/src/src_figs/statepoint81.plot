#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output 'statepoint81.pdf'
set key bottom center
set xlabel 'Slab Position [cm]' 
set ylabel 'Fission Source'
set yrange [0.0:1.6]
set grid
unset key
set title  'Batch number 81'
set style line 1 lt 1 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "blue" lw 3
set style line 3 lt 1 lc rgb "orange" lw 3
plot '-' with lines ls 1

0.0581274047936
0.11111508083
0.163244709495
0.205412004619
0.259240234999
0.306654406433
0.357367859531
0.401817342913
0.440988022219
0.496754554213
0.546548740687
0.580154724289
0.637350582793
0.661319723982
0.714855512224
0.777222792737
0.800529828206
0.84314086919
0.887689784215
0.919612974979
0.960202636772
0.994522435731
1.0355360885
1.05197641994
1.09215688511
1.13948212831
1.14859224751
1.19034292829
1.21099578486
1.22395372664
1.27767036
1.28477794689
1.30669967843
1.32758194072
1.35020057089
1.36855907218
1.38570646233
1.4041652578
1.41935537462
1.42774852343
1.44644105924
1.45015121065
1.46300899586
1.46709489843
1.48276362231
1.48767111931
1.49793254126
1.49527884715
1.49550328504
1.51033566482
1.50113077148
1.49770973179
1.5013076479
1.49079791298
1.49636174793
1.47300457386
1.47192055657
1.46573555524
1.45855342881
1.44465265657
1.4086756685
1.41167334697
1.42515367946
1.3860644863
1.34555863818
1.35272004425
1.32288603936
1.30059141868
1.27705609755
1.25754685451
1.23128335417
1.21160243028
1.18080830099
1.14624613082
1.11037553764
1.09075077454
1.06187640427
1.01693856665
0.998885387814
0.957378956678
0.914012054205
0.877541647554
0.847775543575
0.801097481748
0.772626626049
0.722877811549
0.673535778042
0.64200666313
0.588948238806
0.548200074894
0.500658531523
0.450574543962
0.404877404695
0.35705639253
0.309488870064
0.256594809612
0.204127878549
0.157469097061
0.109728985375
e

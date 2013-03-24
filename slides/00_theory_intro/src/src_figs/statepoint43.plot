#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output 'statepoint43.pdf'
set key bottom center
set xlabel 'Slab Position [cm]' 
set ylabel 'Fission Source'
set yrange [0.0:1.6]
set grid
unset key
set title  'Batch number 43'
set style line 1 lt 1 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "blue" lw 3
set style line 3 lt 1 lc rgb "orange" lw 3
plot '-' with lines ls 1

0.0612868977045
0.127858854421
0.180838061213
0.240842034322
0.288218126362
0.350712461506
0.409367235925
0.461363737186
0.511527792484
0.559172860977
0.61435982874
0.659620782477
0.710299483858
0.750149809076
0.794278613597
0.825464377725
0.873552285155
0.91116525811
0.959302802818
0.973918093719
1.0213461469
1.0463276267
1.07257415314
1.10393897029
1.13854646551
1.16443672526
1.18168709268
1.22284132414
1.23568194862
1.25070392186
1.27082550224
1.28711908003
1.29757993839
1.3168665103
1.34200474311
1.33273388633
1.34939073172
1.35964119426
1.36523490674
1.36131967562
1.37834073596
1.37803702679
1.3674216262
1.39745729023
1.38632368443
1.3923597867
1.38575081355
1.39790195779
1.38498522464
1.39510401057
1.39805934371
1.3883160899
1.40056097631
1.39547609281
1.38659209555
1.38336351809
1.37335212124
1.3806701322
1.3601037483
1.35253676145
1.36389031888
1.34985300485
1.34323042091
1.32807641655
1.32147431291
1.29368836514
1.2865063907
1.27674164484
1.25427817865
1.2472467416
1.23134134801
1.21071124346
1.19476227017
1.17362341935
1.14420885268
1.13617097385
1.08819102624
1.07414012578
1.03367814001
1.00801081364
0.969901514693
0.933656095382
0.901249839205
0.855485073535
0.814669299892
0.771847364939
0.737338504664
0.694296097939
0.655877942103
0.595763753211
0.549013139693
0.489825166312
0.447126049675
0.386709444389
0.337935591316
0.291211808761
0.237523720528
0.183500818106
0.120431789779
e

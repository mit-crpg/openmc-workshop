#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output 'statepoint106.pdf'
set key bottom center
set xlabel 'Slab Position [cm]' 
set ylabel 'Fission Source'
set yrange [0.0:1.6]
set grid
unset key
set title  'Batch number 106'
set style line 1 lt 1 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "blue" lw 3
set style line 3 lt 1 lc rgb "orange" lw 3
plot '-' with lines ls 1

0.0544221976686
0.106776031346
0.157773354088
0.206967051942
0.247263614218
0.301895348825
0.347404893981
0.397209522478
0.436407686952
0.49090373099
0.52601809798
0.572638433387
0.616064350633
0.671716870079
0.708035740885
0.73999128375
0.787171399683
0.83155559006
0.876587235793
0.914378320146
0.941222855196
0.971463003202
1.02034649392
1.04808865602
1.06879195004
1.12291493537
1.15083887161
1.16703767461
1.20386944515
1.22693741246
1.24949157989
1.28555189342
1.315223766
1.32568441879
1.35211039755
1.37613642469
1.3952320581
1.41355866534
1.42234295435
1.44057249266
1.46108812262
1.46251547971
1.49071162567
1.49561151746
1.50324059671
1.50591254928
1.50412082488
1.50860666864
1.51262248588
1.50677160896
1.51053732047
1.51225103402
1.50920103506
1.52063648926
1.51640450296
1.4859688266
1.48887735084
1.46479636795
1.46158319529
1.45128321456
1.43404362618
1.41539602194
1.38873611923
1.40826780607
1.37861570696
1.36351468566
1.3312639778
1.31976290154
1.28786037255
1.26835948232
1.24681950204
1.21430506285
1.18456699549
1.1406730931
1.11190781843
1.08239931888
1.05284319236
1.01550773529
0.990047152415
0.943421922551
0.91308524265
0.874841583097
0.830975275155
0.804687932356
0.760999673875
0.714883499404
0.666023003988
0.639279258523
0.592553427976
0.540006781727
0.498099303625
0.447423123675
0.390151171506
0.348601640976
0.30589066263
0.259242166066
0.21190396368
0.15480086582
0.10690140954
e

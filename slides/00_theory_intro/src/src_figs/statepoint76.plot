#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output 'statepoint76.pdf'
set key bottom center
set xlabel 'Slab Position [cm]' 
set ylabel 'Fission Source'
set yrange [0.0:1.6]
set grid
unset key
set title  'Batch number 76'
set style line 1 lt 1 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "blue" lw 3
set style line 3 lt 1 lc rgb "orange" lw 3
plot '-' with lines ls 1

0.057977427388
0.109889762461
0.159800593208
0.208378668004
0.261380617783
0.316228475749
0.362561964358
0.40597510084
0.451080690627
0.511634627883
0.542669765583
0.587145558044
0.6279741581
0.687831461096
0.722114590564
0.772515567745
0.813979874312
0.841672433431
0.885729765158
0.927750659814
0.964169322732
1.00700031043
1.02514702383
1.06817755522
1.10500498585
1.12634147409
1.16262917933
1.17190331698
1.21488066428
1.24188578403
1.25775067912
1.29111578034
1.30613669656
1.33921185275
1.35100812088
1.37234423775
1.38789907529
1.40151419263
1.4092244278
1.42974544297
1.43234623578
1.458483677
1.46204641684
1.46906786837
1.4758158565
1.48785649814
1.4897869653
1.49246931375
1.49341277213
1.49478084652
1.50331932966
1.49455455633
1.48452516983
1.47708207763
1.47738747385
1.46969514195
1.4669955227
1.44456331353
1.44688135575
1.41367000258
1.41703219704
1.40595189001
1.40065058805
1.37375657212
1.36178127532
1.33855652887
1.31973534452
1.30573026701
1.2652075292
1.25940322081
1.22593247848
1.19510176767
1.17365472023
1.1376842728
1.12724476204
1.08656411671
1.06040881143
1.02239337153
1.00404257156
0.953074054886
0.918741214262
0.894366945194
0.850906015593
0.812040854538
0.773575792676
0.73089211997
0.687869498289
0.64348649939
0.601337417435
0.5538089277
0.514251895248
0.449952478943
0.408219292025
0.354378694433
0.306627047844
0.258508410718
0.216156055212
0.158041075105
0.110811145969
e
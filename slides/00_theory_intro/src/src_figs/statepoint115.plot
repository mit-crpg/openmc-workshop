#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output 'statepoint115.pdf'
set key bottom center
set xlabel 'Slab Position [cm]' 
set ylabel 'Fission Source'
set yrange [0.0:1.6]
set grid
unset key
set title  'Batch number 115'
set style line 1 lt 1 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "blue" lw 3
set style line 3 lt 1 lc rgb "orange" lw 3
plot '-' with lines ls 1

0.0549706044527
0.102831207023
0.156625575639
0.201245468058
0.246983900464
0.291609740463
0.350648960828
0.397239350148
0.437184633293
0.478519731424
0.525316948972
0.5646723395
0.615789773213
0.659269311496
0.713914250093
0.745390243698
0.786690936528
0.817006480782
0.870989017658
0.909566231008
0.940039225243
0.979491650738
1.01767935227
1.04393815182
1.08066775594
1.1220119842
1.14748335202
1.17417116483
1.2045598196
1.23238940657
1.25504016742
1.27622415904
1.29510548845
1.32671470835
1.35971350556
1.37773067216
1.38252952066
1.40739612419
1.43005683281
1.45713244419
1.46250679625
1.47204979726
1.48999693876
1.493612954
1.48626480102
1.49660484945
1.51118650903
1.50371812768
1.53172438021
1.52209892905
1.53239206425
1.52232276441
1.50954694931
1.5103709934
1.51474157425
1.49405921025
1.49277584853
1.4762904609
1.4525587342
1.45636423319
1.44017633374
1.42654758095
1.416760447
1.39000991163
1.38444524304
1.37277668402
1.33501658202
1.31169401573
1.27715951559
1.27204960884
1.23687058821
1.20124996446
1.18063832013
1.15089670809
1.10248009489
1.09100104234
1.06280376948
1.01518687142
0.990316623747
0.941637634441
0.927836660025
0.875363790716
0.839961373476
0.805081187174
0.757799037101
0.715132333054
0.679218419722
0.618896691574
0.582215179551
0.529710635084
0.492167219262
0.440836340702
0.399922791024
0.344401584046
0.296796363906
0.25529666472
0.208887965976
0.155684391044
0.107346725864
e

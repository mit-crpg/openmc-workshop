#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output 'statepoint51.pdf'
set key bottom center
set xlabel 'Slab Position [cm]' 
set ylabel 'Fission Source'
set yrange [0.0:1.6]
set grid
unset key
set title  'Batch number 51'
set style line 1 lt 1 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "blue" lw 3
set style line 3 lt 1 lc rgb "orange" lw 3
plot '-' with lines ls 1

0.0622203510939
0.118532742908
0.178585018888
0.226753737337
0.285461912762
0.32734285403
0.389070323142
0.442395046091
0.480899953118
0.534363408136
0.582817426726
0.635278064142
0.68267650649
0.729694403632
0.766132355105
0.812410244543
0.867408276123
0.887598323418
0.928704907449
0.975769147778
1.01415531345
1.03608368188
1.06245294755
1.10814758512
1.13567209195
1.15579821812
1.17163109748
1.20740970897
1.22573542738
1.23807550644
1.27076105319
1.28578323652
1.30486014615
1.31403276142
1.33479472181
1.34561520992
1.34371762516
1.37258692476
1.38226854594
1.3973081224
1.39171943521
1.38860893433
1.41015494725
1.41000589968
1.41263477966
1.42733543343
1.4157606951
1.42474987458
1.42611722767
1.42932567044
1.43746864877
1.42790688667
1.42274081436
1.41934286655
1.41346625789
1.41418882241
1.41036886815
1.39543014962
1.40896651678
1.38630167331
1.37627295304
1.36330422226
1.34714752311
1.3512153326
1.33652423784
1.31977379164
1.29272654874
1.3002956447
1.26917170208
1.2590380386
1.22861242758
1.20860301441
1.18185595162
1.15950751474
1.13935716954
1.11445672928
1.09062791919
1.04590364999
1.02478905455
0.987795518324
0.949858108329
0.916691283893
0.872085102286
0.835355416835
0.79647007331
0.755972791796
0.713450457851
0.672873638948
0.62921796957
0.581576357121
0.534058488795
0.480774581186
0.440298997564
0.387224404675
0.331824212951
0.2845515642
0.232473118992
0.174069165571
0.12059599188
e

rem ""
rem "starting fasta35_t - protein on win32"
rem ""
..\bin\fasta35_t -q -m 6 -Z 100000 ..\seq\mgstm1.aa:1-100 q > test_m1.ok2_t.html
..\bin\fasta35_t -S -q -z 11 -O test_m1.ok2_t_p25 -s P250 ..\seq\mgstm1.aa:100-218 q
rem "done"
rem "starting fastxy35_t"
..\bin\fastx35_t -m 9c -S -q ..\seq\mgtt2_x.seq q 1 > test_t2.xk1_t
..\bin\fasty35_t -S -q ..\seq\mgtt2_x.seq q > test_t2.yk2_t
..\bin\fastx35_t -m 9c -S -q -z 2 ..\seq\mgstm1.esq a > test_m1.xk2_tz2
..\bin\fasty35_t -S -q -z 2 ..\seq\mgstm1.esq a > test_m1.yk2_tz2
rem "done"
rem "starting fastxy35_t rev"
..\bin\fastx35_t -m 9c -q -m 5 ..\seq\mgstm1.rev q > test_m1.xk2r_t
..\bin\fasty35_t -q -m 5 -M 200-300 -z 2 ..\seq\mgstm1.rev q > test_m1.yk2r_tz2
..\bin\fasty35_t -q -m 5 -z 11 ..\seq\mgstm1.rev q > test_m1.yk2rz11_t
rem "done"
rem "starting ssearch35_t"
..\bin\ssearch35_t -m 9c -S -z 3 -q ..\seq\mgstm1.aa  q > test_m1.ss_tz3
..\bin\ssearch35_t -q -M 200-300 -z 2 -Z 100000 -s P250 ..\seq\mgstm1.aa q > test_m1.ss_t_p25
rem "starting ssearch35_t"
..\bin\ssearch35sse2_t -m 9c -S -z 3 -q ..\seq\mgstm1.aa  q > test_m1.ss_tz3sse2
..\bin\ssearch35sse2_t -q -M 200-300 -z 2 -Z 100000 -s P250 ..\seq\mgstm1.aa q > test_m1.ss_t_p25sse2
rem "done"
rem "starting prss\prfx35"
..\bin\ssearch35_t -q -k 1000 -A ..\seq\mgstm1.aa ..\seq\xurt8c.aa  > test_m1.rss
..\bin\fastx35_t -q -k 1000 -A ..\seq\mgstm1.esq ..\seq\xurt8c.aa > test_m1.rfx
rem "done"
rem "starting fasta35_t - DNA"
..\bin\fasta35_t -S -q -z 2 ..\seq\mgstm1.seq %M 4 > test_m1.ok4_tz2
..\bin\fasta35_t -S -q -r '+1/-2' ..\seq\mgstm1.rev %M 4 > test_m1.ok4r_t
rem "done"
rem "starting tfastxy35_t"
..\bin\tfastx35_t -m 9c -q -i -3 -m 6 ..\seq\mgstm1.aa %M > test_m1.tx2_t.html
..\bin\tfasty35_t -q -i -3 -N 5000 ..\seq\mgstm1.aa %M > test_m1.ty2_t
rem "done"
rem "starting fastf35_t"
..\bin\fastf35_t -q ..\seq\m1r.aa q > test_mf.ff_t
..\bin\fastf35 -q ..\seq\m1r.aa q > test_mf.ff_s
rem "done"
rem "starting tfastf35_t"
..\bin\tfastf35_t -q ..\seq\m1r.aa %m > test_mf.tf_tr
rem "done"
rem "starting fasts35_t"
..\bin\fasts35_t -q -V '*?@' ..\seq\ngts.aa q > test_m1.fs1_t
..\bin\fasts35_t -q ..\seq\ngt.aa q > test_m1.fs_t
..\bin\fasts35_t -q -n ..\seq\mgstm1.nts m > test_m1.nfs_t
rem "done"
rem "starting tfasts35_t"
..\bin\tfasts35_t -q ..\seq\n0.aa %m > test_m1.ts_r
rem "done"
rem "starting fasta35 - protein"
..\bin\fasta35 -q -z 2 ..\seq\mgstm1.aa q 1 > test_m1.ok1z2
..\bin\fasta35 -q -s P250 ..\seq\mgstm1.aa q > test_m1.ok2_p25 
rem "done"
rem "starting fastx3"
..\bin\fastx35 -m 9c -q ..\seq\mgstm1.esq q > test_m1.ok2x 
rem "done"
rem "starting fasty3"
..\bin\fasty35 -q ..\seq\mgstm1.esq q > test_m1.ok2y 
rem "done"
rem "starting fasta35 - DNA "
..\bin\fasta35 -m 9c -q ..\seq\mgstm1.seq M 4 > test_m1.ok4 
rem "done"
rem "starting ssearch3"
..\bin\ssearch35 -S -q -z 2 ..\seq\mgstm1.aa q > test_m1.ss_z2
..\bin\ssearch35 -q -s P250 ..\seq\mgstm1.aa q > test_m1.ss_p25 
..\bin\ssearch35sse2 -S -q -z 2 ..\seq\mgstm1.aa q > test_m1.ss_z2_sse2
..\bin\ssearch35sse2 -q -s P250 ..\seq\mgstm1.aa q > test_m1.ss_p25_sse2 
rem "done"
rem "starting tfastxy3"
..\bin\tfastx35 -q ..\seq\mgstm1.aa M > test_m1.tx2 
..\bin\tfasty35 -m 9c -q ..\seq\mgstm1.aa M > test_m1.ty2 
rem "done"
rem "starting fasts35"
..\bin\fasts35 -q -V '@?*' ..\seq\ngts.aa q > test_m1.fs1
..\bin\fasts35 -q ..\seq\ngt.aa q > test_m1.fs
rem "done"

#!/bin/csh -f
echo ""
setenv FASTLIBS /slib2/info/fast_libs_e.test
echo "FASTLIBS reset: " $FASTLIBS
echo "starting fasta35_t - protein" `date` "on" `hostname`
echo `uname -a`
echo ""
if (! -d results) mkdir results
../bin/fasta35_t -q -m 6 -Z 100000 ../seq/mgstm1.aa:1-100 +nq > results/test_m1.ok2_t.html
../bin/fasta35_t -S -q -z 11 -O results/test_m1.ok2_t_p25 -s P250 ../seq/mgstm1.aa:100-218 q
echo "done"
echo "starting fastxy35_t" `date`
../bin/fastx35_t -m 9c -S -q ../seq/mgtt2_x.seq q 1 > results/test_t2.xk1_t
../bin/fasty35_t -S -q ../seq/mgtt2_x.seq q > results/test_t2.yk2_t
../bin/fastx35_t -m 9c -S -q -z 2 ../seq/mgstm1.esq a > results/test_m1.xk2_tz2
../bin/fasty35_t -S -q -z 2 ../seq/mgstm1.esq a > results/test_m1.yk2_tz2
echo "done"
echo "starting fastxy35_t rev" `date`
../bin/fastx35_t -m 9c -q -m 5 ../seq/mgstm1.rev +ns > results/test_m1.xk2r_t
../bin/fasty35_t -q -m 5 -M 200-300 -z 2 ../seq/mgstm1.rev q > results/test_m1.yk2r_tz2
../bin/fasty35_t -q -m 5 -z 11 ../seq/mgstm1.rev q > results/test_m1.yk2rz11_t
echo "done"
echo "starting ssearch35_t" `date`
../bin/ssearch35_t -m 9c -S -z 3 -q ../seq/mgstm1.aa  q > results/test_m1.ss_tz3
../bin/ssearch35_t -q -M 200-300 -z 2 -Z 100000 -s P250 ../seq/mgstm1.aa q > results/test_m1.ss_t_p25
../bin/ssearch35_t -q -M 200-300 -s BL62 -f -11 -g -1  ../seq/mgstm1.aa q > results/test_m1.ss_t_bl62
../bin/ssearch35_t -q -M 200-300 -s ../data/blosum62.mat ../seq/mgstm1.aa q > results/test_m1.ss_t_bl62m
../bin/ssearch35_t -q -M 200-300 -S -s ../data/blosum62.mat ../seq/mgstm1.aa q > results/test_m1.ss_t_bl62mS
echo "done"
echo "starting ssearch35s_t" `date`
../bin/ssearch35s_t -m 9c -S -z 3 -q ../seq/mgstm1.aa  q > results/test_m1.sss_tz3
../bin/ssearch35s_t -q -M 200-300 -z 2 -Z 100000 -s P250 ../seq/mgstm1.aa q > results/test_m1.sss_t_p25
echo "done"
echo "starting prss35(ssearch/fastx)" `date`
../bin/ssearch35_t -q -k 1000 -A ../seq/mgstm1.aa ../seq/xurt8c.aa  > results/test_m1.rss
../bin/fastx35_t -q -k 1000 -A ../seq/mgstm1.esq ../seq/xurt8c.aa > results/test_m1.rfx
echo "done"
echo "starting ggsearch35/glsearch35" `date`
../bin/ggsearch35_t -q -m 9i -w 80 ../seq/hahu.aa q > results/test_h1.gg_t
../bin/glsearch35_t -q -m 9i -w 80 ../seq/hahu.aa q > results/test_h1.gl_t
../bin/ggsearch35_t -q ../seq/gtt1_drome.aa q > results/test_t1.gg_t
../bin/glsearch35_t -q ../seq/gtt1_drome.aa q > results/test_t1.gl_t
echo "done"
echo "starting fasta35_t - DNA" `date`
../bin/fasta35_t -S -q -z 2 ../seq/mgstm1.seq %M 4 > results/test_m1.ok4_tz2
../bin/fasta35_t -S -q -r '+1/-2' ../seq/mgstm1.rev %M 4 > results/test_m1.ok4r_t
echo "done"
#echo "starting tfasta35_t" `date`
#tfasta35_t -q ../seq/mgstm1.aa %RMB > results/test_m1.tk2_t
#echo "done"
echo "starting tfastxy35_t" `date`
../bin/tfastx35_t -m 9c -q -i -3 -m 6 ../seq/mgstm1.aa %p > results/test_m1.tx2_t.html
../bin/tfasty35_t -q -i -3 -N 5000 ../seq/mgstm1.aa %p > results/test_m1.ty2_t
echo "done"
echo "starting fastf35_t" `date`
../bin/fastf35_t -q ../seq/m1r.aa q > results/test_mf.ff_t
../bin/fastf35 -q ../seq/m1r.aa q > results/test_mf.ff_s
echo "done"
echo "starting tfastf35_t" `date`
../bin/tfastf35_t -q ../seq/m1r.aa %r > results/test_mf.tf_tr
echo "done"
echo "starting fasts35_t" `date`
../bin/fasts35_t -q -V '*?@' ../seq/ngts.aa q > results/test_m1.fs1_t
../bin/fasts35_t -q ../seq/ngt.aa q > results/test_m1.fs_t
../bin/fasts35_t -q -n ../seq/mgstm1.nts m > results/test_m1.nfs_t
echo "done"
echo "starting tfasts35_t" `date`
../bin/tfasts35_t -q ../seq/n0.aa %r > results/test_m1.ts_r
echo "done with threaded section: " `date`
echo "starting lalign35" `date`
../bin/lalign35 -k 1000 -q ../seq/mchu.aa ../seq/mchu.aa > results/test_mc.lal
../bin/lalign35 -z 3 -q ../seq/mchu.aa ../seq/mchu.aa > results/test_mc.lal_z3
../bin/lalign35 -s BL62 -f -11 -g -1  -q ../seq/mchu.aa ../seq/mchu.aa > results/test_mc.lal_bl62
../bin/lalign35 -k 1000 -q ../seq/mwkw.aa ../seq/mwkw.aa > results/test_mw.lal
../bin/lalign35 -z 3 -q ../seq/mwkw.aa ../seq/mwkw.aa > results/test_mw.lal_z3
../bin/lalign35 -s BL62 -f -11 -g -1  -q ../seq/mwkw.aa ../seq/mwkw.aa > results/test_mw.lal_bl62
echo "done"
echo "starting fasta35 - protein" `date`
../bin/fasta35 -q -z 2 ../seq/mgstm1.aa q 1 > results/test_m1.ok1z2
../bin/fasta35 -q -s P250 ../seq/mgstm1.aa q > results/test_m1.ok2_p25 
echo "done"
echo "starting fastx35" `date`
../bin/fastx35 -m 9c -q ../seq/mgstm1.esq q > results/test_m1.ok2x 
echo "done"
echo "starting fasty35" `date`
../bin/fasty35 -q ../seq/mgstm1.esq q > results/test_m1.ok2y 
echo "done"
echo "starting fasta35 - DNA " `date`
../bin/fasta35 -m 9c -q ../seq/mgstm1.seq %RMB 4 > results/test_m1.ok4 
echo "done"
echo "starting ssearch35" `date`
../bin/ssearch35 -S -q -z 2 ../seq/mgstm1.aa a > results/test_m1.ss_z2
../bin/ssearch35 -q -s P250 ../seq/mgstm1.aa a > results/test_m1.ss_p25 
echo "done"
echo "starting ssearch35s" `date`
../bin/ssearch35s -S -q -z 2 ../seq/mgstm1.aa a > results/test_m1.sss_z2
../bin/ssearch35s -q -s P250 ../seq/mgstm1.aa a > results/test_m1.sss_p25 
echo "done"
#echo "starting tfasta35" `date`
#tfasta35 -q ../seq/mgstm1.aa %RMB > results/test_m1.tk2 
#echo "done"
echo "starting tfastxy35" `date`
../bin/tfastx35 -q ../seq/mgstm1.aa %RMB > results/test_m1.tx2 
../bin/tfasty35 -m 9c -q ../seq/mgstm1.aa %RMB > results/test_m1.ty2 
echo "done"
echo "starting fasts35" `date`
../bin/fasts35 -q -V '@?*' ../seq/ngts.aa q > results/test_m1.fs1
../bin/fasts35 -q ../seq/ngt.aa q > results/test_m1.fs
echo "done" `date`

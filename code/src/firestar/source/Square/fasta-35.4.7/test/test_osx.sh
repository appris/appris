#!/bin/csh -f
echo ""
echo "starting fasta35_t - protein" `date` "on" `hostname`
echo `uname -a`
echo ""
fasta35_t -q -m 6 -Z 100000 mgstm1.aa:1-100 q > test_m1.ok2_t.html
fasta35_t -S -q -z 11 -O test_m1.ok2_t_p25 -s P250 mgstm1.aa:100-218 q
echo "done"
echo "starting fastxy35_t" `date`
fastx35_t -m 9 -S -q mgtt2_x.seq q 1 > test_t2.xk1_t
fasty35_t -S -q mgtt2_x.seq q > test_t2.yk2_t
fastx35_t -m 9 -S -q -z 2 mgstm1.esq a > test_m1.xk2_tz2
fasty35_t -S -q -z 2 mgstm1.esq a > test_m1.yk2_tz2
echo "done"
echo "starting fastxy35_t rev" `date`
fastx35_t -m 9 -q -m 5 mgstm1.rev q > test_m1.xk2r_t
fasty35_t -q -m 5 -M 200-300 -z 2 mgstm1.rev q > test_m1.yk2r_tz2
fasty35_t -q -m 5 -z 11 mgstm1.rev q > test_m1.yk2rz11_t
echo "done"
echo "starting ssearch35_t" `date`
ssearch35_t -m 9 -S -z 3 -q mgstm1.aa  q > test_m1.ss_tz3
ssearch35_t -q -M 200-300 -z 2 -Z 100000 -s P250 mgstm1.aa q > test_m1.ss_t_p25
echo "done"
echo "starting ssearch35_t" `date`
ssearch35s_t -m 9 -S -z 3 -q mgstm1.aa  q > test_m1.sss_tz3
ssearch35s_t -q -M 200-300 -z 2 -Z 100000 -s P250 mgstm1.aa q > test_m1.sss_t_p25
echo "done"
echo "starting prss35(ssearch/fastx)" `date`
ssearch35_t -q -k 1000 -A mgstm1.aa xurt8c.aa  > test_m1.rss
fastx35_t -q -k 1000 -A mgstm1.esq xurt8c.aa > test_m1.rfx
echo "done"
echo "starting ggsearch35/glsearch35" `date`
ggsearch35_t -q -m 9i -w 80 hahu.aa q > test_h1.gg_t
glsearch35_t -q -m 9i -w 80 hahu.aa q > test_h1.gl_t
ggsearch35_t -q gtt1_drome.aa q > test_t1.gg_t
glsearch35_t -q gtt1_drome.aa q > test_t1.gl_t
echo "done"
echo "starting fasta35_t - DNA" `date`
fasta35_t -q -z 2 mgstm1.seq %M 4 > test_m1.ok4_tz2
fasta35_t -q mgstm1.rev %M 4 > test_m1.ok4r_t
echo "done"
echo "starting tfastxy35_t" `date`
tfastx35_t -m 9 -q -i -3 -m 6 mgstm1.aa %m > test_m1.tx2_t.html
tfasty35_t -q -3 -N 5000 mgstm1.aa %m > test_m1.ty2_t
echo "done"
echo "starting fastf35_t" `date`
fastf35_t -q m1r.aa q > test_mf.ff_s
echo "done"
echo "starting tfastf35_t" `date`
tfastf35_t -q m1r.aa %m > test_mf.tf_r
echo "done"
echo "starting fasts35_t" `date`
fasts35_t -q n0.aa q > test_m1.fs_s
echo "done"
echo "starting tfasts35_t" `date`
tfasts35_t -q n0.aa %m > test_m1.ts_r
echo "done"
echo "done with threaded section: " `date`
echo "starting lalign35" `date`
lalign35 -k 1000 -q mchu.aa mchu.aa > test_mc.lal
lalign35 -z 3 -q mchu.aa mchu.aa > test_mc.lal_z3
lalign35 -s BL62 -f -11 -g -1  -q mchu.aa mchu.aa > test_mc.lal_bl62
lalign35 -k 1000 -q mwkw.aa mwkw.aa > test_mw.lal
lalign35 -z 3 -q mwkw.aa mwkw.aa > test_mw.lal_z3
lalign35 -s BL62 -f -11 -g -1  -q mwkw.aa mwkw.aa > test_mw.lal_bl62
echo "done"
echo "starting fasta35 - protein" `date`
fasta35 -q -z 2 mgstm1.aa q > test_m1.ok2z2
fasta35 -q -s P250 mgstm1.aa q > test_m1.ok2_p25 
echo "done"
echo "starting fastx35" `date`
fastx35 -m 9 -q mgstm1.esq q > test_m1.ok2x 
echo "done"
echo "starting fasty35" `date`
fasty35 -q mgstm1.esq q > test_m1.ok2y 
echo "done"
echo "starting fasta35 - DNA " `date`
fasta35 -m 9 -q mgstm1.seq %m 4 > test_m1.ok4 
echo "done"
echo "starting ssearch35" `date`
ssearch35 -S -q -z 2 mgstm1.aa q > test_m1.ss_z2
ssearch35 -q -s P250 mgstm1.aa q > test_m1.ss_p25 
echo "done"
echo "starting ssearch35" `date`
ssearchs35 -S -q -z 2 mgstm1.aa q > test_m1.sss_z2
ssearchs35 -q -s P250 mgstm1.aa q > test_m1.sss_p25 
echo "done"
echo "starting tfastxy35" `date`
tfastx35 -q mgstm1.aa %m > test_m1.tx2 
tfasty35 -m 9 -q mgstm1.aa %m > test_m1.ty2 
echo "done"
echo "starting fasts35" `date`
fasts35 -q -V '@?*' ngts.aa q > test_m1.fs1
fasts35 -q ngt.aa q > test_m1.fs
echo "done" `date`

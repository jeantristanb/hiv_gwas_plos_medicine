grep -w -f BPList.info /dataG/cancer/isabel_work/pGWAS/metal/HIV/meta_hiv_africa_loo/meta_hiv_africa_x_noJCS1.tbl > meta_hiv_africa_x_jcs1.tbl.tbl

fileI=/dataG/cancer/isabel_work/pGWAS/metal/HIV/meta_hiv_africa_x/meta_hiv_africa_x1.tbl
fileout=`basename $fileI`
grep -w -f BPList.info $fileI  > $fileout &
exit
for fileI in `ls ../*/*.tbl`
do
fileout=`basename $fileI`
grep -w -f BPList.info $fileI  > $fileout &
done
wait

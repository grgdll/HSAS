liste_w=`ls D_foQ_pa | cut -d_ -f2 | sort -n | awk 'BEGIN{p=-1};{if($1 != p)print; p=$1}'`
liste_s=`ls D_foQ_pa | cut -d_ -f3 | sort -n | awk 'BEGIN{p=-1};{if($1 != p)print; p=$1}'`
liste_c=`ls D_foQ_pa | cut -d_ -f4 | sort -n | awk 'BEGIN{p=-1};{if($1 != p)print; p=$1}'`
liste_n=`head -11 foQ.info | tail -1`
liste_a=`head -13 foQ.info | tail -1`

n_w=`echo $liste_w | awk '{print NF}'`
n_s=`echo $liste_s | awk '{print NF}'`
n_c=`echo $liste_c | awk '{print NF}'`
n_n=`echo $liste_n | awk '{print NF}'`
n_a=`echo $liste_a | awk '{print NF}'`
echo $n_w $liste_w
echo $n_s $liste_s
echo $n_c $liste_c
echo $n_n $liste_n
echo $n_a $liste_a

echo "      INTEGER N_W,N_S,N_C,N_N,N_A" >foq.h
echo "      PARAMETER(N_W=$n_w,N_S=$n_s,N_C=$n_c,N_N=$n_n,N_A=$n_a)" >>foq.h
echo "      real WAVE(N_W),SUN(N_S),CHL(N_C),NADIR(N_N),AZIM(N_A)" >>foq.h
echo "      data WAVE/`echo $liste_w | tr \" \" \",\"`/"  >>foq.h
echo "      data SUN/`echo $liste_s | tr \" \" \",\"`/"  >>foq.h
echo "      data CHL/`echo $liste_c | tr \" \" \",\"`/"  >>foq.h
echo "      data NADIR/`echo $liste_n | tr \" \" \",\"`/"  >>foq.h
echo "      data AZIM/`echo $liste_a | tr \" \" \",\"`/"  >>foq.h

test -s foq.dat && rm  foq.dat
for w in $liste_w
do
for s in $liste_s
do
for c in $liste_c
do
	cat D_foQ_pa/foQ_${w}_${s}_${c} >>foq.dat
done
done
done
awk -F="\000" '{OFS="";if(length($0) < 73)print;else print substr($0,1,72),"\n","     &",substr($0,73)}' foq.h >foq.tmp
mv foq.tmp foq.h

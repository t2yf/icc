V="0 1"
METRICS="L3CACHE ENERGY FLOPS_DP" 
N="64 100 128 " #1024 2000 2048 3000 4096 6000 7000 10000 50000 60000 70000 100000

for v in ${V}
do
    make clean && make V=${v}
    for m in ${METRICS}
    do
        #LIKWID_LOG="V${v}_${m}.log" 
        #rm -f ${LIKWID_LOG}
        for n in ${N}
        do
            LIKWID_OUT="V${v}${m}_${n}.txt"
            likwid-perfctr -C 3 -g ${m} -o ${LIKWID_OUT} -m -- ./prova ${n}
            #cat ${LIKWID_OUT} >> ${LIKWID_LOG}
            #rm -f ${LIKWID_OUT}
        done
    done
done
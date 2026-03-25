#!/bin/bash

METRICA="FLOPS_DP"
CPU=3

if [ "$1" = "-c" ]; then
    rm -f *~ ajustePol *.log
    exit
fi

LIKWID_HOME=/home/soft/likwid
CFLAGS="-I${LIKWID_HOME}/include -DLIKWID_PERFMON"
LFLAGS="-L${LIKWID_HOME}/lib -llikwid"
SCALING_FREQ="/sys/devices/system/cpu/cpufreq/policy${CPU}/scaling_governor"

likwid-topology -c -g >TOPOLOGIA.log

if [ -w ${SCALING_FREQ} ]; then
    echo "performance" > ${SCALING_FREQ}
fi

gcc ${CFLAGS} -O0 ajustePol.c -o ajustePol ${LFLAGS}
for k in $METRICA
do
    likwid-perfctr -C ${CPU} -g ${k} -m ./ajustePol >${k}_SemOtimiz.log
done

rm -f ajustePol
gcc ${CFLAGS} -O3 ajustePol.c -o ajustePol ${LFLAGS}
for k in $METRICA
do
    likwid-perfctr -C ${CPU} -g ${k} -m ./ajustePol >${k}_Otimiz.log
done

if [ -w ${SCALING_FREQ} ]; then
    echo "powersave" > ${SCALING_FREQ}
fi

echo ''
echo '-----------------------------------'
if [ -f FLOPS_DP_Otimiz.log -a -f TOPOLOGIA.log ]; then
    echo 'Configuração LIKWID: SUCESSO'
    echo ''
    ls FLOPS_*.log TOPOLOGIA.log
    echo ''
else
    echo 'Configuração LIKWID: FALHA'
fi
echo '-----------------------------------'
echo ''


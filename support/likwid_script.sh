#!/bin/bash

CMD_DIR=`dirname $0`

# Criando ambiente venv para rodar o script que gera os gráficos
VENV="broyden_venv"

if [ ! -d "${VENV}" ]; then
    echo " Criando '${VENV}' e instalando dependências (pandas, matplotlib)..."

    python3 -m venv ${VENV}
    ./${VENV}/bin/pip install --upgrade pip --quiet
    ./${VENV}/bin/pip install pandas matplotlib --quiet
    
    if [ $? -eq 0 ]; then
        echo "Ambiente Python broyden_venv criado"
    else
        echo " ERRO: Falha ao configurar o ambiente Python."
        exit 1
    fi
fi


PROG=${1:-broyden} # TODO colocar nome do arquivo a ser executado
CPU=${3:-3} # TODO qual core CPU vai rodar 

RESULTS_DIR="results/${PROG}"

METRICS="MEM L2CACHE FLOPS_DP FLOPS_AVX" #L3CACHE FLOPS_AVX ENERGY teste de banda de memória, cache miss l2, operações em ponto flutuante
                                            # se não rolar MEM ver L3
OPTIM="normal loop_unrolling gauss_seidel modified_newton one_vector_jacobian" #
MATRIX_SIZE="32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000" # 
INPUT_DIR="test"

rm -rf results/
mkdir -p ${RESULTS_DIR}

echo "performance" > /sys/devices/system/cpu/cpufreq/policy${CPU}/scaling_governor # manter freq alta de cpu

for optim in ${OPTIM}
do
    OPTIM_DIR="${RESULTS_DIR}/${optim}"
    OUTPUT_DIR="${OPTIM_DIR}/outputs"
    mkdir -p ${OPTIM_DIR}
    mkdir -p ${OUTPUT_DIR}

    OPTIM_MACRO=${optim^^}

    make clean > /dev/null
    make OPTIM=${OPTIM_MACRO}

    # Gera um .log para cada métrica com todos os resultados 
    for m in ${METRICS}
    do 
        LIKWID_LOG="${OPTIM_DIR}/${m}.log" # TODO especificar tipo
        rm -f ${LIKWID_LOG} #junta tudo em um log só, será q faz um para cada??
        INPUT_FILE="${INPUT_DIR}/input.in"
        
        # Cada tamanho
        for n in ${MATRIX_SIZE}
        do 
            LIKWID_OUT="${OPTIM_DIR}/${m}_${n}.txt"
            OUTPUT_FILE="${OUTPUT_DIR}/output_${n}.txt"
            
            #se input file existe
            if [ ! -f ${INPUT_FILE} ]; then
                echo "ERRO: ${INPUT_FILE} não encontrado!"
                continue
            fi
        
            echo "--->>  $m: ./${PROG} $n" >/dev/tty
            likwid-perfctr -O -C ${CPU} -g ${m} -o ${LIKWID_OUT} -m -- ./${PROG} ${n} -o ${OUTPUT_FILE} -v ${optim} < ${INPUT_FILE} 

            echo "N=${n}" >> ${LIKWID_LOG}
            cat ${LIKWID_OUT} >> ${LIKWID_LOG} #juntar todos os testes
            rm -f ${LIKWID_OUT}
        done
    done
done

# plotar os gráficos
source ${VENV}/bin/activate
#${CMD_DIR}/graphplot.py < ${LIKWID_LOG} > ${LIKWID_LOG%%.log}.csv
${CMD_DIR}/graphplot.py ${RESULTS_DIR}
deactivate

echo "powersave" > /sys/devices/system/cpu/cpufreq/policy${CPU}/scaling_governor #freq baixa de cpu
echo "" 








# for v in ${V}
# do
    # make clean && make V=${v}
    # for m in ${METRICS}
    # do
        # #LIKWID_LOG="V${v}_${m}.log" 
        # #rm -f ${LIKWID_LOG}
        # for n in ${N}
        # do
            # LIKWID_OUT="V${v}${m}_${n}.txt"
            # likwid-perfctr -C 3 -g ${m} -o ${LIKWID_OUT} -m -- ./prova ${n}
            # #cat ${LIKWID_OUT} >> ${LIKWID_LOG}
            # #rm -f ${LIKWID_OUT}
        # done
    # done
# done
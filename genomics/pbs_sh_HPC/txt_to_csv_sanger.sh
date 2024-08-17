#!/bin/bash

csv="sanger_stats_compiled_n746.csv"

echo "my_ID,sum,n_sum,ave,largest,N50,n_N50,N60,n_N60,N70,n_N70,N80,n_N80,N90,n_N90,N100,n_N100,N_count,Gaps" > $csv


for f in *.txt; do
    my_ID=$(basename "$f" .txt)
    
    
    while IFS= read -r line; do
        case $line in
            "sum ="*)
                sum=$(echo $line | awk -F ' ' '{print $3}')
                n_sum=$(echo $line | awk -F ' ' '{print $6}')
                ave=$(echo $line | awk -F ' ' '{print $9}')
                largest=$(echo $line | awk -F ' ' '{print $12}')
                ;;
            "N50 ="*)
                N50=$(echo $line | awk -F ' ' '{print $3}')
                n_N50=$(echo $line | awk -F ' ' '{print $6}')
                ;;
            "N60 ="*)
                N60=$(echo $line | awk -F ' ' '{print $3}')
                n_N60=$(echo $line | awk -F ' ' '{print $6}')
                ;;
            "N70 ="*)
                N70=$(echo $line | awk -F ' ' '{print $3}')
                n_N70=$(echo $line | awk -F ' ' '{print $6}')
                ;;
            "N80 ="*)
                N80=$(echo $line | awk -F ' ' '{print $3}')
                n_N80=$(echo $line | awk -F ' ' '{print $6}')
                ;;
            "N90 ="*)
                N90=$(echo $line | awk -F ' ' '{print $3}')
                n_N90=$(echo $line | awk -F ' ' '{print $6}')
                ;;
            "N100 ="*)
                N100=$(echo $line | awk -F ' ' '{print $3}')
                n_N100=$(echo $line | awk -F ' ' '{print $6}')
                ;;
            "N_count ="*)
                N_count=$(echo $line | awk -F ' ' '{print $3}')
                ;;
            "Gaps ="*)
                Gaps=$(echo $line | awk -F ' ' '{print $3}')
                ;;
        esac
    done < "$f"
    
    echo "$my_ID,$sum,$n_sum,$ave,$largest,$N50,$n_N50,$N60,$n_N60,$N70,$n_N70,$N80,$n_N80,$N90,$n_N90,$N100,$n_N100,$N_count,$Gaps" >> $csv
done


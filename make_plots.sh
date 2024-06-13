

# for i in "PO" "FS" "HS" "FC" "GP" "A"
# do
#     python test.py $i
#     Rscript test.R $i
# done



# for i in "PO" "FS" "HS" "FC" "GP" "A" "M" "U" "PO_inbreed";do
# #     for d in "no_deam" "deam"; do
for i in "U" ;do
    for d in "deam" "no_deam"; do
        Rscript test2.R $i $d
    done
done





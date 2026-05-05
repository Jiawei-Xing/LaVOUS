# conda activate simulate

# cd /grid/siepel/home/xing/gene_expression_evolution/SingleCellStochastics/lineage_simulation                                                                                                                             
                                                                                                                                                                                                                        
# # build (once)
# mkdir -p build && cd build && cmake .. && make && cd ..                                                                                                                                                                  
                                                                                                                                                                                                                        
# # simulate
# build/simulate -C 100 -c sim/color.txt -m 1 -mig 1e-2 -K 1000 -f 1 \                                                                                                                                                     
#     > sim/m1_mig1e-2_K1e3_f1.log                                                                                                                                                                                         

# reconstruct without sampling (proportion = 1.0 keeps all extant cells)                                                                                                                                                 
python cell_lineage.py \
    --cells sim/cloneCells.txt \                                                                                                                                                                                         
    --division sim/cellDivisionHistory.txt \                                                                                                                                                                             
    --migration sim/migrationCells.txt \

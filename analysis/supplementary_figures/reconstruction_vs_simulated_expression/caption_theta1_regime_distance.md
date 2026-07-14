# Supplementary figure: reconstruction diagnostics by theta on the original r=5 tree simulations

Panels A-C show the theta1=1, r=5 Gene 500 simulation on the original tree; panels D-F show theta1=3, r=5; panels G-I show theta1=5, r=5; panels J-L show theta1=7, r=5. In each row, the first panel compares LaVOUS posterior mean reconstructed expression with simulated latent expression at matched tree nodes, the second panel shows absolute reconstruction error versus distance from the root, and the third panel shows posterior reconstruction variance versus distance from the root. Black open circles and connecting lines in the error panels show median absolute error within fixed 4-unit root-distance bins. Expression values in the first column are centered by the simulated root expression for each scenario. Points are colored by reconstructed branch regime: gray for regime 0 and red for regime 1.

Summary across all nodes: theta1=1: r=0.74, MAE=0.99; theta1=3: r=0.76, MAE=0.78; theta1=5: r=0.82, MAE=0.98; theta1=7: r=0.81, MAE=0.97.

Source files:
- `analysis/lineage_simulation/tree.nwk`, `analysis/expression_simulation/simulation/sim_history_Base.tsv`, and `analysis/expression_simulation/reconst/history_Base_gene500.tsv`
- `analysis/lineage_simulation/tree.nwk`, `analysis/expression_simulation/simulation/sim_history_Theta1.tsv`, and `analysis/expression_simulation/reconst/history_Theta1_gene500.tsv`
- `analysis/lineage_simulation/tree.nwk`, `analysis/expression_simulation/simulation/sim_history_t5r5.tsv`, and `analysis/expression_simulation/reconst/history_t5r5_gene500.tsv`
- `analysis/lineage_simulation/tree.nwk`, `analysis/expression_simulation/simulation/sim_history_t7r5.tsv`, and `analysis/expression_simulation/reconst/history_t7r5_gene500.tsv`

Generated files:
- `lavous_theta1_regime_distance_reconstruction.pdf`
- `lavous_theta1_regime_distance_reconstruction.png`
- `theta1_regime_distance_source_data.tsv`
- `theta1_regime_distance_summary.tsv`
- `theta1_regime_distance_error_bins.tsv`

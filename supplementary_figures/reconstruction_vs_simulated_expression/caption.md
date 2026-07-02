# Supplementary figure: reconstructed versus simulated expression histories

Panels A-D compare LaVOUS posterior mean reconstructed expression with the simulated latent expression at matched tree nodes for the four Gene_500 OU examples shown in Fig. 2D. Values are centered on the simulated root for each scenario, matching the color centering used in Fig. 2D. Gray points are internal nodes and teal points are leaves; dashed lines indicate equality between simulated and reconstructed expression. Panel E shows residuals, defined as reconstructed minus simulated expression. Panel F summarizes mean absolute error (MAE) and root mean squared error (RMSE) across all nodes.

Summary across all nodes: $\theta_1=1$: r=0.65, MAE=1.21; $\theta_1=3$: r=0.76, MAE=0.70; $\theta_1=5$: r=0.80, MAE=0.89; $\theta_1=7$: r=0.87, MAE=0.89.

Source files:
- `expression_simulation/simulation/sim_history_r0.tsv` and `expression_simulation/reconst/history_r0_gene500.tsv`
- `expression_simulation/simulation/sim_history_t3r0.tsv` and `expression_simulation/reconst/history_t3r0_gene500.tsv`
- `expression_simulation/simulation/sim_history_t5r0.tsv` and `expression_simulation/reconst/history_t5r0_gene500.tsv`
- `expression_simulation/simulation/sim_history_t7r0.tsv` and `expression_simulation/reconst/history_t7r0_gene500.tsv`

Generated files:
- `lavous_reconstruction_vs_simulated_expression.pdf`
- `lavous_reconstruction_vs_simulated_expression.png`
- `source_data.tsv`
- `summary.tsv`

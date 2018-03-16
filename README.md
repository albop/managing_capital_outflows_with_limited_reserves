This code replicates the graphs from the paper *Managing Capital Outflows: The Role of Foreign. Exchange Intervention.* by Suman S. Basu, Atish R. Ghosh, Jonathan D. Ostry and Pablo E. Winant.

The code requires Python 3.6 and depends on two libraries:

- [dolo](https://github.com/EconForge/dolo)
- [backtothetrees](https://github.com/albop/backtothetrees)

An up-to-date version of the code can be found on [github](https://github.com/albop/managing_capital_outflows_with_limited_reserves)


To generate the results run:

- `python compute_simulations.py` : creates files with various simulation results
    - `precomputed_decision_rules.pickle`
    - `precomputed_moving_target.pickle`
    - `precomputed_simulations.pickle`

- `python compute_welfares.py`: create welfare comparison files:
    - `simple_rules_welfare.xlsx` (p=1.0)
    - `simple_rules_welfare_9.xlsx` (p=0.9)
    - `simple_rules_welfare_8.xlsx` (p=0.8)


The graphs are then created in the notebooke `gen_graphs.ipynb` which can be opened and run using Jupyter.

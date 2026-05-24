R-Code for the simulation study and the application of the paper "A robust nonparametric test for spatial isotropy in lattice data".

    Data: an ordner containing the following ordners
        Block: simulated grf with an outlier block for different block sizes and different outlier distributions (needed for section 4.2)
        Iso: simulated grf with isolated outliers for different amounts of outliers and differen outlier distributions (needed for section 4.2)
        Non: simulated grf without outliers (needed for section 4.1)

    Results: an ordner containing the following simulation results 
        res.quant: size corrected quantile (needed for section 4.1)
        res.size: size corrected power values (needed for section 4.1)
        res_block: power and type one error for data with outlier blocks (needed for section 4.2)
        res_iso: power and type one error for data with isolated outliers (needed for section 4.2)
        res_non: power and type one error for data without outliers (needed for section 4.1)

    Runtimes: an ordner containing the the simulated run times
     
    Application-Satellite-Data.R: R-Code for section 5 (Application to the satellite data)

    figure-tables.R: R-Code to obtain the graphs and tables of the paper (expect them of section 5)

    Parametercombinations.xlsx: an excel with all parameter combinations of the simulations for scenarios of the paper

    Parameters.R: R-Code, which defines the parameters for the simulations

    runtimes.R: R-Code, to determine the runtimes of the methods

    simulation_grf_block_outlier.R: R-Code of the simulations for section 4.2 (Contamination with block outliers)

    simulation_grf_isolated_outlier.R: R-Code of the simulations for section 4.2 (Contamination with isolated outliers)

    simulation_grf_non.R: R-Code of the simulations for section 4.1 (Gaussian Data)

    simulation-of-data.R: R-Code for the simulation of the grf (with and without outliers)

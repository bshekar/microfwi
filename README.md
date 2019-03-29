# Full waveform Inversion for microseismic events using sparsity constraints

We have developed an FWI algorithm to invert for a spatio-temporal source function that accounts for micro-seismic events with spatially localized or distributed locations and source signatures. The algorithm does not require assumptions to be made about the number or type of sources; however, it does require that the velocity model is close to the true subsurface velocity. We reformulate the conventional FWI algorithm based on the l2-norm data-misfit function by adding sparsity constraints using a sparsity promoting l1-norm as an additional regularization term to get more focused and less noise-sensitive event locations. The Orthant-Wise Limited-memory quasi-Newton algorithm is used to solve the optimization problem. It inherits the advantageous (fast convergence) properties of the limited memory Broyden- Fletcher-Goldfarb-Shanno method and can easily overcome the nondifferentiability of l1-norm at null positions. We determine the performance of the algorithm on noise-free and noisy synthetic data from the SEG/EAGE overthrust model.

For more details, please refer to: 

Shekar, B., and H. Sethi, 2019, Full waveform Inversion for microseismic events using sparsity constraints: Geophysics, 84, no.2, KS1-KS12

## Citation

    @article{microsparsefwi,
      title={Full waveform Inversion for microseismic events using sparsity constraints},
      author={Shekar, Bharath and Sethi, Harpreet Singh},
      journal={Geophysics},
      year={2019},
      Volume={84},
      Number={2},
      Pages={KS1--KS12}
    }

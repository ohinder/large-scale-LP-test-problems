mkdir -p large-problem-instances

julia --project generate-multicommodity-flow.jl \
    --output_file large-problem-instances/multicommodity-flow-instance.mps.gz \
    --num_commodities 10000 \
    --num_warehouses 100 \
    --num_stores 1000 \
    --seed 1

julia --project generate-heat-source-location.jl \
    --output_file large-problem-instances/heat-source-instance-easy.mps.gz \
    --ground_truth_file large-problem-instances/temperature_ground_truth-easy.hdf5 \
    --grid_size 250 \
    --num_source_locations 5 \
    --num_possible_source_locations 300 \
    --num_measurement_locations 700 \
    --seed 1 \
    --maximum_relative_measurement_error 0.0 \
    --pde_solve_tolerance 1e-12

julia --project generate-heat-source-location.jl \
    --output_file large-problem-instances/heat-source-instance-hard.mps.gz \
    --ground_truth_file large-problem-instances/temperature_ground_truth-hard.hdf5 \
    --grid_size 250 \
    --num_source_locations 10 \
    --num_possible_source_locations 500 \
    --num_measurement_locations 100 \
    --seed 2 \
    --maximum_relative_measurement_error 0.0 \
    --pde_solve_tolerance 1e-12

# Note: we calculate epsilon = 1/sqrt(num_treatment_samples) = 0.00141
julia --project design-matching-synthetic.jl \
    --output_file large-problem-instances/synthetic-design-match.mps.gz \
    --epsilon 0.00141 \
    --num_treatment_samples 500000 \
    --num_control_samples 5000000 \
    --num_covariates 15 \
    --num_edges_per_treatment 10 \
    --control_shift_magnitude 0.1 \
    --seed 1

julia --project generate-production-inventory.jl \
    --output_file large-problem-instances/production-inventory.mps.gz \
    --num_factories 200 \
    --num_stages 400 \
    --uncertainty_level 0.2

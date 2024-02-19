mkdir -p medium-problem-instances

julia --project generate-multicommodity-flow.jl \
    --output_file medium-problem-instances/multicommodity-flow-instance.mps.gz \
    --num_commodities 500 \
    --num_warehouses 100 \
    --num_stores 1000 \
    --seed 1

julia --project generate-heat-source-location.jl \
    --output_file medium-problem-instances/heat-source-instance1.mps.gz \
    --ground_truth_file medium-problem-instances/temperature_ground_truth1.hdf5 \
    --grid_size 150 \
    --num_source_locations 5 \
    --num_possible_source_locations 250 \
    --num_measurement_locations 50 \
    --seed 1

julia --project generate-heat-source-location.jl \
    --output_file medium-problem-instances/heat-source-instance2.mps.gz \
    --ground_truth_file medium-problem-instances/temperature_ground_truth2.hdf5 \
    --grid_size 150 \
    --num_source_locations 5 \
    --num_possible_source_locations 250 \
    --num_measurement_locations 25 \
    --seed 2

# Note: we calculate epsilon = 1/sqrt(num_treatment_samples) = 0.00316
julia --project design-matching-synthetic.jl \
    --output_file medium-problem-instances/synthetic-design-match.mps.gz \
    --epsilon 0.00316 \
    --num_treatment_samples 100000 \
    --num_control_samples 1000000 \
    --num_covariates 15 \
    --num_edges_per_treatment 10 \
    --control_shift_magnitude 0.1 \
    --seed 1

julia --project generate-production-inventory.jl \
    --output_file medium-problem-instances/production-inventory.mps.gz \
    --num_factories 100 \
    --num_stages 300 \
    --uncertainty_level 0.2

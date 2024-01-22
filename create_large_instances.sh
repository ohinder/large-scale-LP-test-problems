mkdir large-problem-instances

julia --project generate-multicommodity-flow.jl \
    --output_file large-problem-instances/multicommodity-flow-small-test-instance.mps.gz \
    --num_commodities 40000 \
    --num_warehouses 100 \
    --num_stores 2000 \
    --seed 1

julia --project generate-heat-source-location.jl \
    --output_file large-problem-instances/heat-source-instance1.mps.gz \
    --ground_truth_file large-problem-instances/temperature_ground_truth1.txt \
    --grid_size 500 \
    --num_source_locations 10 \
    --num_possible_source_locations 500 \
    --num_measurement_locations 100 \
    --seed 1

julia --project generate-heat-source-location.jl \
    --output_file large-problem-instances/heat-source-instance2.mps.gz \
    --ground_truth_file large-problem-instances/temperature_ground_truth2.txt \
    --grid_size 500 \
    --num_source_locations 10 \
    --num_possible_source_locations 500 \
    --num_measurement_locations 50 \
    --seed 2

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
    --num_factories 100 \
    --num_stages 200 \
    --uncertainty_level 0.2 \
    --seed 1

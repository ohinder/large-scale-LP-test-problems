mkdir small-problem-instances

julia generate-multicommodity-flow.jl \
    --output_file small-problem-instances/multicommodity-flow-small-test-instance.mps.gz \
    --num_commodities 100 \
    --num_warehouses 30 \
    --num_stores 100

julia generate-heat-source-location.jl \
    --output_file small-problem-instances/heat-source-instance1.mps.gz \
    --ground_truth_file small-problem-instances/temperature_ground_truth1.txt \
    --grid_size 50 \
    --num_source_locations 3 \
    --num_possible_source_locations 100 \
    --seed 1 \
    --num_measurement_locations 80

julia generate-heat-source-location.jl \
    --output_file small-problem-instances/heat-source-instance2.mps.gz \
    --ground_truth_file small-problem-instances/temperature_ground_truth2.txt \
    --grid_size 50 \
    --num_source_locations 3 \
    --num_possible_source_locations 100 \
    --seed 2 \
    --num_measurement_locations 40
    

# Note: we calculate epsilon = 1/sqrt(num_treatment_samples) = 0.0141
julia design-matching-synthetic.jl \
    --output_file small-problem-instances/synthetic-design-match.mps.gz \
    --epsilon 0.0141 \
    --num_treatment_samples 5000 \
    --num_control_samples 20000 \
    --num_covariates 8 \
    --num_edges_per_treatment 10 \
    --control_shift_magnitude 0.1
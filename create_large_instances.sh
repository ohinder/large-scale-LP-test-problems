mkdir large-problem-instances

julia-1.7 generate-multicommodity-flow.jl \
    --output_file large-problem-instances/multicommodity-flow-small-test-instance.mps.gz \
    --num_commodities 40000 \
    --num_warehouses 100 \
    --num_stores 2000

julia-1.7 generate-heat-source-location.jl \
    --output_file large-problem-instances/heat-source-instance1.mps.gz \
    --ground_truth_file large-problem-instances/temperature_ground_truth.txt \
    --grid_size 500 \
    --num_source_locations 10 \
    --num_possible_source_locations 500 \
    --num_measurement_locations 100

julia-1.7 generate-heat-source-location.jl \
    --output_file large-problem-instances/heat-source-instance2.mps.gz \
    --ground_truth_file large-problem-instances/temperature_ground_truth.txt \
    --grid_size 500 \
    --num_source_locations 10 \
    --num_possible_source_locations 500 \
    --num_measurement_locations 50

julia-1.7 design-matching-synthetic.jl \
    --output_file test_design_match.mps.gz \
    --epsilon 0.00141 \ : '1/sqrt(500000)'
    --num_treatment_samples 500000 \
    --num_control_samples 5000000 \
    --num_covariates 15 \
    --num_edges_per_treatment 10 \
    --control_shift_magnitude 0.1
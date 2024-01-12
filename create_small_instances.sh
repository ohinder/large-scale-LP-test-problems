mkdir small-problem-instances

julia-1.7 generate-multicommodity-flow.jl \
    --output_file small-problem-instances/multicommodity-flow-small-test-instance.mps.gz \
    --num_commodities 100 \
    --num_warehouses 30 \
    --num_stores 100

julia-1.7 generate-heat-source-location.jl \
    --output_file small-problem-instances/heat-source-instance1.mps.gz \
    --ground_truth_file small-problem-instances/temperature_ground_truth.txt \
    --grid_size 50 \
    --num_source_locations 3 \
    --num_possible_source_locations 100 \
    --num_measurement_locations 80

julia-1.7 generate-heat-source-location.jl \
    --output_file small-problem-instances/heat-source-instance2.mps.gz \
    --ground_truth_file small-problem-instances/temperature_ground_truth.txt \
    --grid_size 50 \
    --num_source_locations 3 \
    --num_possible_source_locations 100 \
    --num_measurement_locations 40

julia-1.7 design-matching-synthetic.jl \
    --output_file test_design_match.mps.gz \
    --epsilon 0.0141 \ : '1/sqrt(5000)'
    --num_treatment_samples 5000 \
    --num_control_samples 20000 \
    --num_covariates 8 \
    --num_edges_per_treatment 10 \
    --control_shift_magnitude 0.1
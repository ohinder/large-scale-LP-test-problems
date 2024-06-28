gfortran -O3 -o qap-gen newlp-56.f # this build the https://netlib.sandia.gov/lp/generators/qap/newlp.f file modified so 56 is the maximum input size

instance_name="tai50a"
file_in_tj_format="QAP-TJ-format/${instance_name}.txt"
julia convert-QAP-to-TJ-format.jl \
    --input_file "QAPLIB/${instance_name}.dat" \
    --output_file "$file_in_tj_format"

# When I read this mps file (in julia), I get:
# ERROR: LoadError: Duplicate row name C01A*** at line 62068
./qap-gen < "$file_in_tj_format" | gzip > "QAP-LP/${instance_name}.mps.gz"
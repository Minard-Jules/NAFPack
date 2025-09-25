#!/bin/bash
SUPP_FILE="$(dirname "$0")/testdrive.supp"
TEST_DIR=$(find ./build -type f \( -path "*/test/*" -o -path "./build/bin/*" \) -executable \
    -not -path "*/test-drive/test/*" \
    -not -path "./build/test/*" \
    ! -name "*.log")
    
declare -A errors
declare -A leaks
total_leaks=0

for exe in $TEST_DIR; do
    echo "Testing $exe"
    output=$(valgrind --leak-check=full --suppressions="$SUPP_FILE" "$exe" 2>&1)

    lost=$(echo "$output" | grep "definitely lost:" | awk '{print $4}')
    if [[ "$lost" == "bytes" ]]; then
        lost=0
    fi
    leaks["$exe"]=$lost
    total_leaks=$((total_leaks + $lost))

    err_sum=$(echo "$output" | grep "ERROR SUMMARY" | sed 's/^.*==[0-9]*== //')
    errors["$exe"]="$err_sum"
done

echo ""
echo "Valgrind Memory Leak Summary :"
echo "-----------------------------"
for exe in "${!leaks[@]}"; do
    echo "$exe : ${leaks[$exe]} bytes lost "
    echo "${errors[$exe]}"
    echo "-----------------------------"
done
echo "Total bytes lost : $total_leaks"
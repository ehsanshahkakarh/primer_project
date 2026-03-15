#!/bin/bash
# Verification script for PPDesign pipelines
# Tests all four main pipelines with sample data

set -e  # Exit on first error

echo "=========================================="
echo "PPDesign Pipeline Verification"
echo "=========================================="
echo ""

# Check PYTHONPATH
if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH=src
fi

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Test counter
TESTS_PASSED=0
TESTS_FAILED=0

# Helper function to run test
run_test() {
    local test_name="$1"
    local command="$2"
    local expected_file="$3"

    echo "Testing: $test_name"
    echo "Command: $command"

    if eval "$command" > /dev/null 2>&1; then
        if [ -f "$expected_file" ]; then
            echo -e "${GREEN}✓ PASS${NC}: $test_name"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}✗ FAIL${NC}: Output file not found: $expected_file"
            TESTS_FAILED=$((TESTS_FAILED + 1))
        fi
    else
        echo -e "${RED}✗ FAIL${NC}: Command failed"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    echo ""
}

# Clean up old results
echo "Cleaning up old test results..."
rm -rf results/test_*
echo ""

# Test 1: Unified pipeline (if test data available)
if [ -f "tests/test.fna" ]; then
    run_test "Unified Pipeline" \
        "python -m ppdesign.cli unified main --fasta-input tests/test.fna --conservation 0.3 --output-dir test_unified --threads 2" \
        "results/test_unified/summary.txt"
fi

# Test 2: Nucleotide pipeline (k-mer method)
if [ -f "tests/test.fna" ]; then
    run_test "Nucleotide Pipeline (k-mer)" \
        "python -m ppdesign.cli nucleotide main --fasta-input tests/test.fna --method kmer --conservation 0.7 --output-dir test_nucleotide --threads 2" \
        "results/test_nucleotide/selected_oligos.csv"
fi

# Test 3: Primer pipeline
if [ -f "resources/test_data/primer/test_targets.fna" ]; then
    run_test "Primer Pipeline" \
        "python -m ppdesign.cli primer main --fasta-input resources/test_data/primer/test_targets.fna --conservation 0.7 --amplicon-min 50 --amplicon-max 500 --output-dir test_primer" \
        "results/test_primer/primer_pairs.csv"
else
    echo "⚠ Skipping primer test: test data not found"
    echo ""
fi

# Test 3b: Primer pipeline with DEFAULT parameters (Bug 3 fix verification)
if [ -f "resources/test_data/primer/test_targets.fna" ]; then
    echo "Testing: Primer Pipeline (Default Parameters)"
    echo "Command: python -m ppdesign.cli primer main --fasta-input resources/test_data/primer/test_targets.fna --output-dir test_primer_default"

    if python -m ppdesign.cli primer main --fasta-input resources/test_data/primer/test_targets.fna --output-dir test_primer_default > /dev/null 2>&1; then
        if [ -f "results/test_primer_default/primer_pairs.csv" ]; then
            # Verify we got at least one pair
            pair_count=$(tail -n +2 results/test_primer_default/primer_pairs.csv | wc -l)
            if [ "$pair_count" -ge 1 ]; then
                echo -e "${GREEN}✓ PASS${NC}: Default parameters test ($pair_count pairs)"
                TESTS_PASSED=$((TESTS_PASSED + 1))
            else
                echo -e "${RED}✗ FAIL${NC}: No valid primer pairs with default parameters (Bug 3 not fixed)"
                TESTS_FAILED=$((TESTS_FAILED + 1))
            fi
        else
            echo -e "${RED}✗ FAIL${NC}: Output file not found"
            TESTS_FAILED=$((TESTS_FAILED + 1))
        fi
    else
        echo -e "${RED}✗ FAIL${NC}: Command failed with default parameters"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    echo ""
fi

# Test 4: gRNA pipeline
if [ -f "tests/test.fna" ]; then
    run_test "gRNA Pipeline" \
        "python -m ppdesign.cli grna main --fasta-input tests/test.fna --conservation 0.3 --output-dir test_grna" \
        "results/test_grna/guide_rnas.csv"
fi

# Verify output directories exist
echo "Verifying output directories..."
for dir in test_unified test_nucleotide test_primer test_grna; do
    if [ -d "results/$dir" ]; then
        echo -e "${GREEN}✓${NC} results/$dir exists"

        # List files in directory
        file_count=$(ls -1 results/$dir 2>/dev/null | wc -l)
        echo "   Files: $file_count"
    else
        echo -e "${RED}✗${NC} results/$dir not found"
    fi
done
echo ""

# Summary
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo -e "Tests passed: ${GREEN}$TESTS_PASSED${NC}"
if [ $TESTS_FAILED -gt 0 ]; then
    echo -e "Tests failed: ${RED}$TESTS_FAILED${NC}"
else
    echo -e "Tests failed: $TESTS_FAILED"
fi
echo ""

# Check specific output files
echo "Checking specific output files..."

if [ -f "results/test_primer/primer_pairs.csv" ]; then
    pair_count=$(tail -n +2 results/test_primer/primer_pairs.csv | wc -l)
    echo -e "${GREEN}✓${NC} Primer pairs generated: $pair_count"
fi

if [ -f "results/test_primer/primer_pairs.fasta" ]; then
    echo -e "${GREEN}✓${NC} FASTA file created"
fi

if [ -f "results/test_primer/summary.txt" ]; then
    echo -e "${GREEN}✓${NC} Summary file created"
fi

if [ -f "results/test_grna/guide_rnas.csv" ]; then
    grna_count=$(tail -n +2 results/test_grna/guide_rnas.csv | wc -l)
    echo -e "${GREEN}✓${NC} Guide RNAs generated: $grna_count"
fi

echo ""

# Exit with appropriate code
if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed.${NC}"
    exit 1
fi

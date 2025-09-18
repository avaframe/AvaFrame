#!/bin/bash

# Enhanced script to extract and match simulations with their parameters and timing
# Usage: ./extractProcessTimesAndParams <logfile>

if [ $# -eq 0 ]; then
    echo "Usage: $0 <logfile>"
    echo "Example: $0 runAna4ProbAna_20250917_14h34m46s.log"
    exit 1
fi

logFile="$1"

if [ ! -f "$logFile" ]; then
    echo "Error: File '$logFile' not found!"
    exit 1
fi

echo "Extracting simulation data with parameters from: $logFile"
echo "=================================================================="
echo

# Create temporary files to store the extracted data
processLines=$(mktemp)
timingLines=$(mktemp)
thicknessLines=$(mktemp)
musamosatLines=$(mktemp)

# Extract all relevant data with line numbers
grep -n "runs as process" "$logFile" > "$processLines"
grep -n "cpu time DFA" "$logFile" > "$timingLines"
grep -n "'relTh0':" "$logFile" > "$thicknessLines"
grep -n "Friction model parameter used: musamosat" "$logFile" > "$musamosatLines"

# Functions to extract values
getProcessId() {
    echo "$1" | cut -d',' -f1 | cut -d':' -f2
}

getSimName() {
    echo "$1" | sed 's/.*INFO - \(.*\) runs as process.*/\1/'
}

getCpuTime() {
    echo "$1" | sed 's/.*cpu time DFA = \([0-9.]*\) s.*/\1/'
}

getThickness() {
    echo "$1" | sed "s/.*'relTh0': '\([0-9.]*\)'.*/\1/"
}

getMusamosat() {
    echo "$1" | sed 's/.*musamosat with value \([0-9.]*\).*/\1/'
}

getLineNumber() {
    echo "$1" | cut -d':' -f1
}

# Create arrays to store parameter values for each simulation
declare -A thicknessValues
declare -A musamosatValues
declare -A simCounter

# Count simulations per process to track which parameter set belongs to which simulation
while IFS= read -r processLine; do
    processId=$(getProcessId "$processLine")
    if [[ -z "${simCounter[$processId]}" ]]; then
        simCounter[$processId]=0
    else
        ((simCounter[$processId]++))
    fi
done < "$processLines"

# Reset counters for actual processing
declare -A currentSim

# Process thickness values (relTh0) - these appear sequentially before simulations
thicknessIndex=0
while IFS= read -r thicknessLine; do
    thickness=$(getThickness "$thicknessLine")
    thicknessValues["$thicknessIndex"]="$thickness"
    ((thicknessIndex++))
done < "$thicknessLines"

# Reset simulation counters for musamosat processing
declare -A currentSimMu

# Process musamosat values - these appear before each set of simulations
# We need to map them correctly to the simulation order
musamosatIndex=0
while IFS= read -r musamosatLine; do
    musamosat=$(getMusamosat "$musamosatLine")
    musamosatValues["$musamosatIndex"]="$musamosat"
    ((musamosatIndex++))
done < "$musamosatLines"

# Counter for matched pairs
counter=1
musamosatSimIndex=0
thicknessSimIndex=0

echo "Matched Simulations with Parameters and Timing:"
echo "=============================================="
printf "%-3s %-30s %-8s %-12s %-12s %-10s\n" "No." "Simulation Name" "PID" "μ_sat" "Thickness" "CPU Time"
printf "%-3s %-30s %-8s %-12s %-12s %-10s\n" "---" "------------------------------" "--------" "------------" "------------" "----------"

# Read process lines and find matching data
while IFS= read -r processLine; do
    processId=$(getProcessId "$processLine")
    simName=$(getSimName "$processLine")
    processLineNum=$(getLineNumber "$processLine")

    # Find matching timing line
    matchingTiming=$(awk -F: -v pid="$processId" -v startLine="$processLineNum" '
        $1 > startLine && $2 ~ "^"pid"," && /cpu time DFA/ {
            print $0
            exit
        }
    ' "$timingLines")

    if [ -n "$matchingTiming" ]; then
        cpuTime=$(getCpuTime "$matchingTiming")

        # Get thickness for this simulation (sequential mapping)
        thickness="${thicknessValues[$thicknessSimIndex]}"

        # Get musamosat for this simulation (they cycle through the parameter sets)
        musamosat="${musamosatValues[$musamosatSimIndex]}"

        printf "%-3d %-30s %-8s %12.3f %12.3f %-10ss\n" "$counter" "$simName" "$processId" "$musamosat" "$thickness" "$cpuTime"

        ((counter++))
        ((musamosatSimIndex++))
        ((thicknessSimIndex++))

        # Reset musamosat index if we've used all parameter sets
        if [ "$musamosatSimIndex" -ge "$musamosatIndex" ]; then
            musamosatSimIndex=0
        fi
    fi

done < "$processLines"

echo
echo "Summary:"
echo "========"
echo "Total matched simulations: $((counter-1))"
echo "Parameter sets used:"
echo "  Thickness values: $thicknessIndex"
echo "  μ_sat values: $musamosatIndex"

# Process distribution
echo "Process distribution:"
for pid in $(cut -d',' -f1 "$processLines" | cut -d':' -f2 | sort -u); do
    count=$(grep -c "^[0-9]*:$pid," "$processLines")
    echo "  Process $pid: $count simulations"
done

# Clean up temporary files
rm "$processLines" "$timingLines" "$thicknessLines" "$musamosatLines"

echo
echo "Script completed successfully!"

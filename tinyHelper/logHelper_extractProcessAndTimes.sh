#!/bin/bash

# Script to extract and match "runs as process" lines with their corresponding "cpu time DFA" lines
# Usage: ./extractProcessAndTimes <logfile>

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

echo "Extracting process and timing information from: $logFile"
echo "================================================================"
echo

# Create temporary files to store the extracted data
processLines=$(mktemp)
timingLines=$(mktemp)

# Extract all "runs as process" lines with line numbers
grep -n "runs as process" "$logFile" > "$processLines"

# Extract all "cpu time DFA" lines with line numbers
grep -n "cpu time DFA" "$logFile" > "$timingLines"

# Function to extract process ID from a line
getProcessId() {
    echo "$1" | cut -d',' -f1 | cut -d':' -f2
}

# Function to extract simulation name from process line
getSimName() {
    echo "$1" | sed 's/.*INFO - \(.*\) runs as process.*/\1/'
}

# Function to extract CPU time from timing line
getCpuTime() {
    echo "$1" | sed 's/.*cpu time DFA = \([0-9.]*\) s.*/\1/'
}

# Counter for matched pairs
counter=1

echo "Matched Process and Timing Pairs:"
echo "=================================="

# Read process lines and find matching timing lines
while IFS= read -r processLine; do
    processId=$(getProcessId "$processLine")
    simName=$(getSimName "$processLine")
    processLineNum=$(echo "$processLine" | cut -d':' -f1)

    # Find the next timing line with the same process ID that comes after this process line
    matchingTiming=$(awk -F: -v pid="$processId" -v startLine="$processLineNum" '
        $1 > startLine && $2 ~ "^"pid"," && /cpu time DFA/ {
            print $0
            exit
        }
    ' "$timingLines")

    if [ -n "$matchingTiming" ]; then
        cpuTime=$(getCpuTime "$matchingTiming")
        printf "%2d. %-30s (PID: %s) → CPU time: %ss\n" "$counter" "$simName" "$processId" "$cpuTime"
        ((counter++))
    fi

done < "$processLines"

echo
echo "Summary:"
echo "========"
echo "Total matched pairs: $((counter-1))"

# Count simulations per process ID
echo "Process distribution:"
for pid in $(cut -d',' -f1 "$processLines" | cut -d':' -f2 | sort -u); do
    count=$(grep -c "^[0-9]*:$pid," "$processLines")
    echo "  Process $pid: $count simulations"
done

# Clean up temporary files
rm "$processLines" "$timingLines"

echo
echo "Script completed successfully!"

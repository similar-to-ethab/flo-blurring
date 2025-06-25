#!/bin/bash

# Default values
CFLAGS_VALUE=""
USE_NOHUP=false

# Parse arguments
for arg in "$@"; do
    if [ "$arg" == "bview" ]; then
        CFLAGS_VALUE="-DDISPLAY=-1"
    elif [ "$arg" == "bg" ]; then
        USE_NOHUP=true
    fi
done

# Build command
BUILD_CMD="make blur.tst"
[ -n "$CFLAGS_VALUE" ] && BUILD_CMD="CFLAGS=$CFLAGS_VALUE $BUILD_CMD"

# Run the build
if $USE_NOHUP; then
    echo "Running in background with nohup..."
    nohup bash -c "$BUILD_CMD" > blur.log 2>&1 &
    disown
    echo "Output is being logged to blur.log"
else
    echo "Running in foreground..."
    eval "$BUILD_CMD"
fi


#!/bin/zsh

# Remote server info
USER="ethanps"
HOST="comech-1172.me.sc.edu"
REMOTEPATH="/home/ethanps/basilisk/src/test/blur/blur"
PATTERN="video-*-*-*.mp4"
DEST="$BLUR/group_mp4/newton_mp4"

# Secure copy with wildcard pattern (quoted to avoid zsh expansion)
scp "${USER}@${HOST}:${REMOTEPATH}/${PATTERN}" "$DEST"

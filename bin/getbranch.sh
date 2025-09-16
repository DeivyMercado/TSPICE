#!/bin/bash
#Gets the current Git branch name (removes spaces)
branch=$(git branch | grep "*" | cut -f 2 -d '*')
echo "${branch//[[:space:]]/}"
#!/bin/sh
# ==============================================================================
#title        :my_script.sh
#description  :Description of ouput file(s)/stream
#author       :Forename Surname
#date         :XX/XX/XXXX
#version      :0.1
#usage        :sh SKEL.sh file1 [file2] [output_tag]
#input        :file1 (required)
#               - description
#              file2 (optional)
#               - description
#               - default: default.file
#notes        :
# ================================================================================

# this shell will kill itself if it takes over __ kb of virtual memory
# note: 8388608 kb = 8 gb
# ulimit -v 8388608


# --------------------------- Get and check user input -------------------------
file1=$1
file2=${2-"/path/to/default.file"}

# this will require file1 to be present
if [ -z $file1 ]; then
	echo "ERROR: file1 required"
	exit
fi
# ------------------------------------------------------------------------------


# --------------------------- Functions  ---------------------------------------
# function1() {
	# example function
# 	echo "Function 1"
# }
# ------------------------------------------------------------------------------


# --------------------------- Main ---------------------------------------------
# example if clause
# if [ "$file1" != 'comparison_string' ]; then
# 	something
# else
# 	something
# fi
# ------------------------------------------------------------------------------

#! /bin/sh

# Linux data-gathering commands; adjust as necessary for your platform.
#
# Be sure to remove any information from the output that would violate
# SC's double-blind review policies.

env | sed "s/$USER/USER/g"
set -x
lsb_release -a
uname -a
lscpu || cat /proc/cpuinfo
cat /proc/meminfo
inxi -F -c0
lsblk -a
lsscsi -s
module list
nvidia-smi
(lshw -short -quiet -sanitize || lspci) | cat

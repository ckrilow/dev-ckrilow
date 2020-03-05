#!/usr/bin/env python3

# submits sge jobs

import os
import sys
import math

from snakemake.utils import read_job_properties

# below is a snakemake way to do this
jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)


cluster = job_properties.get('cluster', {})

group = cluster.get('group', 'snk') # job group
queue = cluster.get('queue', None) # job queue
job_name = cluster.get('name', 'snk') # job name

threads = job_properties.get('threads', 1)
memory = cluster.get('memory', 10) # default to 10Gb
if memory:
    # mem is specified as total memory in Gb, but we
    # need to give it to SGE as memory per thread in Mb
    memory = int(math.floor(float(memory * 1e9) / (threads * 1e6)))
resources = 'mem_free={0}M,h_vmem={0}M'.format(
    memory, memory
)

output = cluster.get('output', None)
if output:
    output = os.path.realpath(output)
    tdir = os.path.dirname(output)
    if not os.path.exists(tdir): os.makedirs(tdir)
error = cluster.get('error', None)
if error:
    error = os.path.realpath(error)
    tdir = os.path.dirname(error)
    if not os.path.exists(tdir): os.makedirs(tdir)


# build the command
cmd = 'qsub'
if queue:
    cmd = '{} -q {}'.format(cmd, queue)
cmd = '{} -N "{}" -V -cwd -l "{}"'.format(
    cmd,
    job_name,
    resources
)
if threads > 1:
    cmd = '{} -pe make-dedicated {}'.format(cmd, threads)
if output:
    cmd = '{} -o {}'.format(cmd, output)
if error:
    cmd = '{} -e {}'.format(cmd, error)
cmd = '{} {}'.format(cmd, jobscript)


# run the command
os.system(cmd)

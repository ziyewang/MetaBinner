#!/usr/bin/env python

import operator, subprocess, time, math

############################################################
# util
#
#
############################################################

############################################################
# exec_par
#
# Execute the commands in the list 'cmds' in parallel, but
# only running 'max_proc' at a time.
############################################################
def exec_par(cmds, max_proc):
    total = len(cmds)
    finished = 0
    running = 0
    p = []

    while finished + running < total:
        # launch jobs up to max
        while running < max_proc and finished+running < total:
            p.append(subprocess.Popen(cmds[finished+running], shell=True))
            #print 'Running %d' % p[running].pid
            running += 1

        # are any jobs finished
        new_p = []
        for i in range(len(p)):
            if p[i].poll() != None:
                running -= 1
                finished += 1
            else:
                new_p.append(p[i])

        # if none finished, sleep
        if len(new_p) == len(p):
            time.sleep(1)
        p = new_p

    # wait for all to finish
    for i in range(len(p)):
        p[i].wait()

############################################################
# max_i
#
# Find max and return index and value
############################################################
def max_i(lis):
    max_index = 0
    max_val = lis[0]
    for i in range(1,len(lis)):
        if lis[i] > max_val:
            max_index = i
            max_val = lis[i]

    return (max_val,max_index)

############################################################
# min_i
#
# Find min and return index and value
############################################################
def min_i(lis):
    min_index = 0
    min_val = lis[0]
    for i in range(1,len(lis)):
        if lis[i] < min_val:
            min_index = i
            min_val = lis[i]

    return (min_val,min_index)

############################################################
# mean
#
# Return mean of a list
############################################################
def mean(ls):    
    return float(sum(ls)) / float(len(ls))

############################################################
# sd
#
# Return the standard deviation of a list
############################################################
def sd(ls):
    u = mean(ls)
    dev_sum = 0.0
    for x in ls:
        dev_sum += (x-u)*(x-u)
    return math.sqrt(dev_sum / float(len(ls)))

############################################################
# sort_dict
#
# Sort a dict by the values, returning a list of tuples
############################################################
def sort_dict(hash, reverse=False):
    return sorted(hash.items(), key=operator.itemgetter(1), reverse=reverse)

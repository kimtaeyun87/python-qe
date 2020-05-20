#-*- coding: utf-8 -*-


import subprocess as subp


def grep(key, fpath, *args):
    cmd = ['grep', key, fpath] + [arg for arg in args]
    with subp.Popen(cmd, stdout=subp.PIPE).decode('utf-8') as p:
        result = p.stdout.readlines()
    return result



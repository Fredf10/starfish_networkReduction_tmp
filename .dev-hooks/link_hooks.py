#!/usr/bin/python
#
# Symlinks the hooks to the local git directory
ret_val = 0 
from __future__ import print_function, absolute_import
import os
import stat
import errno
import sys
import git

from subprocess import check_output, CalledProcessError

def root():
    ''' returns the absolute path of the repository root '''
    try:
        base = check_output('git rev-parse --show-toplevel', shell=True)
    except CalledProcessError:
        raise IOError('Current working directory is not a git repository')
    return base.decode('utf-8').strip()


repo = git.Repo(root())
local_hooks_path = os.path.abspath(os.path.split(__file__)[0])
filenames = next(os.walk(local_hooks_path))[2]
git_hooks_path = repo.git_dir+"/hooks"

valid_hooks = ["pre-commit", "post-commit", "prepare-commit-msg", "commit-msg",
               "pre-rebase", "post-rewrite", "post-checkout", "post-merge",
               "pre-push", "post-receive"]
for filename in filenames:
    hook_name = filename.split(".py")[0]
    if hook_name in valid_hooks:
        target_path = os.path.join(local_hooks_path, filename)
       
        st = os.stat(target_path)
        os.chmod(target_path, st.st_mode | stat.S_IEXEC)
        st = os.stat(target_path)
        # TODO: Python ize this ?
        bash_exec = "if [ -x "+  target_path + " ]; then true else false; fi"
        exec_stat = os.system(bash_exec)
        if exec_stat != 0:
            print("Unable to set hook as an executable" )
            ret_val = 1

        source_path = os.path.join(git_hooks_path, hook_name)
        try:
            os.symlink(target_path, source_path)
        except OSError as e:
            if e.errno != errno.EEXIST:
               ret_val = 1
               raise
            else:
                current_target = os.path.realpath(os.readlink(source_path))
                linkIsGood = os.path.samefile(current_target,target_path)

                if linkIsGood is not True:
                    print(hook_name, "hook already exists in ", git_hooks_path)
                    ret_val = 1
sys.exit(ret_val)


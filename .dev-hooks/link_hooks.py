#!/usr/bin/python
#
# Symlinks the hooks to the local git directory

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
        source_path = os.path.join(git_hooks_path, hook_name)
        target_path = os.path.join(local_hooks_path, filename)
        try:
            os.symlink(target_path, source_path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
            else:
                print hook_name, "hook already exists in ", git_hooks_path



# -*- coding: utf-8 -*-
# Calculates the current version number.  If possible, this is the
# output of “git describe”, modified to conform to the versioning
# scheme that setuptools uses.  If “git describe” returns an error
# (most likely because we're in an unpacked copy of a release tarball,
# rather than in a git working copy), then we fall back on reading the
# contents of the RELEASE-VERSION file.
#
# To use this script, simply import it your setup.py file, and use the
# results of getGitVersion() as your package version:
#
# from version import *
#
# setup(
#     version=getGitVersion(),
#     .
#     .
#     .
# )
#
#
# This will automatically update the RELEASE-VERSION file, if
# necessary.  Note that the RELEASE-VERSION file should *not* be
# checked into git; please add it to your top-level .gitignore file.
#
# You'll probably want to distribute the RELEASE-VERSION file in your
# sdist tarballs; to do this, just create a MANIFEST.in file that
# contains the following line:
#
#   include RELEASE-VERSION

from subprocess import Popen, PIPE
import os
import subprocess

def callGitDescribe():
    """Determine a version according to git.

    This calls `git describe`, if git is available.

    Returns:
        version from `git describe --tags --always --dirty` if git is
        available; otherwise, None
    """
    try:
        #  os.chdir(get_project_path())
        cmd = ['git', 'describe', '--tags', '--always', '--dirty']
        out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL)
        if not isinstance(out, str):
            out = out.decode('utf-8')
        ver = out.strip()
    except Exception:
        ver = None
    return ver

def readReleaseVersion():
    '''Read version from text file RELEASE-VERSION'''
    try:
        f = open("RELEASE-VERSION", "r")

        try:
            version = f.readlines()[0]
            return version.strip()

        finally:
            f.close()

    except:
        return None

def writeReleaseVersion(version):
    '''write version to text file RELEASE-VERSION'''
    f = open("RELEASE-VERSION", "w")
    f.write("%s\n" % version)
    f.close()

def getGitVersion():

    # Read in the version that's currently in RELEASE-VERSION.
    releaseVersion = readReleaseVersion()

    # First try to get the current version using “git describe”.
    version = callGitDescribe()

    # If that doesn't work, fall back on the value that's in
    # RELEASE-VERSION.
    if version is None:
        version = releaseVersion
        
    # If we still don't have anything, that's an error.
    if version is None:
        raise ValueError("Cannot find the version number!")

    # If the current version is different from what's in the
    # RELEASE-VERSION file, update the file to be current.
    if version != releaseVersion:
        writeReleaseVersion(version)
#
    # Finally, return the current version.
    return version


if __name__ == "__main__":
    print (getGitVersion())

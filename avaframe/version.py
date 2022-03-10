# -*- coding: utf-8 -*-
# Calculates the current version number.  If possible, this is the
# output of “git describe”, modified to conform to the versioning
# scheme that setuptools uses.  If “git describe” returns an error
# (most likely because we're in an unpacked copy of a release tarball,
# rather than in a git working copy), then we fall back on reading the
# contents of the RELEASE-VERSION file.
#
# To use this script, simply import it your setup.py file, and use the
# results of getVersion() as your package version:
#
# from version import *
#
# setup(
#     version=getVersion(),
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

from subprocess import check_output, run, CalledProcessError, DEVNULL, Popen, PIPE
import os
import subprocess
import sys


def getProjectPath():
    """Determine absolute path to the top-level of the catch project.

    This is assumed to be the parent of the directory containing this script.

        Returns
        -------
        path : string
            to top-level of the avaframe project
    """
    # abspath converts relative to absolute path; expanduser interprets ~
    path = __file__  # path to this script
    path = os.path.expanduser(path)  # interpret ~
    path = os.path.abspath(path)  # convert to absolute path
    path = os.path.dirname(path)  # containing directory
    return path


def getVersionFromGitDescribe():
    """Determine a version according to git.

    This calls `git describe`, if git is available.

        Returns
        -------
        version from `git describe --tags --always --dirty` if git is
        available; otherwise, None
    """
    cwd = os.getcwd()
    try:
        os.chdir(getProjectPath())
        cmd = ['git', 'describe', '--tags', '--always', '--dirty']
        out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL)
        if not isinstance(out, str):
            out = out.decode('utf-8')
        ver = out.strip()
        # Pep 440 compliance
        ver = ver.replace('_', '')
        ver = ver.replace('v', '')
        ver = ver.replace('-', '+', 1)
        ver = ver.replace('-', '.')
    except Exception:
        ver = None
    os.chdir(cwd)
    return ver


def releaseFile():
    """Obtain path to file storing version, according to git.

        Returns
        -------
        path to RELEASE-VERSION file
    """
    return os.path.join(getProjectPath(), 'RELEASE-VERSION')


def readReleaseVersion():
    """Read RELEASE-VERSION file, containing git version.

        Returns
        -------
        if RELEASE-VERSION file exists, version stored in it; otherwise, None
    """
    try:
        with open(releaseFile(), 'rt') as inf:
            version = inf.readlines()[0].strip()
    except Exception:
        version = None
    return version


def writeReleaseVersion(version):
    """Save version, according to git, into VERSION file.

        Parameters
        ----------
        version:  str
            version to save
    """
    with open(releaseFile(), 'wt') as outf:
        outf.write(version + '\n')


def getVersion():
    """ Calculates the current version number.  If possible, this is the
        output of “git describe”, modified to conform to the versioning
        scheme that setuptools uses.  If “git describe” returns an error
        (most likely because we're in an unpacked copy of a release tarball,
        rather than in a git working copy), then we fall back on reading the
        contents of the RELEASE-VERSION file.

        Returns
        -------
        version : string
    """

    # Read in the version that's currently in RELEASE-VERSION.
    releaseVersion = readReleaseVersion()

    # First try to get the current version using “git describe”.
    version = getVersionFromGitDescribe()

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
    print(getVersion())

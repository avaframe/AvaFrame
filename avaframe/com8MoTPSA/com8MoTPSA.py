import os
import subprocess
import platform


def _runAndCheck(command):
    if os.name == "nt":
        useShell = True
    elif platform.system() == "Darwin":
        useShell = False
    else:
        useShell = False

    # This starts the subprocess
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=useShell,
        encoding="utf-8",
        errors="replace",
        universal_newlines=True,
    )

    while True:
        realtimeOutput = process.stdout.readline()

        if realtimeOutput == "" and process.poll() is not None:
            break

        if realtimeOutput:
            line = realtimeOutput.strip()
            print(line)


# Example usage

def cfgToRcf(cfg, fileName):
    with open(fileName, 'w') as f:
        for section in cfg.sections():
            if section in ("FOREST_EFFECTS", "ENTRAINMENT"):
                pass
            elif section in ("GENERAL"):
                continue
            else:
                f.write(f"# {section.replace('_', ' ')}\n")
                f.write("#\n")
            for key, value in cfg.items(section):
                # key = key.replace('_', ' ')
                key = key.strip()
                f.write(f"{key:<40}{value}\n")
            f.write("#\n")


def com8MoTPSAMain(cfgMain, cfgInfo=None):
    # Generate command and run via subprocess.run

    rcfFile = 'outputTemplate.rcf'

    cfgToRcf(cfgInfo, rcfFile)

    command = ['./MoT-PSA', rcfFile]
    _runAndCheck(command)

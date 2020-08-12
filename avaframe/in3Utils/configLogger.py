""" Defining a writable object to write the config to the log file"""
class ConfigLogger(object):
    def __init__(self, log):
        self.__log = log
    def __call__(self, config):
        self.__log.info("Config:")
        config.write(self)
    def write(self, data):
        # stripping the data makes the output nicer and avoids empty lines
        line = data.strip()
        self.__log.info(line)

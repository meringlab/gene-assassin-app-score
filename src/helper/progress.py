import logging
import timeit


class ProgressLogger(object):
    def __init__(self, report):
        self.counter = 0
        self.start = timeit.default_timer()
        self.report = report

    def log(self):
        self.counter += 1
        if self.counter % self.report == 0:
            stop = timeit.default_timer()
            logging.info('%d processed in %dsec', self.counter, stop - self.start)

class Region(object):
    def __init__(self, chromosome, start, end):
        start = int(start)
        end = int(end)
        if start > end:
            raise Exception('invalid location: %d-%d' % (start, end))
        self.start = start
        self.end = end
        self.chromosome = chromosome

    @classmethod
    def parse(cls, str):
        chr = str[:str.index(':')]
        # start, stop = list(map(int, str[str.index(':') + 1:].split('-')))
        start, end = str[str.index(':') + 1:].split('-')
        return cls(chr, start, end)

    def __str__(self):
        return '%s:%d-%d' % (self.chromosome, self.start, self.end)

    def __repr__(self):
        return str(self)

class Dummy:
    def __init__(self):
        self.data = []

    def append(self, value):
        self.data.append(value)

    def __getitem__(self, values):
        return self.data[values]


if __name__ == '__main__':
    dum = Dummy()
    for i in range(100):
        dum.append(i)
    print(dum[1:3])

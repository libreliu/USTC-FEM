x = []
for v in range(2):
    def pv(v):
        print(v)
    x.append(lambda : pv(v))

for xx in x:
    xx()

x = []
for v in range(2):
    def pv(v):
        print(v)
    x.append((lambda v: lambda : pv(v))(v))

for xx in x:
    xx()
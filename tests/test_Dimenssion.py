from pytraj.Dimension import Dimension

d = Dimension()
d.Label = "test"
d.Min = 10
d.Max = 20
d.Bins = 100

print("test instance d")
print(d.MinIsSet())
print(d.MaxIsSet())
d.PrintDim()

print("test instance d2")
d2 = Dimension(d)
print(d2.MinIsSet())
print(d2.MaxIsSet())
print("d and d2 are the same instance? %s " % (d == d2))

print("create instance d3")
d3 = Dimension(100, 10, 10000, "d3")

# test wrong inputs
# d3 = Dimension("test", 100, 10, 10000)

print(d3.Min)
print(d3.Max)
print(d3.Bins)
print(d3.Step)
# print d3.CalcBinsOrStep()
print(d3.Coord(100))

# 
d4 = d3
print("d4 and d3 are the same instance? %s" % (d4 == d3))

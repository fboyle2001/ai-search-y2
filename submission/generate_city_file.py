import random

n = 96
lines = [f"SIZE = {n},"]

upper = 200
lower = 10

for i in range(n, 0, -1):
    line = ""

    for j in range(i):
        if line != "":
            line += ","

        line += str(random.randint(lower, upper))

    lines.append(line)

lines.append("NOTE =")

str_n = str(n)

while len(str_n) < 3:
    str_n = "0" + str_n

with open(f"AISearchfile{str_n}.custom", "w") as f:
    f.write("\n".join(lines))

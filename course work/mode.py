
X = [ [1.0698, 0.94789],[1.0694, 0.94961],[0.96973, 3.7748],[1.0839, 0.90023],
[1.0754, 0.92911]]
# X = [[1, 4], [5, 9], [1.5, 4.5], [6, 9]]
data = []
for a in X:
    for b in a:
        data.append(b)
data = sorted(list(set(data)))
print("X:", data, "\n")


z = []
for i in range(len(data) - 1):
    z.append([data[i], data[i + 1]])
print("Z:", z, "\n")
index = []
for i in range(len(z)):
    count = 0
    for j in range(len(X)):
        if z[i][0] >= X[j][0] and z[i][1] <= X[j][1]:
            count += 1
    index.append(count)
print("$\mu_k$:", index)

import pandas as pd

df = pd.read_csv('~/Desktop/hg38_random_10000_guides.csv')
with open('/home/schmidt73/Desktop/hg38_all_sgrnas.txt', 'r') as f:
    grnas = []
    for i, line in enumerate(f):
        if i > 1000000: break
        grnas.append(line.rstrip()[:10])

def make_trie(*words):
    root = dict()
    for word in words:
        current_dict = root
        for letter in word:
            current_dict = current_dict.setdefault(letter, {})
    return root

def size_trie(t):
    if not t:
        return 0

    count = 0
    for (k, v) in t.items():
        count += 1 + size_trie(v)

    return count

trie = make_trie(*grnas)
print(size_trie(trie) / (len(grnas) * 10))

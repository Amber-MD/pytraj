# adapted
# from http://stackoverflow.com/questions/16075966/find-the-count-of-the-words-in-the-text-file-using-python
import sys
from collections import Counter
from string import punctuation

counter = Counter()
with open(sys.argv[1]) as f:
    for line in f:
        counter.update(word.strip(punctuation) for word in line.split())

result = dict(counter)

for key, value in result.items():
    if value > 1 and key and ".h" not in key:
        print(key, value)

import re

with open('test.txt', 'r') as fh:
    text = ''
    for line in fh:
        text = text + line
    string = r'Accession: [A-Z][A-Z]_[A-Z0-9]+\.\d\. None of the genes of interest were found in the \'genes_fasta\': \[\'dnaA\', \'dnaN\'\]\n        Will not use gene locations in prediction.\n Gene dict contained 3'
    matches = re.findall(string, text)
accs = []
for match in matches:
    acc = re.search(r'[A-Z][A-Z]_[A-Z0-9]+', match)
    accs.append(acc.group())
print(accs)

def readfile(filename):
    result = []
    with open(filename) as infile:
        for line in infile:
            result.append(line)
    
    return result
        
import pickle

database_filename = 'examplePickle'

words = ["One", "Two"]

def storeData():
    words = ["Three", "Four"]
    
    # ab: using binary mode
    dbfile = open(database_filename, 'ab')

    pickle.dump(words, dbfile)
    dbfile.close()

def loadData():
    # rb: read binary
    dbfile = open(database_filename, 'rb')

    words = pickle.load(dbfile)
    print(words)
    dbfile.close()

    return words


storeData()
words = loadData()
print(words)
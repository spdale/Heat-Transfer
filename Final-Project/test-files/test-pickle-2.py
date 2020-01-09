import pickle
import os

database_filename = 'examplePickle.pickle'

nums = [1, 2]
words = ["Apple", "Banana"]

def removePickle():
    if os.path.exists(database_filename):
        os.remove(database_filename)

def storeData():
    removePickle()
    nums = [3, 4]
    words = ["Peach", "Grape"]

    store = [nums, words]
    
    # ab: using binary mode
    dbfile = open(database_filename, 'ab')

    pickle.dump(store, dbfile)
    dbfile.close()

def loadData():
    # rb: read binary
    dbfile = open(database_filename, 'rb')

    read = pickle.load(dbfile)
    # print(words)
    dbfile.close()

    return read


storeData()

output = loadData()
nums = output[0]
words = output[1]
print(output)

print(nums)
print(words)
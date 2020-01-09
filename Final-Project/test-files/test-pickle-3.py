import pickle
import os

database_filename = "test-data.pickle"

def remove_pickle(filename):
    if os.path.exists(filename):
        os.remove(filename)

def store_data(filename, data):
    remove_pickle(filename)
    # ab: using binary mode
    dbfile = open(filename, 'ab')
    pickle.dump(data, dbfile)
    dbfile.close()

def load_data(filename):
    if os.path.exists(filename):
        # rb: read binary
        dbfile = open(filename, 'rb')
        read = pickle.load(dbfile)
        dbfile.close()

        return read
    return False


# nums = [12, 22]
# store_data(database_filename, nums)

nums = load_data(database_filename)
print(nums)
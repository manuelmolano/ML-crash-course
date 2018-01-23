#askdjhaskjdhaksjhdkjashd

from data_analysis import low_pass

def function_to_load_data(filename):
    # simple reading of binary
    data = np.fromfile(filename, dtype='float32')
    # but low pass filtering
    new_data = low_pass(data)
    return new_data

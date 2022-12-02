import numpy as np

class Lumping():
    def __init__(self):
        ...

    def lumy(self, lim, data, col_name):

        '''Function to insert ID index to lumpy function

        Parameters
        ----------
        lim: list of tuples
            ranges limits
        data: set
            column to realize lumping

        Return
        ------
        lumy_set: set
            data set with ID refers to lumping
        '''

        count = -1
        co = np.array([])

        for i in data[col_name]:
            count += 1
            if i in range(lim[count][0],
                          lim[count][1]+1,
                          1):
                co = np.append(co, count)
                count -=1
            else:
                co = np.append(co, count+1)

        data['ID'] = co
        data.set_index([pd.Index(co), 'ID'])

        return data

    def lumpy(self, lim, data, col_name):

        '''Function to determinate carbon ranges
        for lumping

        Parameters
        ----------
        lim: tupla
            ranges limits
        data: set
            distribution data set
        col_name: set
            column to realize lumping

        Return
        ------
        lumping_set: set
            lumping data set (groupy)
        '''

        lp = Lumping()
        df = lp.lumy(lim, data, col_name)

        data_lumping = df.groupby(by=['ID'],
                                  dropna=False).mean()
        data_lumping[col_name] = lim

        return data_lumping